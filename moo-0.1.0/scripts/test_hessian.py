# SPDX-License-Identifier: LGPL-3.0-or-later
#
# This file is part of MOO - Modelica / Model Optimizer
# Copyright (C) 2025 University of Applied Sciences and Arts
# Bielefeld, Faculty of Engineering and Mathematics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

F = [[0, 2], [0, 2], [1], [1, 2], [3]]

# Step 1: Build sparsity map: M[(v1, v2)] → list of function outputs (rows)
M = dict()
for f in range(len(F)):
    for v1 in F[f]:
        for v2 in F[f]:
            key = tuple(sorted((v1, v2)))
            if key not in M:
                M[key] = [f]
            else:
                M[key].append(f)

# Step 2: Define colors and reverse lookup
COLORS = [[0, 1], [2, 3]]

# Step 3: Build Q[(c1, c2)] → [ ([i1, i2, ...], nnz_index) ]
Q = dict()
nnz_lookup = dict()
NNZ = 0

for i in range(5):
    for j in range(i + 1):
        key = (j, i) if j <= i else (i, j)
        if key in M:
            print(key)
            nnz_lookup[key] = NNZ
            NNZ += 1

for c1 in range(len(COLORS)):
    for c2 in range(c1 + 1):  # symmetric
        pair_entries = []
        for i1 in COLORS[c1]:
            for i2 in COLORS[c2]:
                key = tuple(sorted((i1, i2)))
                if key in M:
                    row_indices = tuple(set(M[key]))
                    flat_idx = nnz_lookup[key]
                    pair_entries.append([row_indices, flat_idx])
        Q[(c1, c2)] = pair_entries

print("Q mapping (color pairs → [([rows], flat_idx)]):\n")
for key, val in Q.items():
    print(f"{key}: {val}")

print(f"\nTotal NNZ (Lower) entries: {NNZ}")

import numpy as np

def f(x):
    x0, x1, x2, x3 = x
    return np.array([
        x0**2 + x2**2,
        2*x0**2 + 3*x2**2,
        2*x1**2,
        x1**4 + x2**4 * x1**2,
        x3**3
    ])
def dense_hessian_f_weighted(x, lambd):
    x0, x1, x2, x3 = x
    H = np.zeros((4, 4))

    # f0 = x0^2 + x2^2
    H[0, 0] += 2 * lambd[0]
    H[2, 2] += 2 * lambd[0]

    # f1 = 2*x0^2 + 3*x2^2
    H[0, 0] += 4 * lambd[1]
    H[2, 2] += 6 * lambd[1]

    # f2 = 2*x1^2
    H[1, 1] += 4 * lambd[2]

    # f3 = x1^4 + x2^4 * x1^2
    # ∂²f3/∂x1² = 12*x1^2 + 2*x2^4
    H[1, 1] += (12 * x1**2 + 2 * x2**4) * lambd[3]

    # ∂²f3/∂x2² = 12*x2^2 * x1^2
    H[2, 2] += (12 * x2**2 * x1**2) * lambd[3]

    # ∂²f3/∂x1∂x2 = 8 * x2^3 * x1
    mix = 8 * x2**3 * x1 * lambd[3]
    H[1, 2] += mix
    H[2, 1] += mix  # symmetric

    # f4 = x3^3
    # ∂²f4/∂x3² = 6 * x3
    H[3, 3] += 6 * x3 * lambd[4]

    return H

def jvp(x, s):
    x0, x1, x2, x3 = x
    s0, s1, s2, s3 = s

    return np.array([
        2 * x0 * s0 + 2 * x2 * s2,
        4 * x0 * s0 + 6 * x2 * s2,
        4 * x1 * s1,
        4 * x1**3 * s1 + (4 * x2**3 * x1**2) * s2 + (2 * x2**4 * x1) * s1,
        3 * x3**2 * s3
    ])


def seed_vector(color):
    s = np.zeros(4)
    for i in color:
        s[i] = 1.0
    return s

def compute_hessian_entries(x, lambd, h=1e-4):
    NNZ = max(flat_idx for pairs in Q.values() for _, flat_idx in pairs) + 1
    Hflat = np.zeros(NNZ)

    # Precompute base JVP for each color seed
    J = {}
    for c, color in enumerate(COLORS):
        s = seed_vector(color)
        J[c] = jvp(x, s)

    # Loop over color pairs to do perturbed evaluations and finite differences
    for c1 in range(len(COLORS)):
        s1 = seed_vector(COLORS[c1])
        x1 = x + h * s1  # perturb x in direction of color c1

        for c2 in range(c1 + 1):
            s2 = seed_vector(COLORS[c2])
            J_pert = jvp(x1, s2)

            # For each nnz entry associated with color pair (c1,c2)
            for rows, nnz_index in Q[(c1, c2)]:
                out = 0
                for r in rows:
                    Hflat[nnz_index] += lambd[r] * (J_pert[r] -  J[c2][r]) / h

    return Hflat


def richardson_extrapolation(f, x, lambd, h, n=2, order=1, base=10):
    """
    Perform Richardson extrapolation for a function that returns arrays.
    
    Parameters:
    - f: function to extrapolate (returns np.array)
    - x: point at which to evaluate the function
    - h: initial step size
    - n: number of extrapolation steps
    
    Returns:
    - extrapolated array
    - extrapolation table (3D array: steps x steps x array_dimensions)
    """
    # First evaluation to determine array shape
    sample_output = f(x, lambd, h)
    array_shape = sample_output.shape
    
    # Initialize the extrapolation table as a 3D array
    table = np.zeros((n, n, *array_shape))
    
    # Fill first column with function evaluations at different h values
    for i in range(n):
        table[i, 0] = f(x, lambd, h / (base**i))
    
    # Perform Richardson extrapolation (array operations)
    for j in range(1, n):
        factor = base**(order * j)
        for i in range(j, n):
            table[i, j] = (factor * table[i, j-1] - table[i-1, j-1]) / (factor - 1)
    
    return table[n-1, n-1], table


x = np.array([1, 2, 4, 5])
lambd = np.array([1, 2, 5, 1, -5])

"""
A, B = richardson_extrapolation(compute_hessian_entries, x, lambd, 1e-5, n=2, base=2)
print(A)

A, B = richardson_extrapolation(compute_hessian_entries, x, lambd, 1e-3, n=3, base=2)
print(A)
"""
H = compute_hessian_entries(x, lambd, 1e-5)
print(dense_hessian_f_weighted(x, lambd))
print(H)

# grid serach
import numpy as np

h_values = [1, 0.5, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9]
n_values = [1, 2]
base_values = [1.25, 1.5, 1.75, 2, 3, 5, 10, 20, 50, 100]
exact = np.array([10, 580, 0, 1024, 782, -150])
results = []
for h in h_values:
    for n in n_values:
        for base in base_values:
            A, _ = richardson_extrapolation(compute_hessian_entries, x, lambd, h, n=n, order=1, base=base)
            error = np.linalg.norm(A - exact, np.inf)
            results.append({
                'h': h,
                'n': n,
                'base': base,
                'error': error,
                'result': A
            })


results.sort(key=lambda x: x['error'])

print("Top configurations by error:")
for i, res in enumerate(results):
    print(f"{i+1}. h={res['h']:.0e}, n={res['n']}, base={res['base']}, error={res['error']:.10f}")
    print("   Result: [" + ", ".join(f"{x:.10f}" for x in res['result']) + "]")
    print()

## 1-norm                                      | 2-norm                                      | oo-norm
#  n=4: h=1e+00, base=2,    error=0.0000000000 | n=4: h=1e+00, base=2,    error=0.0000000000 | n=4: h=1e+00, base=2,  error=0.0000000000
#  n=3: h=1e-03, base=1.75, error=0.0000000028 | n=3: h=1e-03, base=1.75, error=0.0000000018 | n=3: h=1e-02, base=10, error=0.0000000012
#  n=2: h=5e-05, base=1.75, error=0.0000000998 | n=2: h=5e-05, base=1.75, error=0.0000000775 | n=2: h=5e-05, base=10, error=0.0000000731
#  n=1: h=5e-08, base=x,    error=0.0000198226 | n=1: h=5e-08, base=x,    error=0.0000131057 | n=1: h=1e-08, base=x,  error=0.0000106304

# Re-execute due to state reset

import numpy as np
from sympy import symbols, exp, diff, lambdify, Matrix, sin, cos

# Define symbols
x1, x2, x3, x4, u = symbols('x1 x2 x3 x4 u')

# Constants
R = 1.9872
T = 700 * u

# Reaction rate expressions
k1 = exp(8.86 - 20300 / R / T)
k2 = exp(24.25 - 37400 / R / T)
k3 = exp(23.67 - 33800 / R / T)
k4 = exp(18.75 - 28200 / R / T)
k5 = exp(20.70 - 31000 / R / T)

# Define the function
F = []
dec = 0
if dec == 0:
    F.append((-k1 * x1) - (k3 + k4 + k5) * x1 * x2)
    F.append(k1 * x1 - k2 * x2 + k3 * x1 * x2)
    F.append(k2 * x2 + k4 * x1 * x2)
    F.append(k5 * x1 * x2)
    F.append((-k1 * x1) - (k3 + k4 + k5) * x1 * x2)
    V = [x1, x2, x3, x4, u]
    x1_val = 1.35266244529518659e-01
    x2_val = 3.55066498452544121e-01
    x3_val = 3.34360423673394747e-01
    x4_val = 1.75306833344556989e-01
    u_val = 9.70337348817454814e-01
elif dec == 1:
    F.append(sin(x1 + x2 - 1))
    F.append(cos(x2 + x3 - 2))
    F.append((x3 + u) / (x3**2 + u**2))
    F.append(-1 + x3*x3)
    V = [x1, x2, x3, u]
    x1_val = 1.478975818445809
    x2_val = 1.470723942456369
    x3_val = 1.11596063078672
    u_val = 4.999999970445199

lambd = [0, 0, 0, 0, 1]

#f = sum(lambd[i] * F[i] for i in range(4))
f = F[4]
print(f)
# Variables vector
vars_vec = Matrix(V)

# Compute Hessian
hessian = f.diff(vars_vec).jacobian(vars_vec)

# Lambdify Hessian for numerical evaluation
hessian_func = lambdify((x1, x2, x3, x4, u), hessian, 'numpy')

# Evaluate at given point



np.set_printoptions(precision=5)

hessian_numeric = np.array(hessian_func(x1_val, x2_val, x3_val, x4_val, u_val), dtype=np.float128)
print()
print(hessian_numeric)



#### test simulation
t = [0, 0.08, 0.16, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88, 0.96, 1.04, 1.12, 1.2, 1.28, 1.36, 1.44, 1.52, 1.6, 1.68, 1.76, 1.84, 1.92, 2, 2.08, 2.16, 2.24, 2.32, 2.4, 2.48, 2.56, 2.64, 2.72, 2.8, 2.88, 2.96, 3.04, 3.12, 3.2, 3.28, 3.36, 3.44, 3.52, 3.6, 3.68, 3.76, 3.84, 3.92, 4, 4, 4.08, 4.16, 4.24, 4.32, 4.4, 4.48, 4.56, 4.64, 4.72, 4.8, 4.88, 4.96, 5.04, 5.12, 5.2, 5.28, 5.36, 5.44, 5.52, 5.6, 5.68, 5.76, 5.84, 5.92, 6, 6.08, 6.16, 6.24, 6.32, 6.4, 6.48, 6.56, 6.64, 6.72, 6.8, 6.88, 6.96, 7.04, 7.12, 7.2, 7.28, 7.36, 7.44, 7.52, 7.6, 7.68, 7.76, 7.84, 7.92]
x = [
  [1, 0.999665, 0.999295, 0.998887, 0.998441, 0.997954, 0.997423, 0.996847, 0.996224, 0.995551, 0.994825, 0.994044, 0.993205, 0.992305, 0.991341, 0.990311, 0.98921, 0.988036, 0.986784, 0.985452, 0.984035, 0.982529, 0.98093, 0.979235, 0.977437, 0.975534, 0.97352, 0.97139, 0.969139, 0.966762, 0.964254, 0.961609, 0.958821, 0.955885, 0.952795, 0.949544, 0.946128, 0.942539, 0.938772, 0.93482, 0.930676, 0.926336, 0.921791, 0.917037, 0.912066, 0.906874, 0.901453, 0.895799, 0.889905, 0.883767, 0.877379, 0.877379, 0.870738, 0.863838, 0.856676, 0.849249, 0.841554, 0.833589, 0.825352, 0.816843, 0.808062, 0.799009, 0.789686, 0.780096, 0.770241, 0.760127, 0.749759, 0.739142, 0.728284, 0.717194, 0.70588, 0.694353, 0.682625, 0.670706, 0.658611, 0.646353, 0.633946, 0.621407, 0.608751, 0.595995, 0.583155, 0.57025, 0.557297, 0.544314, 0.53132, 0.518331, 0.505366, 0.492443, 0.47958, 0.466792, 0.454097, 0.441511, 0.429048, 0.416724, 0.404552, 0.392545, 0.380715, 0.369074, 0.357631, 0.346396, 0.335377],
  [0, 0.000326271, 0.000669734, 0.00103127, 0.00141181, 0.00181232, 0.00223382, 0.00267736, 0.00314407, 0.00363511, 0.00415169, 0.00469508, 0.00526662, 0.00586768, 0.00649973, 0.00716425, 0.00786282, 0.00859708, 0.00936872, 0.0101795, 0.0110313, 0.011926, 0.0128655, 0.0138519, 0.0148873, 0.0159739, 0.017114, 0.0183098, 0.0195637, 0.0208782, 0.0222558, 0.0236991, 0.0252105, 0.026793, 0.028449, 0.0301814, 0.031993, 0.0338866, 0.0358649, 0.0379309, 0.0400872, 0.0423367, 0.0446821, 0.0471262, 0.0496715, 0.0523206, 0.055076, 0.0579399, 0.0609146, 0.0640021, 0.0672042, 0.0672042, 0.0705225, 0.0739584, 0.0775129, 0.0811869, 0.0849809, 0.0888951, 0.0929293, 0.0970827, 0.101355, 0.105743, 0.110247, 0.114863, 0.119589, 0.124421, 0.129355, 0.134387, 0.139513, 0.144725, 0.150019, 0.155389, 0.160826, 0.166323, 0.171873, 0.177468, 0.183097, 0.188754, 0.194428, 0.200109, 0.205788, 0.211456, 0.217101, 0.222715, 0.228287, 0.233807, 0.239266, 0.244654, 0.249961, 0.25518, 0.260301, 0.265316, 0.270218, 0.274998, 0.279652, 0.284171, 0.28855, 0.292785, 0.296871, 0.300804, 0.304579],
  [0, 5.1079e-06, 2.07854e-05, 4.75827e-05, 8.60771e-05, 0.000136874, 0.00020061, 0.000277949, 0.00036959, 0.000476265, 0.000598738, 0.000737811, 0.000894323, 0.00106915, 0.00126321, 0.00147747, 0.00171291, 0.00197059, 0.0022516, 0.00255707, 0.00288819, 0.00324618, 0.00363234, 0.00404798, 0.0044945, 0.00497334, 0.00548597, 0.00603394, 0.00661885, 0.00724235, 0.00790612, 0.00861194, 0.0093616, 0.010157, 0.0109999, 0.0118924, 0.0128365, 0.0138342, 0.0148875, 0.0159986, 0.0171696, 0.0184026, 0.0197, 0.0210638, 0.0224962, 0.0239995, 0.0255759, 0.0272275, 0.0289564, 0.0307648, 0.0326548, 0.0326548, 0.0346282, 0.0366871, 0.0388333, 0.0410685, 0.0433942, 0.0458121, 0.0483234, 0.0509293, 0.0536309, 0.056429, 0.0593242, 0.0623169, 0.0654074, 0.0685956, 0.0718812, 0.0752638, 0.0787426, 0.0823164, 0.0859841, 0.0897439, 0.0935941, 0.0975326, 0.101557, 0.105665, 0.109853, 0.114118, 0.118457, 0.122866, 0.127343, 0.131882, 0.13648, 0.141132, 0.145835, 0.150583, 0.155374, 0.160201, 0.165061, 0.169949, 0.174861, 0.179793, 0.18474, 0.189697, 0.194662, 0.19963, 0.204596, 0.209559, 0.214513, 0.219456, 0.224385],
  [0, 3.62711e-06, 1.47587e-05, 3.37839e-05, 6.11105e-05, 9.71664e-05, 0.0001424, 0.000197279, 0.000262298, 0.000337969, 0.00042483, 0.000523445, 0.000634402, 0.000758315, 0.000895825, 0.0010476, 0.00121434, 0.00139678, 0.00159566, 0.00181178, 0.00204596, 0.00229905, 0.00257192, 0.00286551, 0.00318075, 0.00351863, 0.00388017, 0.00426641, 0.00467843, 0.00511735, 0.00558431, 0.00608048, 0.00660706, 0.00716529, 0.00775642, 0.00838172, 0.00904251, 0.00974009, 0.0104758, 0.011251, 0.0120669, 0.0129251, 0.0138268, 0.0147733, 0.0157659, 0.0168061, 0.0178949, 0.0190337, 0.0202237, 0.0214659, 0.0227615, 0.0227615, 0.0241114, 0.0255166, 0.0269778, 0.0284957, 0.030071, 0.0317041, 0.0333954, 0.035145, 0.0369529, 0.038819, 0.0407431, 0.0427245, 0.0447625, 0.0468564, 0.0490049, 0.0512068, 0.0534605, 0.0557642, 0.0581162, 0.0605141, 0.0629556, 0.0654382, 0.0679591, 0.0705154, 0.0731039, 0.0757214, 0.0783646, 0.0810299, 0.0837137, 0.0864123, 0.0891219, 0.0918387, 0.094559, 0.0972788, 0.0999944, 0.102702, 0.105398, 0.108079, 0.11074, 0.11338, 0.115994, 0.11858, 0.121134, 0.123655, 0.126138, 0.128582, 0.130985, 0.133344, 0.135659],
]
u = [
  [710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710],
]

import matplotlib.pyplot as plt
import numpy as np

# Your current x shape: (n_states, n_time_points)
t = np.array(t)
x = np.array(x)

# Plot each state over time
for i in range(x.shape[0]):
    plt.plot(t, x[i, :], label=f'x[{i+1}]')  # i+1 for 1-based indexing

plt.xlabel("Time [s]")
plt.ylabel("State value")
plt.title("State trajectories over time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
