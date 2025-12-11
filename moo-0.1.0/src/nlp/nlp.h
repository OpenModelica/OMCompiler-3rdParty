// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2025 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MOO_NLP_H
#define MOO_NLP_H

#include <base/fixed_vector.h>
#include <base/export.h>
#include <base/log.h>

#include "nlp_scaling.h"

namespace NLP {

/**
 * @brief Generic NLP base class for optimization problems.
 *
 * This class provides a generic interface for Nonlinear Programming (NLP),
 * allowing for abstraction and interoperability with various NLP solvers.
 *
 * The NLP problem is formulated as:
 * \f{align*}{
 * \min & f(x) \\
 * \text{s.t.} & g^{LB} \le g(x) \le g^{UB} \\
 * & x^{LB} \le x \le x^{UB}
 * \f}
 *
 * Users implementing a specific NLP must derive from this class and implement the
 * pure virtual methods (User Callbacks). NLP solvers interact with the problem
 * through the 'Solver API' methods.
 */
class MOO_EXPORT NLP {
public:
    /**
     * @brief Default constructor.
     */
    NLP() = default;

    /**
     * @brief Virtual destructor.
     */
    virtual ~NLP() = default;

    // ============ User API (getters) ============

    /**
     * @brief Gets the total number of variables in the NLP.
     * @return The number of variables.
     */
    inline int get_number_vars()         const { return number_vars; }

    /**
     * @brief Gets the total number of constraints in the NLP.
     * @return The number of constraints.
     */
    inline int get_number_constraints()  const { return number_constraints; }

    /**
     * @brief Gets the number of non-zero elements in the Jacobian matrix.
     * @return The number of non-zero elements in the Jacobian.
     */
    inline int get_nnz_jac()             const { return nnz_jac; }

    /**
     * @brief Gets the number of non-zero elements in the Hessian matrix.
     * @return The number of non-zero elements in the Hessian.
     */
    inline int get_nnz_hes()             const { return nnz_hes; }

    // ============ User Callbacks (virtuals, are sorted in call order) ============

    /**
     * @brief Pure virtual method to get the number of variables and constraints.
     *
     * This method must be implemented by the user to specify the problem dimensions.
     *
     * @param[out] number_vars The total number of primal variables.
     * @param[out] number_constraints The total number of constraints.
     */
    virtual void get_sizes(
        int& number_vars,
        int& number_constraints
    ) = 0;

    /**
     * @brief Pure virtual method to get the number of non-zero elements in the Jacobian and Hessian.
     *
     * This method must be implemented by the user to provide sparsity information.
     *
     * @param[out] nnz_jac The number of non-zero elements in the Jacobian matrix.
     * @param[out] nnz_hes The number of non-zero elements in the Hessian matrix.
     */
    virtual void get_nnz(
        int& nnz_jac,
        int& nnz_hes
    ) = 0;

    /**
     * @brief Gets the scaling object for the NLP.
     *
     * This method can be overridden by users to provide custom scaling for the NLP
     * problem. By default, it returns a null shared pointer, indicating no scaling.
     *
     * @return A shared pointer to a Scaling object.
     */
    virtual std::shared_ptr<Scaling> get_scaling();


    /**
     * @brief Pure virtual method to get the bounds for variables and constraints.
     *
     * This method must be implemented by the user to define the box constraints
     * on variables and the bounds on the constraint functions.
     *
     * @param[out] x_lb Lower bounds for the primal variables \f$x\f$.
     * @param[out] x_ub Upper bounds for the primal variables \f$x\f$.
     * @param[out] g_lb Lower bounds for the constraint functions \f$g(x)\f$.
     * @param[out] g_ub Upper bounds for the constraint functions \f$g(x)\f$.
     */
    virtual void get_bounds(
        FixedVector<f64>& x_lb,
        FixedVector<f64>& x_ub,
        FixedVector<f64>& g_lb,
        FixedVector<f64>& g_ub) = 0;

    /**
     * @brief Pure virtual method to get the initial guess for the NLP.
     *
     * This method provides the initial values for the primal variables,
     * dual variables, and dual multipliers of the variable bounds.
     *
     * @param[in] init_x A boolean indicating if an initial guess for \f$x\f$ is required.
     * @param[out] x_init The initial guess for the primal variables \f$x\f$.
     * @param[in] init_lambda A boolean indicating if an initial guess for \f$\lambda\f$ is required.
     * @param[out] lambda_init The initial guess for the dual variables \f$\lambda\f$.
     * @param[in] init_z A boolean indicating if an initial guess for \f$z\f$ is required.
     * @param[out] z_lb_init The initial guess for dual multipliers of lower variable bounds \f$z_{LB}\f$.
     * @param[out] z_ub_init The initial guess for dual multipliers of upper variable bounds \f$z_{UB}\f$.
     */
    virtual void get_initial_guess(
        bool init_x,
        FixedVector<f64>& x_init,
        bool init_lambda,
        FixedVector<f64>& lambda_init,
        bool init_z,
        FixedVector<f64>& z_lb_init,
        FixedVector<f64>& z_ub_init) = 0;

    /**
     * @brief Pure virtual method to get the sparsity pattern of the Jacobian matrix.
     *
     * This method must be implemented by the user to specify the row and column
     * indices of the non-zero elements in the Jacobian of the constraints \f$\frac{dg}{dx}\f$.
     * The indices should be in COO (Coordinate) format.
     *
     * @param[out] i_row_jac The row indices of the non-zero elements.
     * @param[out] j_col_jac The column indices of the non-zero elements.
     */
    virtual void get_jac_sparsity(
        FixedVector<int>& i_row_jac,
        FixedVector<int>& j_col_jac) = 0;

    /**
     * @brief Pure virtual method to get the sparsity pattern of the Hessian matrix.
     *
     * This method must be implemented by the user to specify the row and column
     * indices of the non-zero elements in the Hessian of the Lagrangian \f$\frac{d^2\mathcal{L}}{dx^2}\f$.
     * The indices should be in COO (Coordinate) format. The Hessian is symmetric,
     * and typically only the lower or upper triangle needs to be specified.
     *
     * @param[out] i_row_hes The row indices of the non-zero elements.
     * @param[out] j_col_hes The column indices of the non-zero elements.
     */
    virtual void get_hes_sparsity(
        FixedVector<int>& i_row_hes,
        FixedVector<int>& j_col_hes) = 0;

    /**
     * @brief Pure virtual method to evaluate the objective function \f$f(x)\f$.
     *
     * @param[in] new_x A boolean indicating if the current primal variables \f$curr\_x\f$ are new and require re-evaluation.
     * @param[in] curr_x The current primal variables \f$x\f$.
     * @param[out] curr_obj The evaluated objective function value \f$f(x)\f$.
     */
    virtual void eval_f(
        bool new_x,
        const FixedVector<f64>& curr_x,
        f64& curr_obj) = 0;

    /**
     * @brief Pure virtual method to evaluate the constraint functions \f$g(x)\f$.
     *
     * @param[in] new_x A boolean indicating if the current primal variables \f$curr\_x\f$ are new and require re-evaluation.
     * @param[in] curr_x The current primal variables \f$x\f$.
     * @param[out] curr_g The evaluated constraint function values \f$g(x)\f$.
     */
        virtual void eval_g(
        bool new_x,
        const FixedVector<f64>& curr_x,
        FixedVector<f64>& curr_g) = 0;

    /**
     * @brief Pure virtual method to evaluate the gradient of the objective function \f$\frac{df}{dx}\f$.
     *
     * @param[in] new_x A boolean indicating if the current primal variables \f$curr\_x\f$ are new and require re-evaluation.
     * @param[in] curr_x The current primal variables \f$x\f$.
     * @param[out] curr_grad_f The evaluated gradient of the objective function \f$\nabla f(x)\f$.
     */
    virtual void eval_grad_f(
        bool new_x,
        const FixedVector<f64>& curr_x,
        FixedVector<f64>& curr_grad_f) = 0;

    /**
     * @brief Pure virtual method to evaluate the Jacobian of the constraint functions \f$\frac{dg}{dx}\f$.
     *
     * The non-zero elements of the Jacobian should be returned in the order corresponding
     * to the sparsity pattern provided by `get_jac_sparsity`.
     *
     * @param[in] new_x A boolean indicating if the current primal variables \f$curr\_x\f$ are new and require re-evaluation.
     * @param[in] curr_x The current primal variables \f$x\f$.
     * @param[in] i_row_jac The row indices of the non-zero elements of the Jacobian.
     * @param[in] j_col_jac The column indices of the non-zero elements of the Jacobian.
     * @param[out] curr_jac The evaluated non-zero elements of the Jacobian matrix.
     */
    virtual void eval_jac_g(
        bool new_x,
        const FixedVector<f64>& curr_x,
        const FixedVector<int>& i_row_jac,
        const FixedVector<int>& j_col_jac,
        FixedVector<f64>& curr_jac) = 0;

    /**
     * @brief Pure virtual method to evaluate the Hessian of the Lagrangian \f$\sigma_f \frac{d^2f}{dx^2} + \sum_{i=1}^{m} \lambda_i \frac{d^2(g_i)}{dx^2}\f$.
     *
     * The non-zero elements of the Hessian should be returned in the order corresponding
     * to the sparsity pattern provided by `get_hes_sparsity`.
     *
     * @param[in] new_x A boolean indicating if the current primal variables \f$curr\_x\f$ are new and require re-evaluation.
     * @param[in] curr_x The current primal variables \f$x\f$.
     * @param[in] new_lambda A boolean indicating if the current dual variables \f$curr\_lambda\f$ are new and require re-evaluation.
     * @param[in] curr_lambda The current dual variables (Lagrange multipliers) \f$\lambda\f$.
     * @param[in] curr_obj_factor The scalar multiplier \f$\sigma_f\f$ for the objective function's Hessian.
     * @param[in] i_row_hes The row indices of the non-zero elements of the Hessian.
     * @param[in] j_col_hes The column indices of the non-zero elements of the Hessian.
     * @param[out] curr_hes The evaluated non-zero elements of the Hessian matrix.
     */
    virtual void eval_hes(
        bool new_x,
        const FixedVector<f64>& curr_x,
        bool new_lambda,
        const FixedVector<f64>& curr_lambda,
        f64 curr_obj_factor,
        const FixedVector<int>& i_row_hes,
        const FixedVector<int>& j_col_hes,
        FixedVector<f64>& curr_hes) = 0;

    /**
     * @brief Pure virtual method to finalize the solution from the solver.
     *
     * This method is called by the solver after a solution has been found, allowing
     * the user to process or store the optimal values.
     *
     * @param[in] MOO_obj The optimal objective function value.
     * @param[in] MOO_x The optimal primal variables \f$x^*\f$.
     * @param[in] MOO_lambda The optimal dual variables (Lagrange multipliers) \f$\lambda^*\f$.
     * @param[in] MOO_z_lb The optimal dual multipliers for the lower variable bounds \f$z_{LB}^*\f$.
     * @param[in] MOO_z_ub The optimal dual multipliers for the upper variable bounds \f$z_{UB}^*\f$.
     */
    virtual void finalize_solution(
        f64 MOO_obj,
        const FixedVector<f64>& opt_x,
        const FixedVector<f64>& opt_lambda,
        const FixedVector<f64>& opt_z_lb,
        const FixedVector<f64>& opt_z_ub) = 0;

    // ============ Solver API ============

    /**
     * @brief Retrieves problem dimensions for the solver.
     *
     * This method is the first point of interaction for an NLP solver. It queries the user-defined
     * NLP implementation for the total number of variables, constraints, and the number of
     * non-zero elements in the Jacobian and Hessian matrices.
     *
     * It initializes internal buffers and sets up the scaling mechanism based on the user's
     * `get_scaling()` implementation. If `get_scaling()` returns a `nullptr`, it defaults to `NoScaling`.
     *
     * @param[out] solver_number_vars Total number of primal variables \f$n\f$.
     * @param[out] solver_number_constraints Total number of constraints \f$m\f$.
     * @param[out] solver_nnz_jac Number of non-zero elements in the Jacobian matrix \f$\frac{dg}{dx}\f$.
     * @param[out] solver_nnz_hes Number of non-zero elements in the Hessian matrix \f$\frac{d^2\mathcal{L}}{dx^2}\f$.
     */
    void solver_get_info(
        int& solver_number_vars,
        int& solver_number_constraints,
        int& solver_nnz_jac,
        int& solver_nnz_hes
    );

    /**
     * @brief Retrieves the bounds for variables and constraints for the solver.
     *
     * This method obtains the lower and upper bounds for the primal variables \f$x\f$
     * and the constraint functions \f$g(x)\f$ from the user-defined NLP.
     *
     * The bounds are then **scaled** according to the `Scaling` object before being
     * returned to the solver.
     *
     * @param[out] solver_x_lb Pointer to an array where the scaled lower bounds of variables \f$x^{LB}\f$ will be written.
     * @param[out] solver_x_ub Pointer to an array where the scaled upper bounds of variables \f$x^{UB}\f$ will be written.
     * @param[out] solver_g_lb Pointer to an array where the scaled lower bounds of constraints \f$g^{LB}\f$ will be written.
     * @param[out] solver_g_ub Pointer to an array where the scaled upper bounds of constraints \f$g^{UB}\f$ will be written.
     */
    void solver_get_bounds(
        f64* solver_x_lb,
        f64* solver_x_ub,
        f64* solver_g_lb,
        f64* solver_g_ub);

    /**
     * @brief Retrieves the initial guess for the NLP problem for the solver.
     *
     * This method obtains the initial guess for the primal variables \f$x\f$,
     * dual variables \f$\lambda\f$, and dual multipliers for variable bounds \f$z_{LB}, z_{UB}\f$
     * from the user-defined NLP.
     *
     * The initial guess for \f$x\f$ and \f$z\f$ are **scaled**, and \f$\lambda\f$ is **unscaled**
     * before being returned to the solver, to match the solver's expected input format.
     *
     * @param[in] init_x A boolean flag indicating if the solver requires an initial guess for \f$x\f$.
     * @param[out] solver_x_init Pointer to an array where the scaled initial guess for primal variables \f$x\f$ will be written.
     * @param[in] init_lambda A boolean flag indicating if the solver requires an initial guess for \f$\lambda\f$.
     * @param[out] solver_lambda_init Pointer to an array where the unscaled initial guess for dual variables \f$\lambda\f$ will be written.
     * @param[in] init_z A boolean flag indicating if the solver requires an initial guess for \f$z\f$.
     * @param[out] solver_z_lb_init Pointer to an array where the scaled initial guess for dual multipliers of lower variable bounds \f$z_{LB}\f$ will be written.
     * @param[out] solver_z_ub_init Pointer to an array where the scaled initial guess for dual multipliers of upper variable bounds \f$z_{UB}\f$ will be written.
     */
    void solver_get_initial_guess(
        bool init_x,
        f64* solver_x_init,
        bool init_lambda,
        f64* solver_lambda_init,
        bool init_z,
        f64* solver_z_lb_init,
        f64* solver_z_ub_init);

    /**
     * @brief Retrieves the sparsity pattern of the Jacobian matrix for the solver.
     *
     * This method obtains the row and column indices of the non-zero elements
     * of the Jacobian matrix \f$\frac{dg}{dx}\f$ from the user-defined NLP.
     * The indices are provided in COO (Coordinate) format.
     *
     * @param[out] solver_i_row_jac Pointer to an array where the row indices of Jacobian non-zeros will be written.
     * @param[out] solver_j_col_jac Pointer to an array where the column indices of Jacobian non-zeros will be written.
     */
    void solver_get_jac_sparsity(
        int* solver_i_row_jac,
        int* solver_j_col_jac
    );

    /**
     * @brief Retrieves the sparsity pattern of the Hessian matrix for the solver.
     *
     * This method obtains the row and column indices of the non-zero elements
     * of the Hessian matrix \f$\frac{d^2\mathcal{L}}{dx^2}\f$ from the user-defined NLP.
     * The indices are provided in COO (Coordinate) format.
     *
     * @param[out] solver_i_row_hes Pointer to an array where the row indices of Hessian non-zeros will be written.
     * @param[out] solver_j_col_hes Pointer to an array where the column indices of Hessian non-zeros will be written.
     */
    void solver_get_hes_sparsity(
        int* solver_i_row_hes,
        int* solver_j_col_hes
    );

    /**
     * @brief Evaluates the objective function \f$f(x)\f$ for the solver.
     *
     * This method first unscales the primal variables \f$x\f$ received from the solver,
     * then calls the user-defined `eval_f` method to compute the objective value.
     * Finally, the computed objective value is **scaled** before being returned to the solver.
     *
     * @param[in] new_x A boolean flag indicating if the current primal variables \f$solver\_x\f$ are new and require re-evaluation.
     * @param[in] solver_x Pointer to an array containing the current primal variables \f$x\f$ (scaled by the solver).
     * @param[out] solver_obj_value The evaluated objective function value \f$f(x)\f$, scaled for the solver.
     */
    void solver_eval_f(
        bool new_x,
        const f64* solver_x,
        f64& solver_obj_value);

    /**
     * @brief Evaluates the gradient of the objective function \f$\frac{df}{dx}\f$ for the solver.
     *
     * This method first unscales the primal variables \f$x\f$ received from the solver,
     * then calls the user-defined `eval_grad_f` method to compute the gradient.
     * Finally, the computed gradient is **scaled** before being returned to the solver.
     *
     * @param[in] new_x A boolean flag indicating if the current primal variables \f$solver\_x\f$ are new and require re-evaluation.
     * @param[in] solver_x Pointer to an array containing the current primal variables \f$x\f$ (scaled by the solver).
     * @param[out] solver_grad_f Pointer to an array where the scaled gradient of the objective function \f$\nabla f(x)\f$ will be written.
     */
    void solver_eval_grad_f(
        bool new_x,
        const f64* solver_x,
        f64* solver_grad_f);

    /**
     * @brief Evaluates the constraint functions \f$g(x)\f$ for the solver.
     *
     * This method first unscales the primal variables \f$x\f$ received from the solver,
     * then calls the user-defined `eval_g` method to compute the constraint values.
     * Finally, the computed constraint values are **scaled** before being returned to the solver.
     *
     * @param[in] new_x A boolean flag indicating if the current primal variables \f$solver\_x\f$ are new and require re-evaluation.
     * @param[in] solver_x Pointer to an array containing the current primal variables \f$x\f$ (scaled by the solver).
     * @param[out] solver_g Pointer to an array where the scaled constraint function values \f$g(x)\f$ will be written.
     */
    void solver_eval_g(
        bool new_x,
        const f64* solver_x,
        f64* solver_g);

    /**
     * @brief Evaluates the Jacobian of the constraint functions \f$\frac{dg}{dx}\f$ for the solver.
     *
     * This method first unscales the primal variables \f$x\f$ received from the solver,
     * then calls the user-defined `eval_jac_g` method to compute the non-zero elements of the Jacobian.
     * Finally, the computed Jacobian elements are **scaled** before being returned to the solver.
     *
     * @param[in] new_x A boolean flag indicating if the current primal variables \f$solver\_x\f$ are new and require re-evaluation.
     * @param[in] solver_x Pointer to an array containing the current primal variables \f$x\f$ (scaled by the solver).
     * @param[out] solver_jac Pointer to an array where the scaled non-zero elements of the Jacobian matrix will be written.
     */
    void solver_eval_jac(
        bool new_x,
        const f64* solver_x,
        f64* solver_jac);

    /**
     * @brief Evaluates the Hessian of the Lagrangian \f$\sigma_f \frac{d^2f}{dx^2} + \sum_{i=1}^{m} \lambda_i \frac{d^2(g_i)}{dx^2}\f$ for the solver.
     *
     * This method first unscales the primal variables \f$x\f$ and scales the dual variables \f$\lambda\f$
     * and objective factor \f$\sigma_f\f$ received from the solver. It then calls the user-defined
     * `eval_hes` method to compute the non-zero elements of the Hessian.
     * Finally, the computed Hessian elements are **scaled** before being returned to the solver.
     *
     * @param[in] new_x A boolean flag indicating if the current primal variables \f$solver\_x\f$ are new and require re-evaluation.
     * @param[in] solver_x Pointer to an array containing the current primal variables \f$x\f$ (scaled by the solver).
     * @param[in] new_lambda A boolean flag indicating if the current dual variables \f$solver\_lambda\f$ are new and require re-evaluation.
     * @param[in] solver_lambda Pointer to an array containing the current dual variables (Lagrange multipliers) \f$\lambda\f$ (scaled by the solver).
     * @param[in] solver_obj_factor The scalar multiplier \f$\sigma_f\f$ for the objective function's Hessian (scaled by the solver).
     * @param[out] solver_hes Pointer to an array where the scaled non-zero elements of the Hessian matrix will be written.
     */
    void solver_eval_hes(
        bool new_x,
        const f64* solver_x,
        bool new_lambda,
        const f64* solver_lambda,
        const f64 solver_obj_factor,
        f64* solver_hes);

    /**
     * @brief Finalizes the solution received from the solver.
     *
     * This method is called by the solver once an optimal solution has been found.
     * It first **unscales** all the optimal solution components (objective value,
     * primal variables \f$x^*\f$, dual variables \f$\lambda^*\f$, and dual multipliers
     * of variable bounds \f$z_{LB}^*, z_{UB}^*\f$).
     *
     * After unscaling, it calls the user-defined `finalize_solution` method,
     * allowing the user to perform any necessary post-processing or storage of the results
     * in their original, unscaled form.
     *
     * @param[in] solver_obj_value The optimal objective function value from the solver (scaled).
     * @param[in] solver_x Pointer to an array containing the optimal primal variables \f$x^*\f$ (scaled).
     * @param[in] solver_lambda Pointer to an array containing the optimal dual variables (Lagrange multipliers) \f$\lambda^*\f$ (scaled).
     * @param[in] solver_z_lb Pointer to an array containing the optimal dual multipliers for lower variable bounds \f$z_{LB}^*\f$ (scaled).
     * @param[in] solver_z_ub Pointer to an array containing the optimal dual multipliers for upper variable bounds \f$z_{UB}^*\f$ (scaled).
     */
    void solver_finalize_solution(
        const f64  solver_obj_value,
        const f64* solver_x,
        const f64* solver_lambda,
        const f64* solver_z_lb,
        const f64* solver_z_ub);

    /**
     * @brief Allows the user to access the value of the objective function \f$f(x)\f$.
     */
    f64 get_objective_value() const;

private:

    // ============ NLP Structures and Info ============

    int number_vars = 0;        // total number of variables in the NLP
    int number_constraints = 0; // total number of constraints in the NLP
    int nnz_jac = 0;            // nnz Jacobian in the NLP
    int nnz_hes = 0;            // nnz Hessian in the NLP

    // current iterates
    FixedVector<f64> curr_x;       // current NLP primal variables
    FixedVector<f64> curr_lambda;  // current NLP dual variables
    f64              curr_sigma_f; // current objective weight in hessian

    // variable bounds
    FixedVector<f64> x_lb; // lower bounds on NLP variables
    FixedVector<f64> x_ub; // upper bounds on NLP variables

    // optimal dual multipliers (set after finalize)
    FixedVector<f64> z_lb; // dual multipliers of lower bounds on NLP variables
    FixedVector<f64> z_ub; // dual multipliers of upper bounds on NLP variables

    // nlp function data
    f64 curr_obj;               // current NLP objective value
    FixedVector<f64> curr_grad; // current NLP gradient of the objective function
    FixedVector<f64> curr_g;    // current NLP constraint function evaluation
    FixedVector<f64> curr_jac;  // current NLP jacobian of the constraints
    FixedVector<f64> curr_hes;  // current NLP hessian of the lagrangian

    // scaled constraint bounds
    FixedVector<f64> g_lb; // lower bounds on NLP constraints
    FixedVector<f64> g_ub; // upper bounds on NLP constraints

    // COO sparsity patterns
    FixedVector<int> i_row_jac; // row COO of the Jacobian
    FixedVector<int> j_col_jac; // column COO of the Jacobian
    FixedVector<int> i_row_hes; // row COO of the Hessian
    FixedVector<int> j_col_hes; // column COO of the Hessian

    // generic scaling routine
    std::shared_ptr<Scaling> scaling; // scaling object, create custom by implemention get_scaling();

    // TODO: add a generic block BFGS routine, which can calculate blocks of the Lagrangian Hessian $\nabla_{xx} \mathcal{L}_{AA -> BB}$

    // ============ Internal Methods ============

    void allocate_buffers();
    void allocate_sparsity_buffers();

    // ============ Helpers for Scaling ============

    // Attention: this part is super confusing, since the lambda and sigma_f *unscale* is actually a scale with g and f!!
    // this is the case, because the Lagrangian Hessian does not allow for to scale f and g independently after evaluation
    // so the only possibility is to apply this scaling a priori by updating the duals and sigma_f

    inline void unscale_dual_bounds(const f64* z_L, const f64* z_U) {
        scaling->unscale_x(z_L, z_lb.raw(), number_vars);
        scaling->unscale_x(z_U, z_ub.raw(), number_vars);
    }

    inline void update_unscale_curr_x(bool new_x, const f64* x) {
        if (new_x) scaling->unscale_x(x, curr_x.raw(), number_vars);
    }

    inline void update_unscale_curr_lambda(bool new_lambda, const f64* lambda) {
        if (new_lambda) scaling->scale_g(lambda, curr_lambda.raw(), number_constraints);
    }

    inline void unscale_curr_lambda(f64* lambda) {
        scaling->unscale_g(curr_lambda.raw(), lambda, number_constraints);
    }

    inline void unscale_curr_sigma_f(const f64* sigma_f) {
        scaling->scale_f(sigma_f, &curr_sigma_f);
    }

    inline void unscale_objective(const f64* obj) {
        scaling->unscale_f(obj, &curr_obj);
    }
};

} // namespace NLP

#endif  // MOO_NLP_H
