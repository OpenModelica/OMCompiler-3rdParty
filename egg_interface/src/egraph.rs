/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2024, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF AGPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.8.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GNU AGPL
 * VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs:
 * http://www.openmodelica.org or
 * https://github.com/OpenModelica/ or
 * http://www.ida.liu.se/projects/OpenModelica,
 * and in the OpenModelica distribution.
 *
 * GNU AGPL version 3 is obtained from:
 * https://www.gnu.org/licenses/licenses.html#GPL
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

//use egg::{define_language, Id, Symbol, Analysis, DidMerge};
use egg::*;
use num_traits::{Pow, Zero};
use ordered_float::NotNan;

pub type EGraph = egg::EGraph<ModelicaExpr, ConstantFold>;


/// Constant needs to implement `Ord` so we can't just use `f64`
pub type Constant = NotNan<f64>;

define_language! {
    pub enum ModelicaExpr {
        "+" = Add([Id; 2]),
        "-" = Sub([Id; 2]),
        "*" = Mul([Id; 2]),
        "/" = Div([Id; 2]),
        "^" = Pow([Id; 2]),
        "der" = Der(Id),
        "sin" = Sin(Id),
        Constant(Constant),
        Symbol(Symbol),
    }
}

#[derive(Default)]
pub struct ConstantFold;
impl Analysis<ModelicaExpr> for ConstantFold {
    type Data = Option<Constant>;

    fn make(egraph: &EGraph, enode: &ModelicaExpr) -> Self::Data {
        let x = |i: &Id| egraph[*i].data;
        Some(match enode {
            ModelicaExpr::Add([a, b]) => x(a)? + x(b)?,
            ModelicaExpr::Sub([a, b]) => x(a)? - x(b)?,
            ModelicaExpr::Mul([a, b]) => x(a)? * x(b)?,
            ModelicaExpr::Div([a, b]) if !x(b)?.is_zero() => x(a)? / x(b)?,
            ModelicaExpr::Pow([a, b]) => x(a)?.pow(x(b)?),
            ModelicaExpr::Der(a) if x(a).is_some() => Constant::zero(),
            ModelicaExpr::Sin(a) if x(a).is_some() => NotNan::new(x(a)?.sin()).unwrap(),
            ModelicaExpr::Constant(n) => *n,
            _ => return None,
        })
    }

    fn merge(&mut self, to: &mut Self::Data, from: Self::Data) -> DidMerge {
        egg::merge_max(to, from)
    }

    fn modify(egraph: &mut EGraph, id: Id) {
        if let Some(i) = egraph[id].data {
            let added = egraph.add(ModelicaExpr::Constant(i));
            egraph.union(id, added);
            // to not prune, comment this out
            //egraph[id].nodes.retain(|n| n.is_leaf());
        }
    }
}
