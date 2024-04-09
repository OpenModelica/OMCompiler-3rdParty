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

use egg::rewrite;

use crate::egraph::{ConstantFold, ModelicaExpr};

pub type RuleSet = Vec<egg::Rewrite<ModelicaExpr, ConstantFold>>;

/// Return vector of rewrite rules
pub fn make_rules() -> RuleSet {
  vec![
      rewrite!("add-commute";   "(+ ?a ?b)" => "(+ ?b ?a)"),
      rewrite!("add-associate"; "(+ (+ ?a ?b) ?c)" => "(+ ?a (+ ?b ?c))"),
      rewrite!("add-neutral";   "(+ ?a 0)" => "?a"),
      rewrite!("add-inverse";   "(- ?a ?a)" => "0"),

      rewrite!("sub-associate"; "(+ ?a (- ?b ?c))" => "(- (+ ?a ?b) ?c)"),
      rewrite!("sub-associate2"; "(- (+ ?a ?b) ?c)" => "(+ ?a (- ?b ?c))"),

      rewrite!("mul-commute"; "(* ?a ?b)" => "(* ?b ?a)"),
      rewrite!("mul-associate"; "(* (* ?a ?b) ?c)" => "(* ?a (* ?b ?c))"),
      rewrite!("mul-1"; "(* ?a 1)" => "?a"),

      rewrite!("div-associate"; "(* (/ ?a ?b) ?c)" => "(* ?a (/ ?c ?b))"),
      rewrite!("div-inv"; "(/ ?a ?a)" => "1"),

      rewrite!("add-mul-distribute"; "(+ (* ?a ?b) (* ?a ?c))" => "(* ?a (+ ?b ?c))"),
      rewrite!("mul-0"; "(* ?a 0)" => "0"),

      rewrite!("add-same-base"; "(+ ?a ?a)" => "(* ?a 2)"),
      rewrite!("add-same"; "(+ ?a (* ?a ?n))" => "(* ?a (+ ?n 1))"),

      rewrite!("mul-same-base"; "(* ?a ?a)" => "(^ ?a 2)"),
      rewrite!("mul-same"; "(* ?a (^ ?a ?n))" => "(^ ?a (+ ?n 1))"),

      rewrite!("pow-zero"; "(^ ?a 0)" => "1"),
      rewrite!("pow-one"; "(^ 1 ?a)" => "1"),
      rewrite!("pow-distribute"; "(^ (* ?a ?b) ?n)" => "(* (^ ?a ?n) (^ ?b ?n))"),

      rewrite!("sin-0"; "(sin 0)" => "0"),
  ]
}
