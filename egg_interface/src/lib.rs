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

use egg::{AstSize, CostFunction, Extractor, RecExpr};
use simple_logger;
use std::ffi::{CStr, CString};
use std::mem;
use std::os::raw::c_char;
use std::time::{Duration, Instant};

pub mod runner;

pub mod rules;
use crate::rules::*;

pub mod egraph;
use crate::egraph::ModelicaExpr;

mod debugging;

/// Make the vector of rewrite rules.
#[no_mangle]
pub extern "C" fn egg_make_rules() -> Box<RuleSet> {
    // TODO: Move to a better location
    simple_logger::init_with_level(log::Level::Info).unwrap();

    let now = Instant::now();
    let rules = rules::make_rules();
    log::debug!("made rules: {:.2?}", now.elapsed());
    Box::new(rules)
}

#[no_mangle]
pub unsafe extern "C" fn egg_free_rules(_rules: Option<Box<RuleSet>>) {
    // dropped implicitly
    log::debug!("dropped rules");
}

/// Make the runner.
#[no_mangle]
pub extern "C" fn egg_make_runner() -> Box<runner::Runner> {
    let now = Instant::now();
    let runner = runner::make_runner();
    log::debug!("made runner: {:.2?}", now.elapsed());
    Box::new(runner)
}

/// Free runner.
#[no_mangle]
pub unsafe extern "C" fn egg_free_runner(_runner: Option<Box<runner::Runner>>) {
    // dropped implicitly
    log::debug!("dropped runner");
}

/// Simplify expression string `expr_str`.
///
/// Expect `expr_str` to be in prefix notation:
/// `"x + 1"` -> `"(+ x 1)"`.
#[no_mangle]
pub extern "C" fn egg_simplify_expr(rules: Option<&RuleSet>, runner_ptr: Option<&mut runner::Runner>, expr_str: *const c_char) -> *mut c_char {
    let mut times: Vec<(Duration, String)> = Vec::new();

    // parse the expression, the type annotation tells it which language to use
    let now = Instant::now();
    let expr = unsafe { CStr::from_ptr(expr_str).to_string_lossy().into_owned() };
    let expr: RecExpr<ModelicaExpr> = expr.parse().unwrap();
    times.push((now.elapsed(), String::from("expr     ")));

    // simplify the expression using a Runner, which creates an e-graph with
    // the given expression and runs the given rules over it
    let best = simplify_expr(&rules.unwrap(), runner_ptr.unwrap(), expr, &mut times);

    times.sort_by(|(a,_), (b,_)| b.cmp(a));
    log::info!("{}", times.iter().fold(String::new(), |acc, (t,s)| acc + &format!("{}\t{:.2?}", s, t) + "\n"));

    CString::new(best.to_string()).expect("return string error").into_raw()
}

fn simplify_expr(rules: &RuleSet, runner_ptr: &mut runner::Runner, expr: RecExpr<ModelicaExpr>, times: &mut Vec<(Duration, String)>) -> RecExpr<ModelicaExpr> {
    let now = Instant::now();
    let cost = AstSize.cost_rec(&expr);
    times.push((now.elapsed(), String::from("cost     ")));

    // we need a variable for the unwrapped ptr so we can swap back at the end
    let mut runner: runner::Runner = mem::replace(runner_ptr, *egg_make_runner());

    // reset runner so it can run again
    runner.stop_reason = None;

    runner = runner.with_expr(&expr).run(rules);
    times.push((now.elapsed(), String::from("runner   ")));
    match runner.stop_reason {
        Some(ref reason) => log::info!("stop reason: {:?}", reason),
        _ => ()
    }

    // the Runner knows which e-class the expression given with `with_expr` is in
    let now = Instant::now();
    let root = *runner.roots.last().unwrap();
    times.push((now.elapsed(), String::from("root     ")));

    // use an Extractor to pick the best element of the root eclass
    let now = Instant::now();
    let extractor = Extractor::new(&runner.egraph, AstSize);
    times.push((now.elapsed(), String::from("extractor")));

    let now = Instant::now();
    let (best_cost, best) = extractor.find_best(root);
    times.push((now.elapsed(), String::from("best     ")));

    log::info!("cost: {} -> {}", cost, best_cost);
    //log::info!("expr {}\n  -> {}", expr, best);

    // give runner back to metamodelica
    mem::swap(runner_ptr, &mut runner);

    return best;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assert simplification of `expression` equals one of `expected` expressions.
    ///
    /// Provide expression in prefix notation.
    /// Panics if simplification doesn't match any of the expected expressions.
    fn assert_simplify(expression: &str, expected: Vec<&str>) {
        let rules = rules::make_rules();
        let mut runner = runner::make_runner();
        let expr: RecExpr<ModelicaExpr> = expression.parse().unwrap();
        let mut times: Vec<(Duration, String)> = Vec::new();
        let result = simplify_expr(&rules, &mut runner, expr, &mut times);
        let result = result.to_string();
        assert!(expected.iter().any(|ex| result == *ex), "Expected solution not found. Found {:?}, expected one of {:?}", expression, expected);
    }

    /// Test `simplify_expr` with various expressions.
    #[test]
    fn simplify_expr_test() {
        assert_simplify("(+ x (+ y (* 2 x)))",
                        vec!["(+ y (* x 3))", "(+ y (* 3 x))", "(+ (* x 3) y)", "(+ (* 3 x) y)"]);
        assert_simplify("(- (+ x0 (+ x1 (+ x2 x3))) (+ x3 x1))",
                        vec!["(+ x0 x2)", "(+ x2 x0)"]);
    }
}
