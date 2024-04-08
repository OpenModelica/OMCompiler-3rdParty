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

use test_cdylib;
use dlopen2::wrapper::{Container, WrapperApi};
use std::os::raw::c_char;
use std::ffi::{CString, CStr};

use egg_interface::{RuleSet, Runner};

#[derive(WrapperApi)]
struct Api {
    egg_make_rules: unsafe extern "C" fn() -> Box<RuleSet>,
    egg_free_rules: unsafe extern "C" fn(_rules: Option<Box<RuleSet>>),
    egg_make_runner: unsafe extern "C" fn() -> Box<Runner>,
    egg_free_runner: unsafe extern "C" fn(_runner: Option<Box<Runner>>),
    egg_simplify_expr: unsafe extern "C" fn(rules: Option<&RuleSet>, runner_ptr: Option<&mut Runner>, expr_str: *const c_char) -> *mut c_char
}

/// Test API as called from OpenModelica.
#[test]
fn omc_test() {
    let dylib_path = test_cdylib::build_current_project();

    let container: Container<Api> = unsafe { Container::load(dylib_path) }.expect("Could not open library or load symbols");

    let rules = unsafe { container.egg_make_rules() };
    let rules_pt =  Some(&(*rules));
    let mut runner = unsafe { container.egg_make_runner() };
    let runner_pt =  Some(&mut(*runner));

    let expression = CString::new("(+ x (+ y (* 2 x)))").unwrap();
    let expression: *const c_char = expression.as_ptr() as *const c_char;
    let simplified_expression: *mut c_char = unsafe { container.egg_simplify_expr(rules_pt, runner_pt, expression) };

    // TODO: How to test if free functions are working?
    unsafe { container.egg_free_rules(Some(rules)) };
    unsafe { container.egg_free_runner(Some(runner)) };

    let result = unsafe { CStr::from_ptr(simplified_expression).to_string_lossy().into_owned() };
    println!("Simplified expression: {:?}", result);
    assert_eq!(result, "(+ y (* x 3))");
}
