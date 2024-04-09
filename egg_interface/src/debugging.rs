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

/*----------------------------------------------------------------------------*/
/* useful functions between Rust and C                                        */
/*----------------------------------------------------------------------------*/

use std::slice;
use std::os::raw::c_char;
use std::ffi::{c_void, CStr};

/// A Rust struct mapping the C struct
#[repr(C)]
#[derive(Debug)]
pub struct RustStruct {
    pub c: char,
    pub ul: u64,
    pub c_string: *const c_char,
}

macro_rules! create_function {
    // This macro takes an argument of designator `ident` and
    // creates a function named `$func_name`.
    // The `ident` designator is used for variable/function names.
    ($func_name:ident, $ctype:ty) => {
        #[no_mangle]
        pub extern "C" fn $func_name(v: $ctype) {
            // The `stringify!` macro converts an `ident` into a string.
            log::trace!(
                "{:?}() is called, value passed = <{:?}>",
                stringify!($func_name),
                v
            );
        }
    };
}

// create simple functions where C type is exactly mapping a Rust type
//create_function!(rust_char, char);
//create_function!(rust_wchar, char);
create_function!(rust_short, i16);
create_function!(rust_ushort, u16);
create_function!(rust_int, i32);
create_function!(rust_uint, u32);
create_function!(rust_long, i64);
create_function!(rust_ulong, u64);
create_function!(rust_void, *mut c_void);

// for NULL-terminated C strings, it's a little bit clumsier
#[no_mangle]
pub extern "C" fn rust_string(c_string: *const c_char) {
    // build a Rust string from C string
    let s = unsafe { CStr::from_ptr(c_string).to_string_lossy().into_owned() };

    log::trace!("rust_string() is called, value passed = <{:?}>", s);
}

// for C arrays, need to pass array size
#[no_mangle]
pub extern "C" fn rust_int_array(c_array: *const i32, length: usize) {
    // build a Rust array from array & length
    let rust_array: &[i32] = unsafe { slice::from_raw_parts(c_array, length as usize) };
    log::trace!(
        "rust_int_array() is called, value passed = <{:?}>",
        rust_array
    );
}

#[no_mangle]
pub extern "C" fn rust_string_array(c_array: *const *const c_char, length: usize) {
    // build a Rust array from array & length
    let tmp_array: &[*const c_char] = unsafe { slice::from_raw_parts(c_array, length as usize) };

    // convert each element to a Rust string
    let rust_array: Vec<_> = tmp_array
        .iter()
        .map(|&v| unsafe { CStr::from_ptr(v).to_string_lossy().into_owned() })
        .collect();

    log::trace!(
        "rust_string_array() is called, value passed = <{:?}>",
        rust_array
    );
}

// for C structs, need to convert each individual Rust member if necessary
#[no_mangle]
pub unsafe extern "C" fn rust_cstruct(c_struct: *mut RustStruct) {
    let rust_struct = &*c_struct;
    let s = CStr::from_ptr(rust_struct.c_string)
        .to_string_lossy()
        .into_owned();

    log::trace!(
        "rust_cstruct() is called, value passed = <{} {} {}>",
        rust_struct.c, rust_struct.ul, s
    );
}
