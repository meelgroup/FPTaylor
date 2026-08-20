(* ========================================================================== *)
(*      FPTaylor: A Tool for Rigorous Estimation of Round-off Errors          *)
(*                                                                            *)
(*      Author: Alexey Solovyev, University of Utah                           *)
(*                                                                            *)
(*      This file is distributed under the terms of the MIT license           *)
(* ========================================================================== *)

(* -------------------------------------------------------------------------- *)
(* Common optimization functions and types                                    *)
(* -------------------------------------------------------------------------- *)

(* Optimization parameters (not all parameters are supported by individual backends) *)
type opt_pars = {
  f_rel_tol : float;
  f_abs_tol : float;
  x_rel_tol : float;
  x_abs_tol : float;
  x_abs_tol_vars : (string * float) list;
  max_iters : int;
  timeout : int;
}

type opt_result = {
  result : float;
  lower_bound : float;
  iters : int;
  time : float;
}

let empty_result = {
  result = 0.;
  lower_bound = 0.;
  iters = 0;
  time = 0.;
}

let get_tol name default =
  let tol = Config.get_float_option name in
  if tol < 0.0 then
    let () = Log.warning "Bad tolerance value: %s = %e. Using default value: %e"
               name tol default in
    default
  else
    tol

let parse_x_abs_tol_vars () =
  let raw = String.trim (Config.get_string_option "opt-x-abs-tol-vars") in
  if raw = "" then
    []
  else
    let parse_entry entry =
      match Str.bounded_split (Str.regexp "[:=]") entry 2 with
      | [name; value] ->
        let name = String.trim name and value = String.trim value in
        if name = "" then
          failwith "Empty variable name in opt-x-abs-tol-vars";
        let tol =
          try float_of_string value
          with _ ->
            failwith ("Cannot convert opt-x-abs-tol-vars value into a float: " ^ entry) in
        if tol < 0.0 then
          failwith ("Negative tolerance in opt-x-abs-tol-vars: " ^ entry);
        name, tol
      | _ ->
        failwith ("Expected NAME=VALUE in opt-x-abs-tol-vars entry: " ^ entry) in
    Str.split (Str.regexp "[,; \t]+") raw |> List.map parse_entry

let default_opt_pars () = {
  f_rel_tol = get_tol "opt-f-rel-tol" 0.01;
  f_abs_tol = get_tol "opt-f-abs-tol" 0.01;
  x_rel_tol = get_tol "opt-x-rel-tol" 0.0;
  x_abs_tol = get_tol "opt-x-abs-tol" 0.01;
  x_abs_tol_vars = parse_x_abs_tol_vars ();
  max_iters = Config.get_int_option "opt-max-iters";
  timeout = Config.get_int_option "opt-timeout";
}

let x_abs_tol_for_var pars name =
  try List.assoc name pars.x_abs_tol_vars
  with Not_found -> pars.x_abs_tol

let x_tols pars var_names var_bounds =
  let r = Interval.size_max_X var_bounds in
  Array.mapi
    (fun i _ -> r *. pars.x_rel_tol +. x_abs_tol_for_var pars (List.nth var_names i))
    var_bounds

let get_float ?default strs name =
  let pat = name in
  let n = String.length pat in
  let r = Str.regexp pat in
  let rec find = function
    | [] ->
       if Lib.is_none default then raise Not_found
       else Lib.option_value default
    | str :: t ->
      let i = try Str.search_forward r str 0 with Not_found -> -1 in
      if i == 0 then
        float_of_string (Str.string_after str n)
      else
        find t in
        find strs
