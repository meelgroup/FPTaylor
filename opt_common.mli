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

val empty_result : opt_result

val default_opt_pars : unit -> opt_pars

val x_abs_tol_for_var : opt_pars -> string -> float

val x_tols : opt_pars -> string list -> Interval.interval array -> float array
	  
val get_float : ?default:float -> string list -> string -> float
