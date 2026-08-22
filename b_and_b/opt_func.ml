open Interval

let asinh x = log (x +. sqrt (x *. x +. 1.0))

let acosh x = log (x +. sqrt (x *. x -. 1.0))

let atanh x = 0.5 *. log ((1.0 +. x) /. (1.0 -. x))

let asinh_I x = {
  low = 
    (let sqrt = Fpu.fsqrt_low (Fpu.fadd_low (Fpu.fmul_low x.low x.low) 1.0) in
     Fpu.flog_low (Fpu.fadd_low x.low sqrt));
  high =
    (let sqrt = Fpu.fsqrt_high (Fpu.fadd_high (Fpu.fmul_high x.high x.high) 1.0) in
     Fpu.flog_high (Fpu.fadd_high x.high sqrt));
}

let acosh_I x = 
  if x.high < 1.0 then failwith "acosh_I"
  else {
    low =
      (if x.low <= 1.0 then 0.0
       else
	  let sqrt = Fpu.fsqrt_low (Fpu.fsub_low (Fpu.fmul_low x.low x.low) 1.0) in
	  Fpu.flog_low (Fpu.fadd_low x.low sqrt));
    high = 
      (let sqrt = Fpu.fsqrt_high (Fpu.fsub_high (Fpu.fmul_high x.high x.high) 1.0) in
       Fpu.flog_high (Fpu.fadd_high x.high sqrt));
  }

let atanh_I x =
  if x.high < -1.0 || x.low > 1.0 then failwith "atanh_I"
  else {
    low = 
      (if x.low <= -1.0 then neg_infinity
       else
	  let t = Fpu.fdiv_low (Fpu.fadd_low 1.0 x.low) (Fpu.fsub_high 1.0 x.low) in
	  Fpu.fmul_low 0.5 (Fpu.flog_low t));
    high =
      (if x.high >= 1.0 then infinity
       else
	  let t = Fpu.fdiv_high (Fpu.fadd_high 1.0 x.high) (Fpu.fsub_low 1.0 x.high) in
	  Fpu.fmul_high 0.5 (Fpu.flog_high t));
  }	

(* For a given positive floating-point number f,
   returns the largest floating-point number 2^n such that
   2^n < f.
*)
let floor_power2 =
  let p2 f =
    let s, q = frexp f in
    if s = 0.5 then
      ldexp 1.0 (q - 2)
    else
      ldexp 1.0 (q - 1) in
  fun f ->
    match (classify_float f) with
      | FP_zero -> f
      | FP_infinite -> f
      | FP_nan -> f
      | FP_normal | FP_subnormal -> 
	if f < 0.0 then -.p2 (-.f) else p2 f

let floor_power2_I x = {
  low = floor_power2 x.low;
  high = floor_power2 x.high
}

let goldberg_ulp (prec, e_min) =
  let ulp f =
    let _, e = frexp f in
    ldexp 1.0 (max e (e_min + 1) - prec) in
  fun f ->
    match (classify_float f) with
    | FP_zero | FP_infinite | FP_nan -> f
    | FP_subnormal | FP_normal ->
      if f < 0. then -.ulp(-.f) else ulp f

let goldberg_ulp_I pars x = {
  low = goldberg_ulp pars x.low;
  high = goldberg_ulp pars x.high
}

let sub2 (x, y) = 
  if (0.5 *. x <= y && y <= 2.0 *. x) || 
     (2.0 *. x <= y && y <= 0.5 *. x) 
  then 0.0 
  else x -. y

let sub2_I (x, y) = {
  low = if (0.5 *. x.low <= y.high && y.high <= 2.0 *. x.low) ||
           (2.0 *. x.low <= y.high && y.high <= 0.5 *. x.low) 
        then 0.0 
        else Fpu.fsub_low x.low y.high;
  high = if (0.5 *. x.high <= y.low && y.low <= 2.0 *. x.high) ||
            (2.0 *. x.high <= y.low && y.low <= 0.5 *. x.high)
         then 0.0 
         else Fpu.fsub_high x.high y.low;
}

let abs_err (t, x) =
  if x >= t then 1.
  else if x <= -.t then -1.
  else 0. (* should be [-1, 1] but we cannot return an interval here; 
             this function is not important *)

let neg_one_I = {low = -1.; high = -1.}
let neg_one_one_I = {low = -1.; high = 1.}
    
let abs_err_I (t, x) =
  if x.low >= t.high then one_I
  else if x.high <= t.low then neg_one_I
  else neg_one_one_I

(* -------------------------------------------------------------------------- *)
(* Rigorous lgamma / digamma support, for x > 0 only.                         *)
(*                                                                            *)
(* This is a duplicate of the implementation in the top-level func.ml, kept  *)
(* in sync by hand: this file is compiled standalone (see                    *)
(* b_and_b/compile.template) together with each generated *.ml expression    *)
(* file, for the external "bb" branch-and-bound optimizer backend, and       *)
(* cannot simply reference the main fptaylor binary's Func module. See       *)
(* func.ml for the full derivation notes and numeric validation.             *)
(* -------------------------------------------------------------------------- *)

let gamma_stirling_T = 20.0

let gamma_bernoulli = [|
  (1.0, 6.0);
  (-1.0, 30.0);
  (1.0, 42.0);
  (-1.0, 30.0);
  (5.0, 66.0);
  (-691.0, 2730.0);
|]

let gamma_b14 = (7.0, 6.0)

let thin_I (x : float) = { low = x; high = x }

let rat_I (n : float) (d : float) =
  { low = Fpu.fdiv_low n d; high = Fpu.fdiv_high n d }

let pi_I = { low = 0x1.921fb54442d18p+1; high = 0x1.921fb54442d19p+1 }
let two_pi_I = 2.0 *.$ pi_I

let gamma_shift_count (z : float) (t : float) =
  let d = t -. z in
  if d <= 0.0 then 0
  else int_of_float (ceil d)

let lgamma_core (zi : interval) =
  let lnz = log_I zi in
  let base = ((zi -$. 0.5) *$ lnz) -$ zi +$ (0.5 *.$ (log_I two_pi_I)) in
  let sum = ref zero_I in
  Array.iteri (fun i (bn, bd) ->
    let n = i + 1 in
    let k = 2 * n in
    let denom = float_of_int (k * (k - 1)) in
    let zpow = pow_I_i zi (k - 1) in
    let term = (rat_I bn bd) /$ (denom *.$ zpow) in
    sum := !sum +$ term
  ) gamma_bernoulli;
  let result = base +$ !sum in
  let (pn, pd) = gamma_b14 in
  let k = 14 in
  let denom = float_of_int (k * (k - 1)) in
  let zpow = pow_I_i zi (k - 1) in
  let pad_mag = (abs_I ((rat_I pn pd) /$ (denom *.$ zpow))).high in
  result +$ { low = -. pad_mag; high = pad_mag }

let digamma_core (zi : interval) =
  let lnz = log_I zi in
  let base = lnz -$ (inv_I (2.0 *.$ zi)) in
  let sum = ref zero_I in
  Array.iteri (fun i (bn, bd) ->
    let n = i + 1 in
    let k = 2 * n in
    let zpow = pow_I_i zi k in
    let term = (rat_I bn bd) /$ ((float_of_int k) *.$ zpow) in
    sum := !sum +$ term
  ) gamma_bernoulli;
  let result = base -$ !sum in
  let (pn, pd) = gamma_b14 in
  let k = 14 in
  let zpow = pow_I_i zi k in
  let pad_mag = (abs_I ((rat_I pn pd) /$ ((float_of_int k) *.$ zpow))).high in
  result +$ { low = -. pad_mag; high = pad_mag }

let lgamma_scalar_I (z : float) =
  let n = gamma_shift_count z gamma_stirling_T in
  let zi = thin_I z in
  let shifted = zi +$. (float_of_int n) in
  let core = lgamma_core shifted in
  let correction = ref zero_I in
  for i = 0 to n - 1 do
    correction := !correction +$ (log_I (zi +$. (float_of_int i)))
  done;
  core -$ !correction

let digamma_scalar_I (z : float) =
  let n = gamma_shift_count z gamma_stirling_T in
  let zi = thin_I z in
  let shifted = zi +$. (float_of_int n) in
  let core = digamma_core shifted in
  let correction = ref zero_I in
  for i = 0 to n - 1 do
    correction := !correction +$ (inv_I (zi +$. (float_of_int i)))
  done;
  core -$ !correction

let digamma_I (x : interval) =
  if x.low <= 0.0 then failwith "digamma_I: argument must be positive"
  else {
    low = (digamma_scalar_I x.low).low;
    high = (digamma_scalar_I x.high).high;
  }

let gamma_x0_lo = 1.4616321449683622
let gamma_x0_hi = 1.4616321449683625

let lgamma_I (x : interval) =
  if x.low <= 0.0 then failwith "lgamma_I: argument must be positive"
  else if x.high <= gamma_x0_lo then
    { low = (lgamma_scalar_I x.high).low; high = (lgamma_scalar_I x.low).high }
  else if x.low >= gamma_x0_hi then
    { low = (lgamma_scalar_I x.low).low; high = (lgamma_scalar_I x.high).high }
  else
    let min_lo =
      min (lgamma_scalar_I gamma_x0_lo).low (lgamma_scalar_I gamma_x0_hi).low in
    let max_hi =
      max (lgamma_scalar_I x.low).high (lgamma_scalar_I x.high).high in
    { low = min_lo; high = max_hi }

let lgamma (z : float) =
  if z <= 0.0 then nan
  else let iv = lgamma_scalar_I z in 0.5 *. (iv.low +. iv.high)

let digamma (z : float) =
  if z <= 0.0 then nan
  else let iv = digamma_scalar_I z in 0.5 *. (iv.low +. iv.high)
