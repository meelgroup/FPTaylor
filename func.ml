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
(* Rigorous lgamma / digamma / trigamma support                               *)
(*                                                                            *)
(* Domain: only x > 0 is supported (the reflection formula for x <= 0 is not  *)
(* needed for this tool's use cases).                                        *)
(*                                                                            *)
(* Method: argument-shift recurrence to push the argument above a threshold  *)
(* T, then a truncated Stirling/Euler-Maclaurin asymptotic series (Bernoulli  *)
(* numbers B_2 .. B_12), with the returned interval widened by (at least) the *)
(* magnitude of the first omitted term (the B_14 term) to rigorously account  *)
(* for truncation error, on top of the correctly-rounded interval arithmetic  *)
(* used for every step. See REFERENCE.md / task notes for the derivation and  *)
(* numeric validation (values cross-checked against mpmath to ~1e-17 or       *)
(* tighter at the threshold, comfortably below the padding we add).           *)
(* -------------------------------------------------------------------------- *)

(* Threshold above which the asymptotic series is used directly.
   Chosen (and numerically validated against mpmath) so that the truncation
   error after the B_12 term is comfortably below 2^-60 for all z >= T:
   at z = 20 the true truncation error is about 7.8e-20 for lgamma,
   5.1e-20 for digamma, and 3.6e-20 for trigamma -- all well under
   2^-60 = 8.67e-19. (At the previously-suggested T = 12 the error is only
   about 6e-17, i.e. barely below double precision, so T = 20 is used
   instead for a comfortable rigor margin.) *)
let gamma_stirling_T = 20.0

(* Bernoulli numbers B_2, B_4, B_6, B_8, B_10, B_12 as exact (numerator,
   denominator) pairs of doubles (each is exactly representable). Verified
   against standard tables (Abramowitz & Stegun 23.1.2 / NIST DLMF 24.2). *)
let gamma_bernoulli = [|
  (1.0, 6.0);        (* B_2  = 1/6 *)
  (-1.0, 30.0);       (* B_4  = -1/30 *)
  (1.0, 42.0);        (* B_6  = 1/42 *)
  (-1.0, 30.0);       (* B_8  = -1/30 *)
  (5.0, 66.0);        (* B_10 = 5/66 *)
  (-691.0, 2730.0);   (* B_12 = -691/2730 *)
|]

(* B_14 = 7/6, used only to bound/pad the truncation error of the series. *)
let gamma_b14 = (7.0, 6.0)

let thin_I (x : float) = { low = x; high = x }

(* Correctly-rounded interval enclosure of the exact rational n/d
   (n, d are exactly-representable doubles, e.g. small integers). *)
let rat_I (n : float) (d : float) =
  { low = Fpu.fdiv_low n d; high = Fpu.fdiv_high n d }

(* Rigorous enclosure of pi, and of 2*pi (exact doubling, no extra rounding). *)
let pi_I = { low = 0x1.921fb54442d18p+1; high = 0x1.921fb54442d19p+1 }
let two_pi_I = 2.0 *.$ pi_I

(* Number of shift steps n (>= 0) needed so that z + n >= t, for z > 0.
   Since z > 0, n is always <= ceil(t), so this is O(t) in the worst case
   (e.g. at most 20 steps for T = 20), regardless of how small z > 0 is. *)
let gamma_shift_count (z : float) (t : float) =
  let d = t -. z in
  if d <= 0.0 then 0
  else int_of_float (ceil d)

(* Asymptotic series core, valid for zi enclosing some z >= gamma_stirling_T.
   lgamma(z) ~ (z - 1/2) ln z - z + (1/2) ln(2 pi)
               + sum_{n=1}^{6} B_{2n} / ((2n)(2n-1) z^{2n-1})
   padded by the magnitude of the first omitted (B_14, n=7) term. *)
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

(* digamma(z) ~ ln z - 1/(2z) - sum_{n=1}^{6} B_{2n} / (2n z^{2n}) *)
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

(* trigamma(z) ~ 1/z + 1/(2z^2) + sum_{n=1}^{6} B_{2n} / z^{2n+1} *)
let trigamma_core (zi : interval) =
  let base = (inv_I zi) +$ (0.5 *.$ (inv_I (pow_I_i zi 2))) in
  let sum = ref zero_I in
  Array.iteri (fun i (bn, bd) ->
    let n = i + 1 in
    let k = 2 * n + 1 in
    let zpow = pow_I_i zi k in
    let term = (rat_I bn bd) /$ zpow in
    sum := !sum +$ term
  ) gamma_bernoulli;
  let result = base +$ !sum in
  let (pn, pd) = gamma_b14 in
  let k = 15 in
  let zpow = pow_I_i zi k in
  let pad_mag = (abs_I ((rat_I pn pd) /$ zpow)).high in
  result +$ { low = -. pad_mag; high = pad_mag }

(* Rigorous scalar bound (as a possibly-non-thin interval) of lgamma/digamma/
   trigamma at an exact positive float z, via argument-shift + asymptotic
   series. z is assumed to be an exact real value (a concrete double),
   i.e. thin_I z exactly represents it, so every subsequent step performed
   with correctly-rounded interval arithmetic yields a sound enclosure. *)
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

let trigamma_scalar_I (z : float) =
  let n = gamma_shift_count z gamma_stirling_T in
  let zi = thin_I z in
  let shifted = zi +$. (float_of_int n) in
  let core = trigamma_core shifted in
  let correction = ref zero_I in
  for i = 0 to n - 1 do
    correction := !correction +$ (inv_I (pow_I_i (zi +$. (float_of_int i)) 2))
  done;
  core +$ !correction

(* digamma is strictly increasing on (0, infinity) (its derivative, trigamma,
   is a sum of positive terms), so the interval extension is simply the
   scalar function evaluated (rigorously, rounded outward) at the two
   endpoints. *)
let digamma_I (x : interval) =
  if x.low <= 0.0 then failwith "digamma_I: argument must be positive"
  else {
    low = (digamma_scalar_I x.low).low;
    high = (digamma_scalar_I x.high).high;
  }

(* trigamma is strictly decreasing on (0, infinity) (its derivative is
   -2 * sum 1/(x+n)^3 < 0), so the interval extension swaps the endpoints
   compared to digamma_I. *)
let trigamma_I (x : interval) =
  if x.low <= 0.0 then failwith "trigamma_I: argument must be positive"
  else {
    low = (trigamma_scalar_I x.high).low;
    high = (trigamma_scalar_I x.low).high;
  }

(* lgamma is convex on (0, infinity) with a single interior minimum at
   x0 = 1.46163214496836234126... (where digamma(x0) = 0), and
   lgamma(x0) = -0.12148629053584960824...
   (independently verified via mpmath: mp.findroot(mp.digamma, 1.4616) and
   mp.loggamma of the root, to 60 decimal digits).
   x0_lo/x0_hi below are two consecutive doubles that rigorously bracket the
   true x0 (verified: x0_lo < x0 < x0_hi). *)
let gamma_x0_lo = 1.4616321449683622
let gamma_x0_hi = 1.4616321449683625

let lgamma_I (x : interval) =
  if x.low <= 0.0 then failwith "lgamma_I: argument must be positive"
  else if x.high <= gamma_x0_lo then
    (* Entirely left of the minimum: lgamma is decreasing here. *)
    { low = (lgamma_scalar_I x.high).low; high = (lgamma_scalar_I x.low).high }
  else if x.low >= gamma_x0_hi then
    (* Entirely right of the minimum: lgamma is increasing here. *)
    { low = (lgamma_scalar_I x.low).low; high = (lgamma_scalar_I x.high).high }
  else
    (* The interval straddles (or is close enough to) the minimum: bound the
       low end by evaluating at the (rigorously bracketed) minimizer x0
       itself, and the high end by the larger of the two endpoint values.
       This is always a sound (if occasionally slightly loose) bound, since
       lgamma(x0) <= lgamma(y) for every y > 0. *)
    let min_lo =
      min (lgamma_scalar_I gamma_x0_lo).low (lgamma_scalar_I gamma_x0_hi).low in
    let max_hi =
      max (lgamma_scalar_I x.low).high (lgamma_scalar_I x.high).high in
    { low = min_lo; high = max_hi }

(* Non-rigorous point evaluation, used only for float-valued (display /
   debugging) evaluation; accuracy is good (effectively double precision)
   since it reuses the rigorous interval computation and takes its
   midpoint, but no rounding-direction guarantee is claimed here. *)
let lgamma (z : float) =
  if z <= 0.0 then nan
  else let iv = lgamma_scalar_I z in 0.5 *. (iv.low +. iv.high)

let digamma (z : float) =
  if z <= 0.0 then nan
  else let iv = digamma_scalar_I z in 0.5 *. (iv.low +. iv.high)


