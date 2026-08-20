open Interval

let mk_const_interval v = {low = v; high = v}

type dom = {
  bounds : interval array;
  mid : float array;
}

type split_mode =
  | Midpoint
  | Geometric

let split_mode_of_string = function
  | "midpoint" -> Midpoint
  | "geometric" -> Geometric
  | mode -> failwith ("Unknown bb split mode: " ^ mode)

let arithmetic_midpoint b =
  (b.low +. b.high) *. 0.5

let geometric_span b =
  if b.low > 0.0 && b.high > b.low then
    Some (log (b.high /. b.low))
  else if b.high < 0.0 && b.high > b.low then
    Some (log ((-. b.low) /. (-. b.high)))
  else
    None

let geometric_midpoint b =
  let midpoint_abs lo hi sign =
    let m = sign *. exp ((log lo +. log hi) *. 0.5) in
    if b.low < m && m < b.high then Some m else None in
  if b.low > 0.0 && b.high > b.low then
    midpoint_abs b.low b.high 1.0
  else if b.high < 0.0 && b.high > b.low then
    midpoint_abs (-. b.high) (-. b.low) (-1.0)
  else
    None

let split_point split_mode b =
  match split_mode with
  | Midpoint -> arithmetic_midpoint b
  | Geometric ->
    begin match geometric_midpoint b with
    | Some m -> m
    | None -> arithmetic_midpoint b
    end

let split_score split_mode b =
  match split_mode with
  | Midpoint -> b.high -. b.low
  | Geometric ->
    match geometric_span b with
    | Some s -> s
    | None ->
      b.high -. b.low

let domain_small split_mode x_tols geometric_ratio_tol bounds =
  let ratio_tol = log geometric_ratio_tol in
  let interval_small i b =
    let x_tol = x_tols.(i) in
    match split_mode with
    | Midpoint -> b.high -. b.low <= x_tol
    | Geometric ->
      begin match geometric_span b with
      | Some s -> s <= ratio_tol
      | None -> b.high -. b.low <= x_tol
      end in
  let small = ref true in
  for i = 0 to Array.length bounds - 1 do
    small := !small && interval_small i bounds.(i)
  done;
  !small

let mk_dom split_mode a = {
  bounds = a;
  mid = Array.map (split_point split_mode) a;
}

let split_dom split_mode dom =
  let w = Array.mapi (fun i b -> (i, split_score split_mode b)) dom.bounds in
  let (i, _) = Array.fold_left 
    (fun (i, v1) (j, v2) -> if v1 > v2 then (i, v1) else (j, v2)) (0, neg_infinity) w in
  let bi = dom.bounds.(i) in
  let m = split_point split_mode bi in
  let m1 = split_point split_mode {low = bi.low; high = m} and
      m2 = split_point split_mode {low = m; high = bi.high} in
  let d1 = {
    bounds = (let c = Array.copy dom.bounds in 
		      let _ = c.(i) <- {low = c.(i).low; high = m} in c);
    mid = (let c = Array.copy dom.mid in
	   let _ = c.(i) <- m1 in c);
  } in
  let d2 = {
    bounds = (let c = Array.copy dom.bounds in 
	      let _ = c.(i) <- {low = m; high = c.(i).high} in c);
    mid = (let c = Array.copy dom.mid in
	   let _ = c.(i) <- m2 in c);
  } in
  d1, d2


let opt0 split_mode f x_tols geometric_ratio_tol f_rel_tol f_abs_tol max_iters =
  let counter = ref 0 in
  let rec opt upper_bound lower_bound doms acc =
    match doms with
      | [] -> 
	if acc = [] then
	  upper_bound, lower_bound
	else
	  (* Gives "Stack overflow" error for macro2/delta.txt *)
(*	  let doms0 = List.sort (fun (v1, _) (v2, _) -> compare v2 v1) acc in
	  let doms1 = List.map snd doms0 in
*)
	  let doms0 = Array.of_list acc in
	  let _ = Array.sort (fun (v1, _) (v2, _) -> compare v2 v1) doms0 in
	  let doms1 = Array.to_list (Array.map snd doms0) in
	  opt upper_bound lower_bound doms1 []
      | dom :: rest ->
	let v = f dom.bounds in
	if v.high <= lower_bound then
	  opt upper_bound lower_bound rest acc
	else
	  let d_min = Array.map (fun d -> mk_const_interval d.low) dom.bounds and
	      d_max = Array.map (fun d -> mk_const_interval d.high) dom.bounds and
	      d_mid = Array.map mk_const_interval dom.mid in
	  let v2_min = f d_min and
	      v2_max = f d_max and
	      v2_mid = f d_mid in
	  let v2 = max (max v2_min.low v2_max.low) v2_mid.low in
	  let lower_bound = max v2 lower_bound in
	  if abs_float (v.high -. v2) <= f_rel_tol *. abs_float v2 +. f_abs_tol || 
	    domain_small split_mode x_tols geometric_ratio_tol dom.bounds ||
	    (max_iters >= 0 && !counter >= max_iters) then
	    opt (max upper_bound v.high) lower_bound rest acc
	  else
	    let _ = counter := !counter + 1 in
	    let d1, d2 = split_dom split_mode dom in
	    opt (max upper_bound lower_bound) lower_bound rest ((v2, d1) :: (v2, d2) :: acc)
  in
  fun upper_bound lower_bound doms acc ->
    let _ = counter := 0 in
    let upper_bound, lower_bound = opt upper_bound lower_bound doms acc in
    upper_bound, lower_bound, !counter

let opt_with_x_tols ?(split_mode = "midpoint") ?(geometric_ratio_tol = 2.0)
    f a x_tols f_rel_tol f_abs_tol max_iters =
  let split_mode = split_mode_of_string split_mode in
  if geometric_ratio_tol <= 1.0 then
    failwith "bb geometric ratio tolerance must be greater than 1";
  if Array.length x_tols <> Array.length a then
    failwith "opt_with_x_tols: tolerance vector size does not match domain size";
  let upper_bound, lower_bound, counter = 
    opt0 split_mode f x_tols geometric_ratio_tol f_rel_tol f_abs_tol max_iters
      neg_infinity neg_infinity [mk_dom split_mode a] [] in
  (upper_bound, lower_bound, counter)

let opt ?(split_mode = "midpoint") ?(geometric_ratio_tol = 2.0)
    f a x_tol f_rel_tol f_abs_tol max_iters =
  opt_with_x_tols ~split_mode ~geometric_ratio_tol
    f a (Array.make (Array.length a) x_tol) f_rel_tol f_abs_tol max_iters
