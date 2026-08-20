open Interval

val opt : ?split_mode:string -> ?geometric_ratio_tol:float ->
  (interval array -> interval) -> interval array -> 
  float -> float -> float -> int ->
  float * float * int

val opt_with_x_tols : ?split_mode:string -> ?geometric_ratio_tol:float ->
  (interval array -> interval) -> interval array ->
  float array -> float -> float -> int ->
  float * float * int
