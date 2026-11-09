(* Small helpers not supplied by owl-base; see README.md for the tradeoff.
   Root finding uses owl-base directly in the main file. *)
let complex x = { Complex.re = x; im = 0.0 }

(* Horner evaluation; coefficients are in descending powers. *)
let polynomial coefficients s =
  Array.fold_left
    (fun acc a -> Complex.add (Complex.mul acc s) (complex a))
    Complex.zero coefficients

let linspace start stop count =
  Array.init count (fun i ->
    start +. (stop -. start) *. float_of_int i /. float_of_int (count - 1))

let unwrap angles =
  let offset = ref 0.0 in
  Array.mapi (fun i angle ->
    if i > 0 then begin
      let jump = angle -. angles.(i - 1) in
      if jump > Float.pi then offset := !offset -. 2.0 *. Float.pi
      else if jump < -.Float.pi then offset := !offset +. 2.0 *. Float.pi
    end;
    angle +. !offset) angles

(* Same nonuniform central difference and first-order endpoints as
   np.gradient(y, x), whose default edge_order is 1. *)
let gradient y x =
  let n = Array.length y in
  Array.init n (fun i ->
    if i = 0 then (y.(1) -. y.(0)) /. (x.(1) -. x.(0))
    else if i = n - 1 then
      (y.(i) -. y.(i - 1)) /. (x.(i) -. x.(i - 1))
    else
      let left = x.(i) -. x.(i - 1) in
      let right = x.(i + 1) -. x.(i) in
      -.right /. (left *. (left +. right)) *. y.(i - 1)
      +. (right -. left) /. (left *. right) *. y.(i)
      +. left /. (right *. (left +. right)) *. y.(i + 1))

(* One fourth-order Runge-Kutta step for an autonomous ODE y' = f(y). *)
let rk4_step derivative y h =
  let shift y k h = Array.mapi (fun i v -> v +. h *. k.(i)) y in
  let k1 = derivative y in
  let k2 = derivative (shift y k1 (h /. 2.0)) in
  let k3 = derivative (shift y k2 (h /. 2.0)) in
  let k4 = derivative (shift y k3 h) in
  Array.mapi (fun j v ->
    v +. h /. 6.0 *. (k1.(j) +. 2.0 *. k2.(j)
                        +. 2.0 *. k3.(j) +. k4.(j))) y
