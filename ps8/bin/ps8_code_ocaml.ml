(* Bohdan Purtell
   University of Florida
   Problem Set 8 - Code Section

   AI Disclaimer:
   Artificial intelligence tools were used in the construction of these code
   sections, namely in organizational factors and translation to OCaml.
   Supporting calculations were made on paper originally, and have been ported
   into the comments of this code.

   A supporting appendix of the paper calculations has been included.
*)

open Ps8_helpers
open Numerics

(* GIVEN *)
let filter_order = 4
let cutoff_freq = 100.0 (* Hz *)
let gain = 10.0

(* PART A
   f_c = w_c / (2 pi) = 100 Hz
   => 2 pi (100) = 628.32 rad/s
*)

(* Butterworth
   Butterworth General Form (4th Order)
   let x = s / w_c in
   D_B(x) = x^4 + 2.61312593x^3 + 3.41421356x^2 + 2.61312593x + 1

   and with a constant gain term 10 w_c^4
   H_B(s) =
*)

(* Problem 1, parts B-D. Uses owl-base and the Gnuplot executable.
   Run from the project root: dune exec ./bin/ps8_code_ocaml.exe
   By default, plots are saved in the current working directory. An existing output
   directory can optionally be passed as the first command-line argument.
   Open the resulting SVG files in a browser; no Python or GUI is invoked.
*)

let pi = Float.pi
let wc = 2.0 *. pi *. cutoff_freq
type filter = {
  name : string;
  color : string;
  denominator : float array; (* Polynomial in x = s / wc. *)
}

let response filter omega =
  let a = filter.denominator in
  Complex.div (complex (gain *. a.(filter_order)))
    (polynomial a { Complex.re = 0.0; im = omega /. wc })

(* Part A support: exact fourth-order Butterworth coefficients. *)
let butterworth =
  let a3 = sqrt (4.0 +. 2.0 *. sqrt 2.0) in
  { name = "Butterworth"; color = "#2563eb";
    denominator = [| 1.0; a3; 2.0 +. sqrt 2.0; a3; 1.0 |] }

(* Reverse Bessel polynomial: theta_4(z) = z^4+10z^3+45z^2+105z+105.
   Find q such that |105/theta_4(jq)| = 1/sqrt(2), then substitute
   z = q*s/wc and make the denominator monic. This matches SciPy's
   bessel(4, wc, analog=True, norm="mag"), not delay/phase normalization. *)
let bessel =
  let theta = [| 1.0; 10.0; 45.0; 105.0; 105.0 |] in
   let q = Owl_base.Root.bisec ~xtol:1e-14 (fun q ->
     105.0 /. Complex.norm (polynomial theta { Complex.re = 0.0; im = q })
     -. 1.0 /. sqrt 2.0) 0.0 10.0 in
  { name = "Bessel"; color = "#ea580c";
    denominator = Array.mapi (fun i a -> a /. q ** float_of_int i) theta }

let filters = [ butterworth; bessel ]
let samples = 2000

let frequencies = Array.map (fun x -> 10.0 ** x) (linspace 0.0 5.0 samples)
let angular_frequencies = Array.map (fun f -> 2.0 *. pi *. f) frequencies

(* Part D: solve the unit-step ODE in dimensionless time u = wc*t:
   y'''' + a3*y''' + a2*y'' + a1*y' + a0*y = gain*a0.
   State = [y; dy/du; d2y/du2; d3y/du3], initially zero.
   RK4 replaces scipy.signal.step's matrix-exponential solver. *)
let times = linspace 0.0 0.04 samples
let step_response filter =
  let a = filter.denominator in
  let derivative y =
    [| y.(1); y.(2); y.(3);
       gain *. a.(4) -. a.(1) *. y.(3) -. a.(2) *. y.(2)
       -. a.(3) *. y.(1) -. a.(4) *. y.(0) |]
  in
  let state = ref (Array.make 4 0.0) in
  Array.mapi (fun i t ->
    if i > 0 then begin
      let h = wc *. (t -. times.(i - 1)) in
      state := rk4_step derivative !state h
    end;
    (!state).(0)) times

(* Adapt a filter's label/color to the plotting module's generic curve format. *)
let curve filter xs ys =
  ({ Plot.name = filter.name; color = filter.color }, xs, ys)

let () =
  let directory =
    if Array.length Sys.argv > 1 then Sys.argv.(1)
    else "."
  in
  List.iter (fun filter ->
    let a = filter.denominator in
    Printf.printf "\n%s\nNumerator: [%.12g]\nDenominator: ["
      filter.name (gain *. a.(4) *. wc ** 4.0);
    Array.iteri (fun i v -> Printf.printf "%s%.12g"
      (if i = 0 then "" else "; ") (v *. wc ** float_of_int i)) a;
    Printf.printf "]\n") filters;
  let frequency_data = List.map (fun filter ->
    let h = Array.map (response filter) angular_frequencies in
    let magnitude = Array.map (fun z -> 20.0 *. log10 (Complex.norm z)) h in
    let phase = unwrap (Array.map Complex.arg h) in
    let degrees = Array.map (fun p -> p *. 180.0 /. pi) phase in
    let delay_ms = Array.map (fun d -> -.d *. 1000.0)
        (gradient phase angular_frequencies) in
    (filter, magnitude, degrees, delay_ms)) filters
  in

  (* Part B *)
  Plot.save_svg directory "ps8_bode.svg" ~rows:2 (fun out ->
    let panel = Plot.plot_panel out ~logx:true ~xmin:1.0 ~xmax:1e5
        ~xlabel:"Frequency (Hz)" ~reference:(Some (`Vertical cutoff_freq)) in
    panel ~title:"Fourth-Order LPF Bode Plots" ~ylabel:"Magnitude (dB)"
      ~ymin:(-120.0) ~ymax:40.0 ~yticks:[-120.; -80.; -40.; 0.; 40.]
      (List.map (fun (f, m, _, _) -> curve f frequencies m) frequency_data);
    panel ~title:"Phase Response" ~ylabel:"Phase (degrees)"
      ~ymin:(-360.0) ~ymax:0.0 ~yticks:[-360.; -270.; -180.; -90.; 0.]
      (List.map (fun (f, _, p, _) -> curve f frequencies p) frequency_data));

  (* Part C *)
  let max_delay = List.fold_left (fun m (_, _, _, d) ->
    Array.fold_left max m d) 0.0 frequency_data in
  let delay_limit = ceil (1.1 *. max_delay) in
  Plot.save_svg directory "ps8_group_delay.svg" ~rows:1 (fun out ->
    Plot.plot_panel out ~title:"Group Delay Comparison"
      ~ylabel:"Group Delay (ms)" ~xlabel:"Frequency (Hz)"
      ~logx:true ~xmin:1.0 ~xmax:1e5 ~ymin:0.0 ~ymax:delay_limit
      ~yticks:(Array.to_list (linspace 0.0 delay_limit 6)) ~reference:None
      (List.map (fun (f, _, _, d) -> curve f frequencies d) frequency_data));
  (* Part D *)
  Plot.save_svg directory "ps8_step_response.svg" ~rows:1 (fun out ->
    Plot.plot_panel out ~title:"Step Response Comparison (1 V input)"
      ~ylabel:"Output Voltage (V)" ~xlabel:"Time (ms)"
      ~logx:false ~xmin:0.0 ~xmax:40.0 ~ymin:0.0 ~ymax:12.0
      ~yticks:[0.; 2.; 4.; 6.; 8.; 10.; 12.]
      ~reference:(Some (`Horizontal gain))
      (List.map (fun f ->
        curve f (Array.map (fun t -> t *. 1000.0) times) (step_response f)) filters))
