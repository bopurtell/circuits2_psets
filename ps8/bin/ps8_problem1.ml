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

(* Problem 1, parts B-D. Uses owl-base and the Gnuplot executable.
   Run from the project root: dune exec ./bin/ps8_problem1.exe
   By default, plots are saved in the current working directory. An existing output
   directory can optionally be passed as the first command-line argument.
   Open the resulting SVG files in a browser; no Python or GUI is invoked.
*)

(* Helpers + Constants *)
let pi = Float.pi
let wc = 2.0 *. pi *. cutoff_freq
type filter = {
  name        : string;
  color       : string;
  denominator : float array; (* Polynomial in x = s / wc *)
}

let response filter omega =
  let a = filter.denominator in
  Complex.div (complex (gain *. a.(filter_order)))
    (polynomial a { Complex.re = 0.0; im = omega /. wc })

(* Part A support: exact fourth-order Butterworth coefficients. 
This constructs the raw filter object and it's useful items. *)
let butterworth =
  (* Const: 2.61312593, the cubic term *)
  let a3 = sqrt (4.0 +. 2.0 *. sqrt 2.0) in

  { name = "Butterworth"
  ; color = "#2563eb" (* chosen rendering colour *)
  ; denominator = 
      [| 1.0                (* x^4 *)
      ; a3                  (* x^3 *)
      ; 2.0 +. sqrt 2.0     (* x^2 *)
      ; a3                  (* x^1 *)
      ; 1.0                 (* x^0 *)
      |] 
  }

(* Reverse Bessel polynomial: theta_4(z) = z^4+10z^3+45z^2+105z+105.
   Find q such that |105/theta_4(jq)| = 1/sqrt(2), then substitute
   z = q*s/wc and make the denominator monic. This matches SciPy's
   bessel(4, wc, analog=True, norm="mag"), not delay/phase normalization. *)
let bessel =
  let theta = 
    [| 1.0
    ; 10.0
    ; 45.0
    ; 105.0
    ; 105.0 
    |] in

   let q = Owl_base.Root.bisec ~xtol:1e-14 (fun q ->
     105.0 /. Complex.norm (polynomial theta { Complex.re = 0.0; im = q })
     -. 1.0 /. sqrt 2.0) 0.0 10.0 in

  { name = "Bessel"
  ; color = "#ea580c"
  ; denominator = Array.mapi (fun i a -> a /. q ** float_of_int i) theta 
  }

(* Construct the filters together. *)
let filters : filter list = [ butterworth; bessel ]
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
  (* Output directory: use the command-line argument, or the current directory. *)
  let directory =
    if Array.length Sys.argv > 1 then
      Sys.argv.(1)
    else
      "."
  in

  (* Part A is done on paper and included in the appendix;
     print the coefficients here as a numerical check.
     Convert the normalized polynomial in s/wc into a polynomial in s. *)
  List.iter
    (fun filter ->
      let denominator = filter.denominator in
      let numerator = gain *. denominator.(4) *. wc ** 4.0 in

      Printf.printf "\n%s\nNumerator: [%.12g]\nDenominator: ["
        filter.name
        numerator;

      Array.iteri
        (fun index coefficient ->
          let separator = if index = 0 then "" else "; " in
          let scaled_coefficient = coefficient *. wc ** float_of_int index in

          Printf.printf "%s%.12g" separator scaled_coefficient)
        denominator;

      Printf.printf "]\n")
    filters;

  (* Shared calculations for parts B and C.
     frequency_data stores calculated plot data; it does not draw any plots.
     List.map runs the calculation below once for each filter, producing a
     list with two entries: one for Butterworth and one for Bessel.

     Each entry is a tuple:
       (filter, magnitude (dB), phase (degrees), delay (ms))

     The three result arrays each contain 2000 values. Their positions match
     the frequency arrays: result.(i) belongs to frequencies.(i) in Hz,
     or equivalently angular_frequencies.(i) in rad/s.
     Parts B and C reuse these results instead of recalculating them. *)
  let frequency_data =
    List.map
      (fun filter ->
        (* Evaluate H(jw) at every sampled angular frequency.
           Each result is complex: it contains both magnitude and phase.
           Array.map processes every frequency for this one filter. *)
        let frequency_response =
          Array.map (response filter) angular_frequencies
        in

        (* Complex.norm gives |H(jw)|. Convert this amplitude ratio into dB
           using 20 * log10(|H|); a gain of 10 corresponds to 20 dB. *)
        let magnitude_db =
          Array.map
            (fun value -> 20.0 *. log10 (Complex.norm value))
            frequency_response
        in

        (* Complex.arg extracts the phase in radians, wrapped at +/-pi.
           unwrap removes artificial 2*pi jumps so the phase curve is
           continuous and can be differentiated for group delay. *)
        let phase_radians =
          unwrap (Array.map Complex.arg frequency_response)
        in

        (* Degrees are easier to read on the Bode phase plot.
           Keep phase_radians separately for the group-delay calculation. *)
        let phase_degrees =
          Array.map
            (fun phase -> phase *. 180.0 /. pi)
            phase_radians
        in

        (* Group delay = -d(phase)/d(angular frequency); convert s to ms. *)
        (* gradient estimates the slope between the sampled points, using
           their actual spacing: the frequencies are logarithmically spaced,
           not separated by a constant number of rad/s. *)
        let phase_derivative = gradient phase_radians angular_frequencies in
        let group_delay_ms =
          Array.map (fun slope -> -.slope *. 1000.0) phase_derivative
        in

        (* Bundle this filter and its arrays into one list entry.
           Keeping the filter here preserves its name/color for plotting. *)
        (filter, magnitude_db, phase_degrees, group_delay_ms))
      filters
  in

  (* Part B
     Bode magnitude and phase plots of the frequency responses. *)
  (* Extract the magnitude array from each stored tuple.
     Each underscore ignores a field we do not need for this plot.
     curve combines the filter's name/color, x values (Hz), and y values (dB).
     This prepares plotting inputs; Plot.save_svg below performs the plotting. *)
  let magnitude_curves =
    List.map
      (fun (filter, magnitude_db, _, _) ->
        curve filter frequencies magnitude_db)
      frequency_data
  in

  (* Same frequency x values, but select phase in degrees as the y values. *)
  let phase_curves =
    List.map
      (fun (filter, _, phase_degrees, _) ->
        curve filter frequencies phase_degrees)
      frequency_data
  in

  (* Save one SVG containing two stacked panels: magnitude, then phase.
     out is the command channel to Gnuplot, not a numerical result array. *)
  Plot.save_svg directory "ps8_bode.svg" ~rows:2
    (fun out ->
      (* Both panels use the same logarithmic frequency axis and cutoff marker. *)
      let panel =
        Plot.plot_panel out
          ~logx:true
          ~xmin:1.0
          ~xmax:1e5
          ~xlabel:"Frequency (Hz)"
          ~reference:(Some (`Vertical cutoff_freq))
      in

      panel
        ~title:"Fourth-Order LPF Bode Plots"
        ~ylabel:"Magnitude (dB)"
        ~ymin:(-120.0)
        ~ymax:40.0
        ~yticks:[-120.; -80.; -40.; 0.; 40.] (* . denotes floats explicitly *)
        magnitude_curves;

      panel
        ~title:"Phase Response"
        ~ylabel:"Phase (degrees)"
        ~ymin:(-360.0)
        ~ymax:0.0
        ~yticks:[-360.; -270.; -180.; -90.; 0.]
        phase_curves);

  (* Part C
     Plot group delays, with 10 percent headroom above the largest delay. *)
  (* The inner fold finds the largest delay within one filter's array.
     The outer fold carries that maximum across both filters. This only
     chooses the vertical-axis limit; it does not change the delay data. *)
  let max_delay =
    List.fold_left
      (fun largest_delay (_, _, _, group_delay_ms) ->
        Array.fold_left max largest_delay group_delay_ms)
      0.0
      frequency_data
  in

  (* Add headroom, round the upper limit up to a whole millisecond,
     and place six evenly spaced tick marks on the vertical axis. *)
  let delay_limit = ceil (1.1 *. max_delay) in
  let delay_ticks = Array.to_list (linspace 0.0 delay_limit 6) in

  (* Reuse the fourth field of each frequency_data tuple: group delay in ms. *)
  let delay_curves =
    List.map
      (fun (filter, _, _, group_delay_ms) ->
        curve filter frequencies group_delay_ms)
      frequency_data
  in

  Plot.save_svg directory "ps8_group_delay.svg" ~rows:1
    (fun out ->
      Plot.plot_panel out
        ~title:"Group Delay Comparison"
        ~ylabel:"Group Delay (ms)"
        ~xlabel:"Frequency (Hz)"
        ~logx:true
        ~xmin:1.0
        ~xmax:1e5
        ~ymin:0.0
        ~ymax:delay_limit
        ~yticks:delay_ticks
        ~reference:None
        delay_curves);

  (* Part D
     Plot the output for a 1 V step input; the final value is the DC gain. *)
  (* Part D uses separate time-domain data, not frequency_data.
     Convert the time samples from seconds to milliseconds for the x axis. *)
  let times_ms = Array.map (fun time -> time *. 1000.0) times in

  (* Calculate one output-voltage array per filter, then pair it with time.
     The horizontal reference line below marks the expected final value. *)
  let step_curves =
    List.map
      (fun filter ->
        let output_voltage = step_response filter in
        curve filter times_ms output_voltage)
      filters
  in

  Plot.save_svg directory "ps8_step_response.svg" ~rows:1
    (fun out ->
      Plot.plot_panel out
        ~title:"Step Response Comparison (1 V input)"
        ~ylabel:"Output Voltage (V)"
        ~xlabel:"Time (ms)"
        ~logx:false
        ~xmin:0.0
        ~xmax:40.0
        ~ymin:0.0
        ~ymax:12.0
        ~yticks:[0.; 2.; 4.; 6.; 8.; 10.; 12.]
        ~reference:(Some (`Horizontal gain))
        step_curves)
