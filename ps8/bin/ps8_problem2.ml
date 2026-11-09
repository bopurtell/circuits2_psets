(* Bohdan Purtell
   University of Florida
   Problem Set 8 - Code Section

   AI Disclaimer:
   Artificial intelligence tools were used in the construction of these code
   sections, including organization, implementation, and supporting calculations.
   Supporting calculations are included in the comments of this code.
*)

open Ps8_helpers
open Numerics

(* GIVEN *)
let prototype_order = 5
let lower_cutoff_hz = 250e3 (* Hz *)
let bandwidth_hz = 1.5e6    (* Hz *)
let gain = 1.0

(* Problem 2, parts A-D. Uses the Gnuplot executable.
   Run from the project root: dune exec ./bin/ps8_problem2.exe
   By default, plots are saved in the current working directory. An existing output
   directory can optionally be passed as the first command-line argument.
   Open the resulting SVG file in a browser; no MATLAB or Python is invoked.

   Order convention: n = 5 is the LOW-PASS PROTOTYPE order. The standard
   low-pass to band-pass substitution produces a degree-10 denominator
   (five band-pass sections), as in butter(5, [w1 w2], 's').
*)

(* Helpers + Constants *)
let pi = Float.pi
let samples = 3000

(* Part C: the bandwidth is the difference between the two cutoff frequencies.
   Part D: the center frequency is their geometric mean. *)
let upper_cutoff_hz = lower_cutoff_hz +. bandwidth_hz
let center_hz = sqrt (lower_cutoff_hz *. upper_cutoff_hz)

(* Convert Hz to rad/s for the transfer function. *)
let bandwidth_rad = 2.0 *. pi *. bandwidth_hz
let omega0 = 2.0 *. pi *. center_hz

(* Part A support: normalized Butterworth polynomial.
   B5(p) = (p+1)(p^2 + a*p + 1)(p^2 + b*p + 1),
   a = (sqrt(5)-1)/2, b = (sqrt(5)+1)/2.

   H(s) = gain / B5((s^2 + omega0^2)/(bandwidth_rad*s)).
   For a polynomial transfer function, let A = s^2 + omega0^2,
   D = bandwidth_rad*s. Then
   H(s) = gain*D^5 / [(A+D)(A^2+a*A*D+D^2)(A^2+b*A*D+D^2)].
   This factored form avoids unwieldy dimensional coefficients. *)
let a = (sqrt 5.0 -. 1.0) /. 2.0
let b = (sqrt 5.0 +. 1.0) /. 2.0

(* Evaluate the frequency response H(jw) for one frequency in Hz. *)
let response frequency_hz =
  (* At s=j*w the substitution is p=j*Omega,
     Omega = (w^2-omega0^2)/(bandwidth_rad*w).
     Work in Hz here: the factors of 2*pi cancel. *)
  let omega =
    (frequency_hz -. center_hz *. (center_hz /. frequency_hz))
    /. bandwidth_hz
  in

  let p = { Complex.re = 0.0; im = omega } in

  let quadratic coefficient =
    polynomial [| 1.0; coefficient; 1.0 |] p
  in

  let denominator =
    Complex.mul (Complex.add p Complex.one)
      (Complex.mul (quadratic a) (quadratic b))
  in

  Complex.div (complex gain) denominator

let () =
  (* Output directory: use the command-line argument, or the current directory. *)
  let directory =
    if Array.length Sys.argv > 1 then
      Sys.argv.(1)
    else
      "."
  in

  (* Parts C and D: print the cutoff and center frequencies. *)
  Printf.printf "PS8 Problem 2: fifth-order Butterworth prototype BPF\n";
  Printf.printf "Band-pass denominator degree: %d\n" (2 * prototype_order);
  Printf.printf "fc1 = %.9g Hz\nBandwidth = %.9g Hz\nfc2 = %.9g Hz\nf0 = %.9g Hz\n"
    lower_cutoff_hz bandwidth_hz upper_cutoff_hz center_hz;

  (* Part A: print the transfer function in factored form. *)
  Printf.printf "\nH(s) = %.9g D^5 / [(A+D)(A^2+%.12g A D+D^2)(A^2+%.12g A D+D^2)]\n"
    gain a b;

  Printf.printf "A = s^2 + %.12g; D = %.12g s (s in rad/s)\n"
    (omega0 *. omega0) bandwidth_rad;

  (* Part B
     Bode magnitude and phase plots of the frequency response. *)
  (* Sample logarithmically from 1 kHz to 100 MHz, including both cutoffs
     and the center exactly in the frequency grid. *)
  let frequencies =
    Array.append
      (Array.map (fun exponent -> 10.0 ** exponent) (linspace 3.0 8.0 samples))
      [| lower_cutoff_hz; center_hz; upper_cutoff_hz |]
  in
  Array.sort Float.compare frequencies;

  (* Each complex response contains the magnitude and phase at one frequency. *)
  let responses = Array.map response frequencies in

  (* Convert the amplitude ratio to dB: unity gain corresponds to 0 dB. *)
  let magnitude_db =
    Array.map
      (fun value -> 20.0 *. log10 (Complex.norm value))
      responses
  in

  (* Remove artificial 2*pi jumps before converting radians to degrees. *)
  let phase_radians = unwrap (Array.map Complex.arg responses) in

  (* Choose the continuous branch +450 -> 0 -> -450 degrees, with zero phase
     at f0. Phase is otherwise defined only modulo 360 degrees. *)
  let phase_degrees =
    Array.map
      (fun phase -> phase *. 180.0 /. pi +. 360.0)
      phase_radians
  in

  (* Both panels use the same filter label and rendering colour. *)
  let style =
    { Plot.name = "Butterworth BPF (n=5 prototype)"
    ; color     = "#2563eb"
    }
  in

  let magnitude_curves = [style, frequencies, magnitude_db] in
  let phase_curves = [style, frequencies, phase_degrees] in

  (* Save one SVG containing two stacked panels: magnitude, then phase.
     Both panels share a logarithmic frequency axis and center-frequency marker. *)
  Plot.save_svg directory "ps8_problem2_bode.svg" ~rows:2
    (fun out ->
      let panel =
        Plot.plot_panel out
          ~logx:true
          ~xmin:1e3
          ~xmax:1e8
          ~xlabel:"Frequency (Hz)"
          ~reference:(Some (`Vertical center_hz))
      in

      panel
        ~title:"PS8 Problem 2: BPF (fc1=250 kHz, fc2=1.75 MHz; dashed line=f0)"
        ~ylabel:"Magnitude (dB)"
        ~ymin:(-160.0)
        ~ymax:5.0
        ~yticks:[-160.; -120.; -80.; -40.; -3.0103; 0.]
        magnitude_curves;

      panel
        ~title:"Continuous Phase (0 degrees at center frequency)"
        ~ylabel:"Phase (degrees)"
        ~ymin:(-450.0)
        ~ymax:450.0
        ~yticks:[-450.; -360.; -270.; -180.; -90.; 0.; 90.; 180.; 270.; 360.; 450.]
        phase_curves)

(* MATLAB version
   Requires the Signal Processing Toolbox for butter and freqs.
   Copy the code below into a MATLAB script to run it.

% GIVEN
n = 5;             % Low-pass prototype order
fc1 = 250e3;       % Hz
beta = 1.5e6;      % Hz
K = 1;

% Helpers + Constants
fc2 = fc1 + beta;
f0 = sqrt(fc1 * fc2);

% Part A: analog Butterworth band-pass transfer function
% A fifth-order prototype gives a tenth-degree band-pass denominator.
[num, den] = butter(n, 2*pi*[fc1 fc2], 'bandpass', 's');
num = K * num;

% Part B: frequency response and Bode plots
f = sort([logspace(3, 8, 3000), fc1, f0, fc2]);
H = freqs(num, den, 2*pi*f);

magnitude_db = 20 * log10(abs(H));
phase_deg = unwrap(angle(H)) * 180/pi + 360;

figure;

subplot(2, 1, 1);
semilogx(f, magnitude_db, 'LineWidth', 1.5);
grid on;
xline(f0, '--');
xlim([1e3 1e8]);
ylim([-160 5]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('PS8 Problem 2: Butterworth Band-Pass Filter');

subplot(2, 1, 2);
semilogx(f, phase_deg, 'LineWidth', 1.5);
grid on;
xline(f0, '--');
xlim([1e3 1e8]);
ylim([-450 450]);
yticks(-450:90:450);
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Continuous Phase (0 degrees at center frequency)');
*)
