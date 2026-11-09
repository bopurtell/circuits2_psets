# Circuits 2 — Problem Set 8

Problem 1 compares fourth-order Butterworth and Bessel low-pass filters with
a 100 Hz cutoff and DC gain of 10. Bessel normalization places its −3 dB
frequency at 100 Hz, matching SciPy's `norm="mag"`.

Problem 2 designs a Butterworth band-pass filter using a fifth-order low-pass
prototype (a tenth-degree band-pass denominator), with unity center gain,
250 kHz lower cutoff and 1.5 MHz bandwidth. Both programs use OCaml, not MATLAB.

## Layout

```text
circuits2_psets/
└── ps8/
    ├── bin/
     │   ├── ps8_problem1.ml    # Problem 1: filters and parts A–D
     │   ├── ps8_problem2.ml    # Problem 2: band-pass transfer function and Bode plots
    │   └── dune
    ├── lib/
    │   ├── numerics.ml        # Small numerical helpers
    │   ├── plot.ml            # Thin Gnuplot adapter
    │   └── dune              # Ps8_helpers library
    ├── dune-project
    ├── circuits2_ps8.opam # OCaml package dependencies
    └── ps8_problems.pdf
```

## Dependencies

- OCaml 4.12+ and Dune 3.16+.
- `owl-base` 1.2: Owl's pure OCaml subset, used for root finding.
- The **Gnuplot executable**: handles axes, grids, legends and SVG rendering.
- `Complex` and `Unix` come with OCaml. No Python is used.

On Debian/Ubuntu, install Gnuplot with `sudo apt install gnuplot-nox`.
On macOS, use `brew install gnuplot`.

From `circuits2_psets/ps8/`, create a project-local switch to keep
dependencies separate from the existing Opam switch (requires an installed
OCaml compiler):

```sh
opam switch create . ocaml-system
eval $(opam env)
opam install . --deps-only
dune build
dune exec ./bin/ps8_problem1.exe
dune exec ./bin/ps8_problem2.exe
```

For subsequent runs, use either `dune exec` command above from this directory.
If using an already suitable switch instead, skip switch creation and just
install the dependencies. Review Opam's proposed changes before approving;
the existing OxCaml switch on the development machine proposed a large
unrelated rebuild, so testing used a separate temporary switch.

Problem 1 prints transfer-function coefficients and writes
`ps8_bode.svg`, `ps8_group_delay.svg`, and `ps8_step_response.svg`
to the current working directory. Open these files in a browser.
An optional argument selects another existing output directory:

```sh
dune exec ./bin/ps8_problem1.exe -- /path/to/output
dune exec ./bin/ps8_problem2.exe -- /path/to/output
```

If Gnuplot is not on PATH, set `GNUPLOT` to its executable path.
The OCaml `gnuplot` package is not needed: `lib/plot.ml` sends commands directly
to the installed Gnuplot program through `Unix`.

## Problem 2

The second executable writes **`ps8_problem2_bode.svg`**, with stacked magnitude
(dB) and continuous phase (degrees) plots over 1 kHz–100 MHz. The dashed vertical
line marks the center frequency. It also prints the factored transfer function,
without MATLAB or Python.

- Upper cutoff: `fc2 = fc1 + bandwidth = 1.75 MHz`.
- Center: `f0 = sqrt(fc1 * fc2) ≈ 661.438 kHz` (geometric, not arithmetic mean).
- Bandwidth in rad/s: `B = 2π * 1.5e6`.
- Substitute `p = (s² + ω0²)/(B*s)` into
  `B5(p) = (p+1)(p²+0.61803398875p+1)(p²+1.61803398875p+1)`;
  then `H(s) = 1/B5(p)`.

The standard prototype-order convention doubles the polynomial degree under
the band-pass transformation. The phase branch is +450° to −450°, passing
through 0° at the center; adding any multiple of 360° gives equivalent phase.
The assignment requests MATLAB code, but this implementation intentionally
uses OCaml as requested.

## Why a few numerical helpers remain

Full Owl includes array sampling utilities, but also requires native
OpenBLAS/LAPACKE libraries. The lighter `owl-base` package does not expose
`linspace`, so a three-line helper remains. Phase unwrapping, the
NumPy-compatible nonuniform gradient, polynomial evaluation, and RK4 are
also kept together in `lib/numerics.ml`; the filter-specific ODE remains in
`bin/ps8_problem1.ml`. Step responses use RK4 rather than SciPy's matrix-exponential
solver.

Submit the source files, build configuration, and any desired plot outputs.
Do not submit `_build/` or a local `_opam/` switch. Run through Dune rather
than invoking the main `.ml` file alone, so its modules are compiled together.
