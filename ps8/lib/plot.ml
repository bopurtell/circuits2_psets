(* Thin Gnuplot adapter: Gnuplot handles axes, legends, grids, and SVG output. *)
type style = { name : string; color : string }

let quote text =
  "\"" ^ String.concat "\\\"" (String.split_on_char '"' text) ^ "\""

let plot_panel out ~title ~ylabel ~xlabel ~logx ~xmin ~xmax
    ~ymin ~ymax ~yticks ~reference curves =
  let emit fmt = Printf.fprintf out fmt in
  emit "unset arrow\nunset logscale\n";
  if logx then emit "set logscale x 10\nset mxtics 10\n"
  else emit "set mxtics default\n";
  emit "set title %s\nset xlabel %s\nset ylabel %s\n"
    (quote title) (quote xlabel) (quote ylabel);
  emit "set xrange [%g:%g]\nset yrange [%g:%g]\n" xmin xmax ymin ymax;
  emit "set ytics (%s)\n"
    (String.concat ", " (List.map (Printf.sprintf "%g") yticks));
  emit "set grid xtics ytics mxtics\nset key top right opaque\n";
  (match reference with
   | Some (`Vertical x) ->
       emit "set arrow from first %g, graph 0 to first %g, graph 1 nohead dt 2 lc rgb 'gray'\n" x x
   | Some (`Horizontal y) ->
       emit "set arrow from graph 0, first %g to graph 1, first %g nohead dt 2 lc rgb 'gray'\n" y y
   | None -> ());
  List.iteri (fun i (_, xs, ys) ->
    emit "$curve%d << EOD\n" i;
    Array.iteri (fun i x -> emit "%.17g %.17g\n" x ys.(i)) xs;
    emit "EOD\n") curves;
  emit "plot %s\n" (String.concat ", " (List.mapi (fun i (style, _, _) ->
    Printf.sprintf "$curve%d using 1:2 with lines lw 2 lc rgb %s title %s"
      i (quote style.color) (quote style.name)) curves))

let save_svg directory filename ~rows draw =
  let path = Filename.concat directory filename in
  let executable = Option.value (Sys.getenv_opt "GNUPLOT") ~default:"gnuplot" in
  let out = Unix.open_process_out (Filename.quote executable) in
  (try
     Printf.fprintf out "set terminal svg size 900,%d noenhanced\nset output %s\nset multiplot layout %d,1\n"
       (rows * 380) (quote path) rows;
     draw out;
     output_string out "unset multiplot\nunset output\n"
   with exn -> ignore (Unix.close_process_out out); raise exn);
  match Unix.close_process_out out with
  | Unix.WEXITED 0 -> Printf.printf "Saved %s\n" path
  | _ -> failwith "Gnuplot failed. Install gnuplot or set GNUPLOT to its executable path."
