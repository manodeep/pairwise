compile_opt idl2, strictarrsubs

timings_file = 'xx.txt'
output_file = 'timings_with_Numpart.txt'
savefile = str_replace(output_file, '.txt', '.sav')
numpart =  [100, 200L, 500, 1000, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L]
exes = ['./pairwise_3d_histogram', './pairwise_3d']
max_functions = 10
generate_eps = 1

if (findfile(saveFile))[0] eq '' then begin
   if n_elements(alltimings) eq 0 then begin
      get_lun, lun
      alltimings = dblarr(n_elements(exes), n_elements(numpart), max_functions)
      function_names = ptrarr(n_elements(exes))
      num_functions = intarr(n_elements(exes))
      for i = 0, n_elements(exes)-1 do begin
         if ptr_valid(function_names[i]) eq 0 then begin
            for j = 0, n_elements(numpart)-1 do begin
               print, exes[i], numpart[j], format = '("Working on benchmarking ",A0," with ",I07," particles")'
               if j eq 0 then this_functions = []
               
               execstring = exes[i] + " " + string(numpart[j], format = '(I0)') + " 1>" + timings_file
               spawn, execstring
               
               ifunction = 0
               line = ''
               openr, lun, timings_file
               readf, lun, line
               while(strpos(line, "naive") eq -1) do begin
                  readf, lun, line
               endwhile
               xx = strsplit(line, ' ', /extract)
               if j eq 0 then begin
                  this_functions = [this_functions, xx[1]]
               endif
               alltimings[i, j, ifunction] = double(xx[2])
               ifunction++
               while ~EOF(lun) do begin
                  readf, lun, line
                  xx = strsplit(line, ' ', /extract)
                  if j eq 0 then begin
                     this_functions = [this_functions, xx[1]]
                  endif
                  alltimings[i, j, ifunction] = double(xx[2])
                  ifunction++
                  if ifunction gt max_functions then stop
               endwhile
               if j eq 0 then begin
                  num_functions[i] =  n_elements(this_functions)
                  function_names[i] = ptr_new(this_functions, /no_copy)
               endif
               close, lun
            endfor
         endif
      endfor
      spawn, "rm -f " + timings_file
      save, /all, file = savefile
      openw, lun, output_file
      printf, lun, " ", " ", format = '(A30," & ",A30," ",$)'
      for k =  0,  n_elements(numpart)-1 do begin
         printf, lun, numpart[k], format = '(" & ", I12,$)'
      endfor
      printf, lun, " \\"
      
      for i = 0, n_elements(exes)-1 do begin
         this_functions = *function_names[i]
         for j = 0, num_functions[i]-1 do begin
            printf, lun, exes[i], this_functions[j], format = '(A30, " & ", A30," ",$)'
            for k = 0, n_elements(numpart)-1 do begin
               printf, lun, alltimings[i,  k, j], format = '(" & ", F12.3,$)'
      endfor
            printf, lun, " \\"
         endfor
      endfor
      close, lun
      free_lun, lun
   endif 
endif else begin
   restore, savefile
endelse


if generate_eps eq 0 then begin
   set_plot, 'x'
   @idl_reset_graphics
   thick =   3
   xthick =   4
   ythick =   4
   size =   1000
   window,   0,   xsize =   size,   ysize =   size
endif
 
line_colors = cgcolor(['dodger blue', 'BLU6','BLU8', $
                       'RYB3', 'RYB2', 'RYB1', $
                       'Green Yellow', 'Lime Green', 'Olive'])


xrange = minmax(numpart)
xrange[0] *= 0.8
xrange[1] *= 1.2
yrange = [0.0, 100.0]
xticklen = 0.04
yticklen = 0.04
position = [0.2, 0.2, 0.9, 0.9]
symsize = 3
charsize = 4
legend_charsize = 2
xtitle = 'Number of points'
ytitle = '% of Peak FLOPS'
;; titles = ['Binned Pairwise', 'Pairwise']

max_cpu_frequency = 3.2d9 ;;; when cpu hits turbo mode -> max cpu frequency (even though the base is 2.6 GHz)
max_performance = max_cpu_frequency*4 ;; sandy-bridge theoretical max is 4 doubles /cycle (or 8 floats/cycle)
flop_per_pair = [12.0d, 8]
ytickinterval = 25
for i = 0, n_elements(exes)-1 do begin
   if generate_eps eq 1 then begin
      set_plot,   'ps'
      size =   12
      thick =   12
      xthick =   8
      ythick =   8
      device,    decomposed =    1
      !p.font =    1
      psfname = 'scaling_with_numpart_for_'+str_replace(exes[i], './', '')+'.eps'
      device,    filename =    psfname,    /encap,    /color,  /cmyk,  bits_per_pixel =    8,    xsize =    size,    ysize =    size,    /in,    /helvetica,   /tt_font,    /bold,    $
                 xoffset =    0,    yoffset =    0

   endif

   this_functions = *function_names[i]
   if n_elements(line_colors) lt n_elements(this_functions) then stop
   cgplot, [0], /nodata, xrange = xrange, yrange = yrange, $
           /xs, /ys, xticklen = xticklen, yticklen = yticklen, $
           xthick = xthick, ythick = ythick, /xlog, position = position, $
           charsize = charsize, xtitle = xtitle, ytitle = ytitle, $
           ;; title = titles[i], $
           ytickinterval = ytickinterval


   colors = []
   totnflop = (double(numpart))^2*flop_per_pair[i]
   ideal_time_in_ms = totnflop/max_performance*1d3 ;; convert to milli-seconds
   for j = 0, num_functions[i]-1 do begin
      timings = reform(alltimings[i, *, j])
      color = line_colors[j]
      colors =  [colors, color]
      ratio = ideal_time_in_ms/timings*100.0 ;; convert to percent
      if min(ratio) lt 0 then stop
      cgplots, numpart, ratio, color = color, noclip = 0
   endfor
   if i eq 1 then begin
      al_legend, *function_names[i], textcolor = colors, charsize = legend_charsize, box = 0,  /top, /right
   endif else begin
      al_legend, *function_names[i], textcolor = colors, charsize = legend_charsize, box = 0,  /top, /right
   endelse

   cgplot, [0], /nodata, xrange = xrange, yrange = yrange, $
           /xs, /yst, xticklen = xticklen, yticklen = yticklen, $
           xthick = xthick, ythick = ythick, /xlog, position = position, $
           charsize = charsize, xtitle = xtitle, ytitle = ytitle, $
           ;; title = titles[i], $
           /noerase, ytickinterval = ytickinterval

   if generate_eps eq 0 then  begin
      xx = tvread(/true)
      pngfile = 'scaling_with_numpart_for_'+str_replace(exes[i], './', '')+'.png'
      write_png, pngfile, xx
   endif else begin
      device, /close
   endelse
endfor

if generate_eps eq 1 then begin
   set_plot, 'x'
   @idl_reset_graphics
endif

end
