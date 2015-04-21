compile_opt idl2, strictarrsubs

timings_file = 'xx.txt'
output_file = 'timings_with_Numpart_binned.txt'
numpart =  [100, 200L, 500, 1000, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L]
exes = ['./pairwise_3d_histogram']
max_functions = 10
get_lun, lun

if n_elements(alltimings) eq 0 then begin

   alltimings = dblarr(n_elements(exes), n_elements(numpart), max_functions)
   function_names = ptrarr(n_elements(exes))
   num_functions = intarr(n_elements(exes))
   for i = 0, n_elements(exes)-1 do begin
      if ptr_valid(function_names[i]) eq 0 then begin
         for j = 0, n_elements(numpart)-1 do begin
            print, exes[i], numpart[j], format = '("Working on benchmarking ",A0," with ",I07," particles")'
            if j eq 0 then this_functions = []
            
            execstring = exes[i] + " " + string(numpart[j], format = '(I0)') + " 2>/dev/null 1>" + timings_file
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
endif   

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


restore, 'timings_with_Numpart.sav'
line_colors = [ $
               [15, 22, 122], $
               [61, 177, 250], $
               [48, 113, 8], $
               [241, 148, 38], $
               [230, 46, 37], $
               [120, 24, 19], $
               [178, 49, 179]]

xrange = minmax(numpart)
xrange[0] *= 0.8
xrange[1] *= 1.2
yrange = [0.001, 1d3]
size = 800
thick = 4
xticklen = 0.04
yticklen = 0.04
position = [0.2, 0.2, 0.7, 0.7]
symsize = 3
charsize = 4
legend_charsize = 2
xtitle = 'Number of points'
ytitle = 'Time [milli-seconds]'
titles = ['Binned Pairwise', 'Pairwise']

max_cpu_frequency = 3.2d9 ;;; when cpu hits turbo mode -> max cpu frequency (even though the base is 2.6 GHz)
max_performance = max_cpu_frequency*8 ;; sandy-bridge theoretical max is 8 doubles /cycle
flop_per_pair = [12.0d, 8]
for i = 0, n_elements(exes)-1 do begin
   window, xsize = size, ysize = size
   this_functions = *function_names[i]
   if n_elements(line_colors) lt n_elements(this_functions) then stop
   cgplot, [0], /nodata, xrange = xrange, yrange = yrange, $
           /xs, ystyle = 9, xticklen = xticklen, yticklen = yticklen, $
           xthick = xthick, ythick = ythick, /xlog, /ylog, position = position, $
           charsize = charsize, ytickformat = 'logticks_exp', xtitle = xtitle, ytitle = ytitle, $
           title = titles[i]
   colors = []
   for j = 0, num_functions[i]-1 do begin
      timings = reform(alltimings[i, *, j])
      color = reform(line_colors[*, j])
      color = color[0]+256L*color[1]+256L^2*color[2]
      cgplots, numpart, timings, color = color, thick = thick, line = 0, noclip = 0 , psym = -(15+j), symsize = symsize
      colors = [colors, color]
   endfor
   cgplot, [0], /nodata, xrange = xrange, yrange = yrange, $
           /xs, ystyle = 9, xticklen = xticklen, yticklen = yticklen, $
           xthick = xthick, ythick = ythick, /xlog, /ylog, position = position, $
           charsize = charsize, ytickformat = 'logticks_exp', /noerase

   al_legend, *function_names[i], textcolor = colors, charsize = legend_charsize, /box,  $
           position = [position[2]-0.04, position[3]+0.25], /normal

   ;;; Now draw the alternate x-axis
   axis, /yaxis, yrange = [0.0, 100.0], /save, yticklen = yticklen, ytitle = '% of Peak FLOPS', charsize = charsize, $
         ythick = ythick, ylog = 0, ytickinterval = 50, yminor = 5

   totnflop = numpart^2*flop_per_pair[i]
   ideal_time_in_ms = totnflop/max_performance*1d3 ;; convert to milli-seconds
   for j = 0, num_functions[i]-1 do begin
      timings = reform(alltimings[i, *, j])
      color = reform(line_colors[*, j])
      color = color[0]+256L*color[1]+256L^2*color[2]
      ratio = ideal_time_in_ms/timings*100.0 ;; convert to percent
      if min(ratio) lt 0 then stop
      cgplots, numpart, ratio, color = color, line = 2, noclip = 0
      colors = [colors, color]
   endfor

   xx = tvread(/true)
   pngfile = 'scaling_with_numpart_for_'+str_replace(exes[i], './', '')+'.png'
   write_png, pngfile, xx
   
endfor


end
