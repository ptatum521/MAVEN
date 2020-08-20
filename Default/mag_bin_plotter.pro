PRO mag_bin_plotter, bin_vals, ysplitter, zsplitter, title 



;  if isa(title, /string) ne 1 then begin
;    print, 'Error: TITLE must be a string'
;    return
;  endif

;  if ((n_elements(ysplitter)-1)*(n_elements(zsplitter)-1) ne n_elements(bin_vals)) then begin
;    print, 'Error: Length of BIN_VALS does not match Y and Z division dimentions'
;    return
;  endif  
  
  blue_to_white = [[indgen(256)],[indgen(256)],[fltarr(256)+255]]
  white_to_red = [[fltarr(256)+255],[indgen(256, start=255, increment=-1)],[indgen(256, start=255,  increment=-1)]]
  color_table_70 = [blue_to_white,white_to_red]
  white_index= floor(n_elements(color_table_70[*,0])/2)


  y = [0]
  z = [0]
  margin = [0.3,0.3,0.3,0.3]    ;[left, bottom, right, top]
  t = plot(y,z, xrange=[-3,3], yrange=[-3,3], title=title, xtitle = '$Y_{MSO}$ , RM', ytitle = '$Z_{MSO}$ , RM',dim = [800,800],  margin=margin, font_size = 15)
  one = text(.81,.7,'1.0', target = t, font_size = 15)
  zero = text(.81,.5,'0.0', target = t, font_size = 15)
  neg_one = text(.81,.29,'-1.0', target = t, font_size = 15)
  c_bar_title1 = text(.88,.525,'$B_{X}$',target = t, font_size = 15)
  c_bar_title2 = text(.88,.52,'__', font_size= 15, font_style = 1,target = t)
  c_bar_title3 = text(.879,.48,'|B|',target = t, font_size = 15)

  ;draw colorbar

  cytop = .7
  cybottom =.3
  cstep = (cytop-cybottom)/n_elements(color_table_70[*,0])
  c_outline = polygon([[.749,cytop+.001],[.801,cytop+.001],[.801,cybottom-.001],[.749,cybottom-.001]], color = 'black', target = t, thick = 0)
  reversed_70 = reverse(color_table_70)
  for i = 0, (n_elements(color_table_70[*,0])-1) do begin
    c_coord = [[.75,cytop-i*cstep],[.8,cytop-i*cstep],[.8,cytop-(i+1)*cstep],[.75,cytop-(i+1)*cstep]]
    c = polygon(c_coord, /fill_background, fill_color = [reversed_70[i,0],reversed_70[i,1],reversed_70[i,2]],$
      color = [reversed_70[i,0],reversed_70[i,1],reversed_70[i,2]], thick = 0,target = t)
  endfor

  bin_num = 0

  for i = (n_elements(zsplitter)-1),1,-1 do begin
    for j = (n_elements(ysplitter)-1),1,-1 do begin
      if bin_vals[bin_num] ne 0 then begin
	      ;                 tl                             tr                                 br                             bl
	      coords = [[ysplitter[j],zsplitter[i]],[ysplitter[j-1],zsplitter[i]],[ysplitter[j-1],zsplitter[i-1]],[ysplitter[j],zsplitter[i-1]]]
	      distance_from_white = round((white_index-1)*bin_vals[bin_num]) ;positive or negative #of indeces to move normalized to colorbar
	      sample_70 = color_table_70[white_index + distance_from_white,*]
	      color = [sample_70[0,0],sample_70[0,1],sample_70[0,2]]

	      grid = polygon(coords, /fill_background, fill_color = color ,$
		color = color,thick = 0, /data, target = t)
	      
      endif
      bin_num = bin_num +1
    endfor
  endfor
  mars = ELLIPSE(0, 0, /DATA, MAJOR=1, COLOR='black', FILL_BACKGROUND=0,  TARGET=t)

return
END
