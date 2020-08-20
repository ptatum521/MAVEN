PRO bin_dens_plotter, dens_vals, ysplitter, zsplitter, title

  y = [0]
  z = [0]
  margin = [0.3,0.3,0.3,0.3]    ;[left, bottom, right, top]
  t = plot(y,z, xrange=[-3,3], yrange=[-3,3], title=title, xtitle = '$Y_{MSO}$ , RM', ytitle = '$Z_{MSO}$ , RM',dim = [800,800],  margin=margin, font_size = 15)
  max = text(.81,.69,strtrim(ceil(max(dens_vals)),1), target = t, font_size = 15)
  mid = text(.81,.49,strtrim(ceil(max(dens_vals)/2d),1), target = t, font_size = 15)
  zero = text(.81,.29,'0', target = t, font_size = 15)
  c_bar_title1 = text(.88,.510,'Samples',target = t, font_size = 15)
  c_bar_title2 = text(.885,.485,'Per Bin',target = t, font_size = 15)

  ;Black and White Color Bar
  black_to_white = [[indgen(256)],[indgen(256)],[indgen(256)]]
  
  cytop = 0.7
  cybottom = 0.3
  cstep = (cytop-cybottom)/n_elements(black_to_white[*,0])
  c_outline = polygon([[.749,cytop+.001],[.801,cytop+.001],[.801,cybottom-.001],[.749,cybottom-.001]], color = 'black', target = t, thick = 0)
  reversed_bw = reverse(black_to_white)
  for i = 0, (n_elements(black_to_white[*,0])-1) do begin
    c_coord = [[.75,cytop-i*cstep],[.8,cytop-i*cstep],[.8,cytop-(i+1)*cstep],[.75,cytop-(i+1)*cstep]]
    c = polygon(c_coord, /fill_background, fill_color = [black_to_white[i,0],black_to_white[i,1],black_to_white[i,2]],$
      color = [black_to_white[i,0],black_to_white[i,1],black_to_white[i,2]], thick = 0, target = t)
  endfor
  bin_num = 0
  for i = (n_elements(zsplitter)-1),1,-1 do begin
    for j = (n_elements(ysplitter)-1),1,-1 do begin
      ;                 tl                             tr                                 br                             bl
      coords = [[ysplitter[j],zsplitter[i]],[ysplitter[j-1],zsplitter[i]],[ysplitter[j-1],zsplitter[i-1]],[ysplitter[j],zsplitter[i-1]]]
      if dens_vals[bin_num] ne 0 then begin
        distance_bw = floor((n_elements(black_to_white[*,0])-1)*(double(dens_vals[bin_num])/double(max(dens_vals))))
        sample_bw = reversed_bw[distance_bw,*]
        grid = polygon(coords, /fill_background, fill_color = [sample_bw[0,0],sample_bw[0,1],sample_bw[0,2]] ,$
          color = [sample_bw[0,0],sample_bw[0,1],sample_bw[0,2]],thick = 0, /data, target = t)
      endif 
      bin_num = bin_num +1
    endfor
  endfor
  mars = ELLIPSE(0, 0, /DATA, MAJOR=1, COLOR='black', FILL_BACKGROUND=0, TARGET=t)
  return
END
