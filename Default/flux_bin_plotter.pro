PRO flux_bin_plotter, bin_vals, mini, maxi ,ysplitter, zsplitter, title 

  y = [0]
  z = [0]
  margin = [0.3,0.3,0.3,0.3]    ;[left, bottom, right, top]
  t = plot(y,z, xrange=[-3,3], yrange=[-3,3], title=title, xtitle = '$Y_{MSO}$ , RM',$
  ytitle = '$Z_{MSO}$ , RM',dim = [800,800],  margin=margin, font_size = 15)
  maxe = strtrim(floor(maxi),1)
  mine = strtrim(mini,1)
  mide = strtrim(floor((maxi+mini)/2),1)
  max = text(.81,.69,'$10^{'+maxe+'}$', target=t, font_size=15)
  mid = text(.81,.5,'$10^{'+mide+'}$', target=t, font_size=15)
  min = text(.81,.29,'$10^{'+mine+'}$', target=t, font_size=15)
  c_bar_title1 = text(.88,.5,'Flux!C[$cm^{-2}s^{-1}$]',target=t, font_size=15)


  loadct, 33, rgb_table = cb

  cytop = 0.7
  cybottom = 0.3
  cstep = (cytop-cybottom)/n_elements(cb[*,0])
  c_outline = polygon([[.749,cytop+.001],[.801,cytop+.001],[.801,cybottom-.001],[.749,cybottom-.001]], color = 'black', target = t, thick = 0)
  rcb = reverse(cb)
  for i = 0, (n_elements(cb[*,0])-1) do begin
    c_coord = [[.75,cytop-i*cstep],[.8,cytop-i*cstep],[.8,cytop-(i+1)*cstep],[.75,cytop-(i+1)*cstep]]
    c = polygon(c_coord, /fill_background, fill_color = [rcb[i,0],rcb[i,1],rcb[i,2]],$
      color = [rcb[i,0],rcb[i,1],rcb[i,2]], thick = 0, target = t)
  endfor
  
  bin_num = 0
  for i = (n_elements(zsplitter)-1),1,-1 do begin
    for j = (n_elements(ysplitter)-1),1,-1 do begin
      ;                 tl                             tr                                 br                             bl
      coords = [[ysplitter[j],zsplitter[i]],[ysplitter[j-1],zsplitter[i]],[ysplitter[j-1],zsplitter[i-1]],[ysplitter[j],zsplitter[i-1]]]
      if bin_vals[bin_num] ne 0 then begin
        if alog10(bin_vals[bin_num]) gt maxi then begin
          distance_cb = n_elements(cb[*,0])-1
        endif else begin
          if alog10(bin_vals[bin_num]) lt mini then begin
            distance_cb = 0
          endif else begin
            distance_cb = floor((n_elements(cb[*,0])-1)*(alog10(bin_vals[bin_num])-mini)/(maxi-mini))
          endelse
        endelse  
        sample_cb = reform(cb[distance_cb,*])
        color = [sample_cb[0],sample_cb[1],sample_cb[2]]
        grid = polygon(coords, /fill_background, fill_color=color ,$
          color=color,thick = 0, /data, target = t)
      endif 

      bin_num = bin_num +1
    endfor
  endfor
  mars = ELLIPSE(0, 0, /DATA, MAJOR=1, COLOR='black', FILL_BACKGROUND=0, TARGET=t)
  return 


END
