pro np

@tplot_com.pro

str_element,tplot_vars,'options.trange',trange
if n_elements(trange) ne 2 then timespan

aa = tplot_vars.options

curt = aa.trange[0]
curorb = round(mvn_orbit_num(t=curt))
newt = mvn_orbit_num(orb=curorb+1)


;;; 15 min on either side is good

buffer = 15.*60.

starttime = newt - buffer
endtime = newt + buffer

tlimit, time_string([starttime, endtime])


end