pro zop, long=long 
;;zoom_on_peri

;;;; use to zoom in on periapsis pass

ctime,time,tmp1,tmp2,npoints=1,prompt='Select a periapsis pass', /silent

peri = round(mvn_orbit_num(t=time))
peritime = mvn_orbit_num(orb=peri)

;;; 15 min on either side is good

buffer = 15.*60.

if keyword_set(long) then buffer = 30.*60.

starttime = peritime - buffer
endtime = peritime + buffer 

tlimit, time_string([starttime, endtime])


end