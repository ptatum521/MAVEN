;+
;
;PROCEDURE:       MVN_STA_UPDOWN_DIR
;
;PURPOSE:         Computes the MAVEN STATIC upward and downward ion
;                 population relative to the planet.
;                 The results are returned as the tplot variables. 
;
;INPUTS:          None. 
;                 But APID used in this procedure should be often specified.
;
;KEYWORDS:
;
;  TRANGE:        Specifies the time range to compute.
;
;   ANGLE:        Specifies the cone angle relative to the radial
;                 direction to the planet to determine the upward/downward
;                 ion populations. Default is 75 deg.
;
;    APID:        Specifies the APID data product. Default is 'D0'.
;
;    UNIT:        Specifies the data unit. Default is 'flux'.
;
;    MASS:        Specifies the mass range to compute. Default is all.  
;
;CREATED BY:      Takuya Hara on 2015-08-24.
;
;LAST MODIFICATION:
; $LastChangedBy: hara $
; $LastChangedDate: 2016-07-26 17:02:29 -0700 (Tue, 26 Jul 2016) $
; $LastChangedRevision: 583 $
; $URL: svn+ssh://hara@mojo.ssl.berkeley.edu/home/hara/work/svn/mypro/trunk/mars/maven/sta/mvn_sta_updown_dir.pro $
;
;-
PRO mvn_sta_updown_dir_pt, trange=tspan, angle=angle, heavy=heavy, proton=proton, $
                        verbose=verbose, unit=unit, apid=apid, mass=mass

  nan = !values.f_nan
  dnan = !values.d_nan
  IF SIZE(mass, /type) EQ 0 THEN mflg = 0 ELSE mflg = 1  ;type 0 is 'undefined'
  IF SIZE(heavy, /type) EQ 0 THEN hflg = 0 ELSE hflg = FIX(heavy)
  IF SIZE(proton, /type) EQ 0 THEN pflg = 0 ELSE pflg = FIX(proton)
  IF SIZE(angle, /type) EQ 0 THEN angle = 75.
  IF SIZE(unit, /type) EQ 0 THEN unit = 'flux'
  IF SIZE(apid, /type) EQ 0 THEN apid = 'd0'
  apid = apid[0]
  status = EXECUTE('time = mvn_sta_get_' + apid + '(/times)')
  ;IF status EQ 0 THEN GOTO, no_data
  ;;IF keyword_set(tspan) THEN trange = time_double(tspan) $
  ;;ELSE get_timespan, trange
  ;------------------------------------WHERE RETURNS -1 EVEN THOUGH ITS THERE!================
;  IF KEYWORD_SET(tspan) THEN BEGIN
;     trange = time_double(tspan)
;     w = WHERE(time GE MIN(trange) AND time LE MAX(trange), ndat)
;     IF ndat EQ 0 THEN BEGIN
;        no_data:
;        dprint, 'STATIC ' + STRUPCASE(apid) + ' data is not available in the specified time interval.', dlevel=2, verbose=verbose
;        RETURN
;     ENDIF ELSE time = time[w]
; -------------------------------------------------------------------------------------    
     ;undefine, idx
  ;ENDIF ELSE BEGIN
     ndat = N_ELEMENTS(time)
     w = INDGEN(ndat)
  ;ENDELSE 
  ;--------------------------------------------SPICE BROKEN------------------------------------
 
;  mso = TRANSPOSE(spice_body_pos('MAVEN', 'Mars', utc=time, frame='MAVEN_MSO')) 
;  ;IF status LT 2 THEN BEGIN
;  IF N_ELEMENTS(mso[*, 0]) NE ndat THEN BEGIN
;     dprint, 'SPICE/kernels are insufficient.', dlevel=2, verbose=verbose
;     RETURN
;  ENDIF
  ;----------------------------------------------------------------------------------------------
  
  get_mvn_eph, time, pos, verbose=verbose
  mso = [ [pos.x_ss], [pos.y_ss], [pos.z_ss] ]
  
  
  mso = mso / REBIN(SQRT(TOTAL(mso*mso, 2)), ndat, 3)
  mso_sta = TRANSPOSE(spice_vector_rotate(TRANSPOSE(mso), time, 'MAVEN_MSO', 'MAVEN_STATIC', check_obj='MAVEN_SPACECRAFT', verbose=verbose))

  FOR i=0, ndat-1 DO BEGIN
     ;sta = CALL_FUNCTION('mvn_sta_get_' + apid, time[i])
     sta = CALL_FUNCTION('mvn_sta_get_' + apid, index=w[i])
     sta = conv_units(sta, unit)
     IF i EQ 0 THEN BEGIN
        dflux = FLTARR(ndat, sta.nenergy)
        uflux = dflux
        energy = dflux
        denergy = dflux
        dsang = FLTARR(ndat)
        usang = dsang
        mode = INTARR(ndat)
        tstep = DBLARR(ndat)
     ENDIF 
     ;energy[i, *] = sta.energy[*, 0]
     ;denergy[i, *] = sta.denergy[*, 0]
     IF ndimen(sta.data) EQ 3 THEN BEGIN
        energy[i, *] = average(average(sta.energy, 3), 2) ;averages over mass and angle bins to give 32 energy bins
        denergy[i, *] = average(average(sta.denergy, 3), 2)       
     ENDIF ELSE BEGIN
        energy[i, *] = average(sta.energy, 2)
        denergy[i, *] = average(sta.denergy, 2)
     ENDELSE 
     mode[i] = sta.mode
     tstep[i] = 0.5 * (sta.time + sta.end_time)
     domega = REFORM(((REFORM(sta.domega, [sta.nenergy, sta.nbins, sta.nmass]))[sta.nenergy-1, *, 0])) ;gives 64 domega bins
     scbins = REFORM(sta.data, [sta.nenergy, sta.nbins, sta.nmass])
     scbins[*] = 1.

     FOR j=0, sta.nbins-1 DO scbins[*, j, *] = sta.bins_sc[j] 
     scbins = REFORM(scbins) ;[32,64,8]

     IF (mflg) THEN BEGIN
        idx = WHERE(sta.mass_arr LT MIN(mass) OR sta.mass_arr GT MAX(mass), nidx)
        IF nidx GT 0 THEN sta.data[idx] = 0.
        undefine, idx, nidx
     ENDIF 
     IF (hflg) THEN BEGIN
        IF ((apid EQ 'd0') OR (apid EQ 'd1')) THEN BEGIN
           sta.data[*, *, 0:3] = 0.
           sta.data[*, *, 7] = 0.
        ENDIF 
        IF (apid EQ 'd4') THEN sta.data[*, *, 0] = 0.
     ENDIF 
     IF (pflg) AND ((apid EQ 'd0') OR (apid EQ 'd1') OR (apid EQ 'd4')) THEN sta.data[*, *, 1:*] = 0.
     
     vector = REFORM(mso_sta[i, *])                                        ;grabs current 3 components as row vector
     sphere_to_cart, 1.d0, sta.theta, sta.phi, x, y, z

     theta = ACOS(x * vector[0] + y * vector[1] + z * vector[2]) * !RADEG
     idx = WHERE(theta LE angle, nidx)                                     ;finding region inside 'cone' 
     IF nidx GT 0 THEN BEGIN
        weight = sta.data * 0.
        weight[idx] = 1.                                                   ;used in integral to selct inside cone
        ibins = WHERE((REFORM(weight, [sta.nenergy, sta.nbins, sta.nmass]))[sta.nenergy-1, *, 0] EQ 1.)
        IF STRLOWCASE(unit) NE 'counts' THEN $
                                                ;data units 'flux'
           uflux[i, *] = TOTAL( TOTAL( ( REFORM(sta.data * sta.domega * weight * scbins * ABS(COS(theta*!DTOR)), [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan) $
                         / TOTAL(TOTAL(( REFORM(sta.domega * weight * scbins, [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan) $
        ELSE uflux[i, *] = TOTAL(TOTAL(( REFORM(sta.data * weight * scbins * ABS(COS(theta*!DTOR)), [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan)
        usang[i] = TOTAL( (domega * sta.bins_sc)[ibins])
     ENDIF ELSE uflux[i, *] = nan
     undefine, idx, nidx, ibins

     idx = WHERE(theta GE 180. - angle, nidx)  ;flipping the angle to 'downward'
     IF nidx GT 0 THEN BEGIN
        weight = sta.data * 0.
        weight[idx] = 1.
        ibins = WHERE((REFORM(weight, [sta.nenergy, sta.nbins, sta.nmass]))[sta.nenergy-1, *, 0] EQ 1.)
        IF STRLOWCASE(unit) NE 'counts' THEN $
                                                ;data units 'flux'
           dflux[i, *] = TOTAL( TOTAL( ( REFORM(sta.data * sta.domega * weight * scbins * ABS(COS(theta*!DTOR)), [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan) $
                         / TOTAL(TOTAL(( REFORM(sta.domega * weight * scbins, [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan) $
        ELSE dflux[i, *] = TOTAL(TOTAL(( REFORM(sta.data * weight * scbins * ABS(COS(theta*!DTOR)), [sta.nenergy, sta.nbins, sta.nmass]) ), 3, /nan), 2, /nan)
        dsang[i] = TOTAL( (domega * sta.bins_sc)[ibins])
     ENDIF ELSE dflux[i, *] = nan
     undefine, idx, nidx, ibins
     undefine, theta, x, y, z, vector, sta
  ENDFOR 

  cang = 2. * !PI * (1. - COS(angle*!DTOR))
  zrange = minmax([minmax(uflux, /pos), minmax(dflux, /pos)])
  status = EXECUTE("dlim = mvn_sta_etspec_get_dlim(apid, unit=unit)")
  IF STRLOWCASE(apid) EQ 'ca' THEN str_element, dlim, 'datagap', /delete
  IF STRLOWCASE(apid) EQ 'd4' THEN $
     dlim={ytitle: 'STA D4', ylog: 1, yrange: zrange, ystyle: 1, ytickformat: 'exponent', ysubtitle: 'Total Flux'}
  
  store_data, 'mvn_sta_' + apid + '_upward', data={x: tstep, y: uflux, v: energy, dv: denergy, mode: mode}, $
              dlim=dlim, lim={ysubtitle: 'Energy [eV]!CUp'}
  
  store_data, 'mvn_sta_' + apid + '_downward', data={x: tstep, y: dflux, v: energy, dv: denergy, mode: mode}, $
              dlim=dlim, lim={ysubtitle: 'Energy [eV]!CDown'}

  store_data, 'mvn_sta_' + apid + '_upward_fov', data={x: tstep, y: usang/cang}, $
              dlim={ytitle: 'STA', ysubtitle: 'Up!CFOV'}
  store_data, 'mvn_sta_' + apid + '_downward_fov', data={x: tstep, y: dsang/cang}, $
              dlim={ytitle: 'STA', ysubtitle: 'Down!CFOV'}

  IF keyword_set(mass) THEN $
     options, 'mvn_sta_' + apid + ['_upward', '_downward'], mass=minmax(mass)

  RETURN
END
