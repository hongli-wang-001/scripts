load "/scratch2/BMC/wrfruc/hwang/wf1/data_monet/WRFUserARW.ncl"

; read obs lat and lon
;  fname  = "AIRNOW_20190801_20190831_b.nc"
;  a      = addfile(fname,"r")

;  lat2    = a->latitude
;  lon2    = a->longitude
;  obs_lat = lat2(0,:)
;  obs_lon = lon2(0,:)
;  printVarSummary(obs_lat)

; read obs value

;  obs_val = asciiread("pm2p5.txt",-1,"float")
;  obs_val@_FillValue= -1.0
;  printVarSummary(obs_val)
;  obs_val_pre = asciiread("pm2p5_pre.txt",-1,"float")
;  obs_val_pre@_FillValue= -1.0
;  printVarSummary(obs_val_pre)

;;;; read ascii pm2.5 converted from bufr files


  nrow = numAsciiRow("pm2p5.txt")
  ncol = numAsciiCol("pm2p5.txt")

  data_tmp = asciiread("pm2p5.txt",(/nrow,ncol/),"float")
  obs_lat = data_tmp(:,0)
  obs_lon = data_tmp(:,1)
  obs_val = data_tmp(:,4)
  obs_val@_FillValue= -1.0
  printVarSummary(data_tmp)
  print(obs_lat)
  print(obs_lon)
  print(obs_val)
  printMinMax (obs_lat, False)
  printMinMax (obs_lon, False) 
  printMinMax (obs_val, False)

  nobs = dimsizes(obs_val)

  N = num(.not.ismissing(obs_val))
  print("Num and Non-Missing: "+nobs+" "+N)
  avg_o = avg(obs_val)
;  print("avgo= "+avg_o)
  std_o = stddev(obs_val)
;  print("std_f= "+std_o)

  a=addfile("./rcmaq.pm25tot_fcst.nc","r")
  wlat=a->lat 
  wlon=a->lon
  pm2p53d=a->PM25_TOT
  printVarSummary(pm2p53d)
  pm2p5=pm2p53d(0,63,:,:)
  print("pm2p5_model: "+avg(pm2p5)+" "+stddev(pm2p5))

  ;print(min(wlat)+","+max(wlat))  ; 7.97334,55.7424
  ;print(min(wlon)+","+max(wlon))  ; 98.3065,177.803
  ;print(wlat(0,0)+ " "+wlon(0,0))
  ;print(wlat(224,291)+ " "+wlon(224,291))
; A sample observation site
  ;olat=22.57418
  ;olon=238.2977
  olat=50.4641
  olon=279.642
; Find the closest point to the site in WRF grids
  obsij=wrf_latlon_to_ij(wlat,wlon,olat,olon)
  print(obsij)    ; => 63,29
; extract model data at the obs stn
;    grid	: 	Lambert Conformal
;    cen_lon	: 	-97.5
;    cen_lat	: 	38.5
;    stdlat1	: 	38.5
;    stdlat2	: 	38.5
;      nx : 393
;      ny : 225
;      lon1 : -121.70231097
;      lat1 : 22.57417972
;      dx : 13000
;      dy : 13000 
   opt = True
   opt@MAP_PROJ          = 1
   opt@TRUELAT1          =  38.5 
   opt@TRUELAT2          =  38.5 
   opt@STAND_LON         =  -97.5 
;   opt@REF_LAT           =  22.57417972 
;   opt@REF_LON           = -121.70231097 
   opt@REF_LAT           = wlat(0,0) ;22.57417972
   opt@REF_LON           = wlon(0,0) ;-121.70231097 
   opt@KNOWNJ            =  1
   opt@KNOWNI            =  1 

   opt@DX                = 13000.
   opt@DY                = 13000.

   loc = wrf_ll_to_ij( olon, olat, opt )
   print(loc)

  nlat = dimsizes(wlat(:,0))
  nlon = dimsizes(wlon(0,:)) 

  print("grid: nlat_nlon=  "+nlat+" "+nlon)

  oflag = obs_val 
  oflag = 0.0 
  iflag = 0.0
  iflag = where(obs_val.le.60.0,1,0)
  fc_val = obs_val 
  fc_val = -1.0 
  vobs = new ( (/6,nobs/), float ) 
  vobs@_FillValue= -9999.0
  vobs = -9999.0
 
  mobs = 0
  do k = 0, nobs-1
;    print("obs= "+obs_val(k)+ " "+obs_lon(k)+" "+obs_lat(k))

    if (iflag(k).eq.1 .and. .not.ismissing(obs_val(k)) .and. .not.ismissing(obs_lon(k)) .and. .not. ismissing(obs_lat(k) )) then
    ;loc  =wrf_latlon_to_ij(wlat,wlon,olat,olon)
    loc = wrf_ll_to_ij( tofloat(obs_lon(k)), tofloat(obs_lat(k)), opt )
;    print("X/Y location is: " + k + " " + obs_lon(k) + " " + obs_lat(k) )
;
; only verify obs inside model model and set oflag to 1
;
    if (loc(0).ge.1.0 .and. loc(0) .le. 1.0*nlon-1. .and. loc(1).ge.1.0 .and. loc(1).le. 1.0*nlat-1.0) then
     mobs=mobs+1
     oflag(k)=1
     ii=floattoint(loc(0))
     jj=floattoint(loc(1))
     print("X/Y location is: " + k + " " + loc(0) + " "+loc(1))
     print("X/Y location is: " + k + " " + ii + " "+jj)
     dx =loc(0)-ii
     dx1=1.0-dx
     dy =loc(1)-jj
     dy1=1.0-dy
;     print("X/Y location is: " + k + " " + dx + " "+dy)
;
;    convert to ncl ij
;
     ii = ii - 1
     jj = jj - 1
     wrf_val = pm2p5(jj,ii)*dx1*dy1 + pm2p5(jj+1,ii+1)*dx*dy + pm2p5(jj,ii+1)*dx*dy1 + pm2p5(jj+1,ii)*dx1*dy
    print("obs= "+obs_val(k)+ " "+wrf_val+" "+obs_lon(k)+" "+obs_lat(k)+" "+ii+" "+jj+" "+ pm2p5(jj,ii)+" "+pm2p5(jj+1,ii+1))
     vobs(0,k) = tofloat(obs_lat(k))
     vobs(1,k) = tofloat(obs_lon(k)) 
     if(oflag(k).ge.1)then
     vobs(2,k) = obs_val(k)
     vobs(5,k) = wrf_val-obs_val(k)
     end if
     vobs(3,k) = wrf_val
     vobs(4,k) = oflag(k) 
;     if (ob_val(k).lt.243.)then
;      print("Attetions obs T2m is " + obs_val(k))
;     end if
     ;dob_val(k) =  wrf_val - obs_val(k)
    end if ; in domain
    end if ; obs missing
  end do
  print("Num of valid obs in domain is: " + mobs )
  printMinMax (vobs(5,:), False)

  avg_o = avg(vobs(2,:))
  avg_f = avg(vobs(3,:))
  print("avgo_f= "+avg_o+" "+avg_f)

  std_o = stddev(vobs(2,:))
  std_f = stddev(vobs(3,:))
  print("stdo_f= "+std_o+" "+std_f)

  print("avefmo= "+avg(vobs(5,:)))
  print("stdfmo= "+stddev(vobs(5,:)))
  print("rmsfmo= "+sqrt(variance(vobs(5,:))))

   wopt = True
   wopt@fout  = "vobs.txt"
   write_matrix (transpose(vobs), "6f13.4" , wopt)
