begin
;---Open files

  a        = addfile("bg_fv3_tracer.nc","r")
  b        = addfile("fv3_grid_spec.nc","r")
  c        = addfile("fv3_tracer.nc","r")

   dstr  = asciiread("var.list", -1, "string")
  printVarSummary(dstr)
  nvar = dimsizes(dstr)

  lat = b->grid_lat ;(0:960,0:959) 
  lon = b->grid_lon ;(0:960,0:959) 

  nd_thin = 1
  it  = 0
  nl  = 63

  do n=0,nvar-1
  temp   = a->$dstr(n)$(it,:,:,:)    ; 3D U at mass points
  temp2  = c->$dstr(n)$(it,:,:,:)    ; 3D V at mass points  nl  = 63; 0; 63 
  temp   = temp2-temp   

  u10 = temp(nl,::nd_thin,::nd_thin)
  printVarSummary(u10) 
  print(max(ndtooned(u10)))
  print(min(ndtooned(u10)))
  xmax=max(ndtooned(abs(u10)))
  print("xmax= "+xmax)
  xmax=(1.0*toint(10000*xmax))/10000.0
  if(xmax.lt.0.0001)then
  xmax=0.0001
  end if



  u10@units       = "ug/kg"
  u10@lat2d = lat(::nd_thin,::nd_thin)
  u10@lon2d = lon(::nd_thin,::nd_thin)

  wks = gsn_open_wks("png","png/diff_"+dstr(n)+"llmp_lev"+nl)

  res                    = True
  res@gsnMaximize        = True
;  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
;  res@mpCenterLonF      = 260.             ; choose center lon
;  res@mpCenterLatF      = 40.               ; choose center lat
  res@pmTickMarkDisplayMode       = "Always"         ; turn on tickmarks
  res@mpFillOn              = True                    ; no map fill
  res@mpFillColors               = (/"white","lightcyan","lightgray","lightcyan"/); colors to draw

  res@mpMinLatF          = 22.0
  res@mpMaxLatF          = 52
  res@mpMinLonF          = 225
  res@mpMaxLonF          = 300

;  res@mpMinLatF          = min(lat)-1
;  res@mpMaxLatF          = max(lat)+1
;  res@mpMinLonF          = min(lon)-1
;  res@mpMaxLonF          = max(lon)+1

 res@cnFillOn                    = True             ; turn on color 
 res@cnLinesOn                   = False            ; no contour lines
; res@cnFillPalette               = "BlGrYeOrReVi200"; set color map
 res@cnFillPalette = "testcmap"
  res@stLineThicknessF   = 3.0
;  res@cnLevels             = fspan(-5,5,11) 
;  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;  res@cnMinLevelValF  = -1.			; set the minimum contour level
;  res@cnMaxLevelValF  =  1.			; set the maximum contour level
;  res@cnLevelSpacingF = .1			; set the interval between contours
  res@cnLevelSelectionMode = "ManualLevels"     ; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = -1*xmax                 ; set the minimum contour level
  res@cnMaxLevelValF  =  1*xmax                 ; set the maximum contour level
  res@cnLevelSpacingF = 2.0*xmax/21    
;  res@cnLevelSelectionMode = "ManualLevels"  
;  res@cnMinLevelValF  = -0.5                    ; set the minimum contour level
;  res@cnMaxLevelValF  =  0.5                    ; set the maximum contour level
;  res@cnLevelSpacingF = .1       
  nl2=nl+1
  res@tiMainString       = "DIFF at Level: "+nl2
  res@tiMainFontHeightF  = 0.015
  res@mpDataBaseVersion  = "HighRes" ;"MediumRes"    ; better map outlines
  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets    = "GeophysicalAndUSStates" ;"National" ;"AllBoundaries"
  res@mpDataBaseVersion        = "MediumRes"
;  res@mpDataSetName            = "Earth..4"      ; U.S. counties

  res@gsnAddCyclic       = False          ; don't add longitude cyclic point
  res@stLineOpacityF              = 0.5       ; make streamlines more transparent
  res@mpUSStateLineColor          = "black"
  res@mpNationalLineColor         = "black"
  res@mpGeophysicalLineColor      = "black"
  res@mpUSStateLineThicknessF     = 2.        ; default is 1.0
  res@mpNationalLineThicknessF    = 2.
  res@mpGeophysicalLineThicknessF = 2.
  plot = gsn_csm_contour_map(wks,u10(::1,::1),res)
  end do
end

