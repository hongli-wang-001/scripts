begin


  fname="latlon.tm00"
  loc = asciiread(fname,(/24840,2/),"float") 

  lat  = loc(:,0)
  lon  = loc(:,1)
  print(lat)
  print(lon)
  
  wks = gsn_open_wks("pdf",fname)

  res                             = True
  res@gsnFrame                    = False   ; So we can draw markers
  res@gsnMaximize                 = True    ; affects ps, eps, pdf only
  res@gsnPaperOrientation         = "portrait"  ; force portrait

  res@tiMainString                = "PM 2.5 Observation Distribution 20190806"

  ;res@tiXAxisString    = "Location"
  res@tiXAxisFontHeightF    =0.02

  res@mpFillOn                    = False
  res@mpOutlineDrawOrder          = "PostDraw"
  res@mpFillDrawOrder             = "PreDraw"
  res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
  res@mpUSStateLineColor          = "Gray10"
  res@mpUSStateLineDashPattern    = 2
  res@mpMaxLatF                   = 90           ; choose subregion
  res@mpMinLatF                   = 0
  res@mpMaxLonF                   = 0
  res@mpMinLonF                   = -180


  res@pmTickMarkDisplayMode       = "Always"

; Draw the map (frame won't get advanced because gsnFrame was set to False).
  map = gsn_csm_map(wks,res)


; Draw markers on the plot in the lat/lon locations.
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.01
  mkres@gsMarkerColor = (/"red"/)
  map@locations = gsn_add_polymarker(wks,map,lon,lat,mkres)

  draw(map)
  frame(wks)   ; Now advance the frame.

end
