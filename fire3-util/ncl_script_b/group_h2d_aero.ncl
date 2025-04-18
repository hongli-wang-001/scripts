;----------------------------------------------------------------------
; wrf_gsn_8.ncl
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Using gsn_csm scripts to plot WRF-ARW data
;   - Drawing streamlines colored by another field over a map
;   - Setting the correct WRF map projection using wrf_map_resources
;   - Subsetting a color map
;   - Using stLevelPalette resource to assign a color palette
;   - Using opacity to emphasize or subdue overlain features
;   - Increasing the thickness of map outlines
;----------------------------------------------------------------------
; This script is meant to show the difference between plotting WRF
; data using wrf_xxx scripts, and using gsn_csm_xxx scripts.
;
; The first plot draws the streamlines in a basic lat/lon map
; projection.
;
; The second plot draws the streamlines in the native projection 
; provided on the WRF file.
;----------------------------------------------------------------------
; In NCL Versions 6.3.1 and earlier, you will get these warnings which
; you can safely ignore:
;
;   warning:start_lat is not a valid resource in wrf_gsn_streamline at this time
;   warning:start_lon is not a valid resource in wrf_gsn_streamline at this time
;   warning:end_lat is not a valid resource in wrf_gsn_streamline at this time
;   warning:end_lon is not a valid resource in wrf_gsn_streamline at this time
;   warning:mpNestTime is not a valid resource in map at this time
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;---Open WRF output file

 vNames = (/ "aso4i","aso4j","aso4k","ano3i","ano3j","ano3k","acli", "aclj", "aclk", "anh4i","anh4j","anh4k","anai", "anaj", "amgj", "akj",  "acaj","aeci",       "aecj", "afej", "aalj", "asij","atij","amnj","aothri","aothrj","axyl1j","axyl2j","axyl3j","atol1j","asoil","acors","aseacat","alvpo1i","alvpo1j", "asvpo1i",      "asvpo1j","asvpo2i","asvpo2j","asvpo3j","aivpo1j","alvoo1i","alvoo2i","asvoo1i","asvoo2i" /) ; "T", "P", "U", "V"
 nNames = dimsizes (vNames)   ; Number of variables on the file

 cgroup = (/"DU","SS","OC_Hydrophobic","OC_Hydrophilic","BC_Hydrophobic","BC_Hydrophilic","SU","NI"/)
; igroup = new(45,integer)
 igroup = (/7,7,7,8,8,8,7,2,2,7,7,7,2,2,1,1,1,5,6,1,1,1,1,1,3,3,3,3,3,3,1,1,2,3,3,3,3,3,3,3,3,4,4,4,4/)
 ndim = dimsizes(igroup)

 print(ndim)

  a        = addfile("./fv3_dynvars.nc","r")
  b        = addfile("./fv3_grid_spec.nc","r")
  c        = addfile("./fv3_tracer","r")

  grav = 9.8
  it  = 0
  delp = a->delp
  print(sum(delp(0,:,1,1)))

  tmp = c->aso4i
  printVarSummary(tmp)
  tmp = 0

  ndims = dimsizes(tmp)
  print(ndims)

  ndimx = ndims(3)
  ndimy = ndims(2)

  temp   = a->T(it,0:7,:,:)  

  do n=0, nNames-1              ; Loop through each variable
    va = c->$vNames(n)$        ; Read the variable to memory
    va = va*delp/grav
    print(vNames(n)+" "+igroup(n))
    do i=0,ndimx-1
    do j=0,ndimy-1
    temp(igroup(n)-1,j,i)=temp(igroup(n)-1,j,i) + sum(va(it,:,j,i))
    end do 
    end do 
  end do

  lat = b->grid_latt 
  lon = b->grid_lont 
  printVarSummary(temp)
  printVarSummary(lat)
;---Get the lowest (bottommost) level

;  nl  = 0; 0; 63 
  do nl=0,7
  u10 = temp(nl,:,:)
  v10 = u10
;---Change the metadata
  u10@units       = " "
  u10@long_name = " "
;---Required for plotting over map (not using WRF's map projection here)
  u10@lat2d = lat
  u10@lon2d = lon
  v10@lat2d = lat
  v10@lon2d = lon

  wks = gsn_open_wks("png","aero_h2d_group_"+cgroup(nl))

  res                    = True
  res@gsnMaximize        = True
;  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
;  res@mpCenterLonF      = -106.             ; choose center lon
;  res@mpCenterLatF      = 54.               ; choose center lat
  res@pmTickMarkDisplayMode       = "Always"         ; turn on tickmarks
  res@mpFillOn              = True                    ; no map fill
  res@mpFillColors               = (/"white","lightcyan","lightgray","lightcyan"/); colors to draw

  res@mpMinLatF          = min(lat)-1
  res@mpMaxLatF          = max(lat)+1
  res@mpMinLonF          = min(lon)-1
  res@mpMaxLonF          = max(lon)+1

 res@cnFillOn                    = True             ; turn on color 
 res@cnLinesOn                   = False            ; no contour lines
; res@cnFillPalette               = "BlGrYeOrReVi200"; set color map
 res@cnFillPalette = "spread_15lev" ;"testcmap"
  res@stLineThicknessF   = 3.0
;  res@cnLevels             = fspan(-5,5,11) 
;  res@cnLevelSelectionMode = "ManualLevels"	; manually set the contour levels with the following 3 resources
;  res@cnMinLevelValF  = -1.			; set the minimum contour level
;  res@cnMaxLevelValF  =  1.			; set the maximum contour level
;  res@cnLevelSpacingF = .1			; set the interval between contours

  res@tiMainString       = "Vertically Integrated: "+cgroup(nl)
  res@tiMainFontHeightF  = 0.015
  res@mpDataBaseVersion  = "MediumRes"    ; better map outlines

  res@gsnAddCyclic       = False          ; don't add longitude cyclic point
;
; We like the "NCV_rainbow2" colormap, but don't want to use the 
; whole thing. Use read_colormap_file to read the colormap as 
; an N x 4 array, and then subscript as desired.  Here we are
; starting at color 19 and ending at color 240. Uncomment the
; draw_color_palette call if you want to see what this colormap
; looks like. 
;
;  draw_color_palette(wks,"NCV_rainbow2",0)

;  colormap = read_colormap_file("NCV_rainbow2")
;  res@stLevelPalette = colormap(19:240,:)

;---Customize the some resources so we can see map outlines better
  res@stLineOpacityF              = 0.5       ; make streamlines more transparent
  res@mpUSStateLineColor          = "black"
  res@mpNationalLineColor         = "black"
  res@mpGeophysicalLineColor      = "black"
  res@mpUSStateLineThicknessF     = 2.        ; default is 1.0
  res@mpNationalLineThicknessF    = 2.
  res@mpGeophysicalLineThicknessF = 2.
  plot = gsn_csm_contour_map(wks,u10,res)

 end do ; nl
end

