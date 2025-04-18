
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;  A = addfile("./fv3_tracer","r")

   diri = "./"
   fili = systemfunc("cd "+diri+" ; ls pe*aod*_01.nc4")
   nfili= dimsizes(fili)
   print("nfili="+nfili)
   A    = addfiles(diri+fili,"r")
   printVarSummary(A)
 vNames = (/ "aso4i","aso4j","aso4k","ano3i","ano3j","ano3k","acli", "aclj", "aclk", "anh4i","anh4j","anh4k","anai", "anaj", "amgj", "akj",  "acaj","aeci", "aecj", "afej", "aalj", "asij","atij","amnj","aothri","aothrj","axyl1j","axyl2j","axyl3j","atol1j","asoil","acors","aseacat","alvpo1i","alvpo1j", "asvpo1i","asvpo1j","asvpo2i","asvpo2j","asvpo3j","aivpo1j","alvoo1i","alvoo2i","asvoo1i","asvoo2i" /) ; "T", "P", "U", "V" 
 nNames = dimsizes (vNames)   ; Number of variables on the file

; igroup = new(45,integer)
 igroup = (/7,7,7,8,8,8,7,2,2,7,7,7,2,2,1,1,1,5,6,1,1,1,1,1,3,3,3,3,3,3,1,1,2,3,3,3,3,3,3,3,3,4,4,4,4/)
 ndim = dimsizes(igroup)

 print(ndim)

 ave = new((/8,64/),double)
 ave = 0.0
 do n=0, nNames-1              ; Loop through each variable
    va = A[:]->$vNames(n)$        ; Read the variable to memory
    print(vNames(n)+" "+igroup(n))
    printVarSummary(va)
    do i=0,63
       ave(igroup(n)-1,i)=ave(igroup(n)-1,i)+avg(va(:,i))
    end do
    print(vNames(n)+" "+avg(va))
    ;gName = vNames(n) 
    ;A->$gName$ = va
 end do

   ave!0 = "x"
   ave!1 = "y"
   ave&y = ispan(1,64,1) 

;-- open workstation
   wks_type = "png"                                       ;-- plot output type
   wks  = gsn_open_wks (wks_type, "group_geos5_obs_space" )

;-- resources pressure height plot
  res                     = True          ; plot mods desired

 res@trYReverse        = True

  res@cnFillOn            = True          ; turn on color
;  res@cnFillPalette       = "gui_default" ; set color map
  res@cnFillPalette       = "sunshine_9lev"
  res@cnLinesOn       = False
  res@lbOrientation        = "Vertical"
  res@gsnXAxisIrregular2Log = True

;  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;if (n.eq.0)then
  ;res@cnLevels    = 100.0*ispan(4,56,4)
  ;else
;  res@cnLevels    = 1000.0*ispan(4,56,4) ;(/ -12, -10.,-8.,-6.,-4.,-2.,-1.,1.,2.,4.,6.,8.,10.,12./)   ; set levels
  ;end if
;  res@vpWidthF            = 0.5           ; change aspect ratio of plot
;  res@vpHeightF           = 0.3

; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend
 
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .90                   ; move units right
 res@pmLegendOrthogonalPosF = -0.8                  ; more neg = down
 
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.25                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@trYMaxF = 64
 res@trYMinF = 1 
 res@xyExplicitLegendLabels = (/"DU","SS","OC_Hydrophobic","OC_Hydrophilic","BC_Hydrophobic","BC_Hydrophilic","SU","NI"/)
 res@xyDashPatterns       =  (/0,0,1,2,1,2,3,4/)                  ; dashed line for 2nd
 res@xyLineColors          = (/"orange","skyblue","red","red","black","black","blue","gray"/)  ; 3 different colors
 res@xyLineThicknessF     = 2.0

  res@gsnMaximize         = True          ; maximize plot size
;            description        :       Mass extinction coefficient for aerosol scatterers
;            units      :       Metres squared per kilogram (m^2.kg^-1)
;  res@tiMainString             = "Mass Extinction Coef: "+nam(n,0:16)
  ;plot = gsn_csm_contour(wks,u,res)       ; contour the variable
;  plot = gsn_csm_contour(wks,var,res)
 plot  = gsn_csm_xy (wks,ave, ave&y,res)


end
