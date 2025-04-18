;----------------------------------------------------------------------
; scatter_7.ncl
;
; Concepts illustrated:
;   - Drawing a scatter plot with markers of different colors and sizes
;   - Annotating a plot with labelbars, text, and markers
;   - Creating a color map using named colors
;   - Drawing polymarkers and text in NDC space
;----------------------------------------------------------------------
; The original version of this script was contributed by Larry
; McDaniel of IMAGe/NCAR. It had several lines of data processing
; calls before the plotting section.
;
; This script is a pared down version of the original one, with 
; the data processing calls removed. The values are read from
; two simple ASCII files. 
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
;----------------------------------------------------------------------
; Data is hard-coded here. It was originally calculated by running
; a compute-intensive script.
;----------------------------------------------------------------------
;  avdataT = (/1.965796,1.777351,1.915377,1.818644,1.766584,1.56958, \
;              2.101717,1.859297,1.662189/)
;  avdataP = (/-8.902548,-9.490757,-9.976177,-1.000492,2.248538,\
;              -2.698187,9.177564,3.478979,4.059988/)
   diri = "./"
   fili = systemfunc("cd "+diri+" ; ls pe*pm2_5_01.nc4")
   nfili= dimsizes(fili)
   print("nfili="+nfili)
   a    = addfiles(diri+fili,"r")
   lat    = a[:]->Latitude
   lon    = a[:]->Longitude
   obs      = a[:]->Observation ;Obs_Minus_Forecast_unadjusted ;Observation
   R = a[:]->Obs_Minus_Forecast_unadjusted
   aflag  = a[:]->Analysis_Use_Flag
   errinvf = a[:]->Errinv_Final
   errinvi = a[:]->Errinv_Input
   ;R = obs-R
   printMinMax (R, False)
   R = mask(R,(aflag.eq.1), True)
   errinvf = mask(errinvf,(aflag.eq.1), True)
   errinvi = mask(errinvi,(aflag.eq.1), True) 
   printMinMax (R, False)
   avdataT = R
   avdataP = 1.0/errinvi
 
;----------------------------------------------------------------------
; Open a PNG file and define the color table
;----------------------------------------------------------------------
;  colors  = (/"white","black","dodgerblue2","firebrick","forestgreen",\
;              "dodgerblue2","darkorange2","forestgreen","deeppink2",\
;              "darkorange2","purple"/)

  wks = gsn_open_wks ("x11","scatter")
;  gsn_define_colormap(wks,colors)
  gsn_define_colormap(wks,"BlueRed")

;----------------------------------------------------------------------
; Set some plot resources
;----------------------------------------------------------------------
  res                    = True

  res@gsnDraw            = False       ; Don't draw plot
  res@gsnFrame           = False       ; Don't advance frame

;
; Don't use gsnMaximize. Instead, we will control 
; exactly where plot goes.
;
  res@vpXF              = 0.2
  res@vpYF              = 0.85
  res@vpWidthF          = 0.7
  res@vpHeightF         = 0.6

;---Some calculations to get min/max for Y axis.
  miny = min(avdataP) 
  maxy = max(avdataP) 
  minx = min(avdataT) 
  maxx = max(avdataT) 

  amaxy = maxy - abs(miny)
  amaxy = amaxy / 2.
  amaxy = amaxy + 25.
  
  aminy = maxy - abs(miny)
  aminy = aminy / 2.
  aminy = aminy - 25.

;---Resources for X/Y axes.
  res@trYMinF            =  1 ;aminy               ; min value on y-axis
  res@trYMaxF            =  4 ;amaxy               ; max value on y-axis
  res@trXMinF            =  -20.                  ; min value on x-axis
  res@trXMaxF            =  20.                  ; max value on x-axis

;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.03
  res@xyMarkerThicknessF = 3.
  res@xyMonoMarkerColor  = False          ; Allow different colors for markers
;  res@xyMarkerColors     = colors(2:)     ; It's okay to list more than
                                           ; you need here
;  res@xyMarkers          = (/4,4,4,7,7,7,3,6,6/)
  res@xyMarkers          =  1 
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "NavyBlue"               ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)v
 
;---Tickmark resources
  res@tmXBMode         = "Manual"
  res@tmYLMode         = "Manual"
  res@tmXBTickSpacingF = 5.
  res@tmYLTickSpacingF = 1.

;---Title resources
  res@tiMainString       = "OBS Error .vs. OMB: PM2.5"
  res@tiXAxisString      = "OMB "
  res@tiYAxisString      = "OBS Error"
  res@tiMainFontHeightF  = 0.03

;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------
  plot = gsn_csm_xy (wks,transpose((/avdataT,avdataT/)), \
                         transpose((/avdataP,avdataP/)),res)

;----------------------------------------------------------------------
; Draw time stamp first, then the plot.
;----------------------------------------------------------------------
  infoTimeStamp(wks,0.005,"")
  draw(plot)


  frame(wks)    ; Now advance the frame.
end
