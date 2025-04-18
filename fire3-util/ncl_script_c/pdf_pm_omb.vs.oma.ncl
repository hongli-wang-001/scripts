;***********************************************************
; pdf_2.ncl
;
; Concepts illustrated:
;   - Generating univariate distributions of different sizes
;   - Using the returned 'bin_center' attribute to assign the appropriate x-axis coordinates
;   - Using bar chart resources to give the plot a different look
;***********************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   

begin
;  a   = random_normal (  0,75, (/ 1000/))
;  b   = random_normal ( 25,20, (/10,40/))
;  c   = random_normal (  5,50, (/  500/))

  op2         = True
  op2@bin_min = -40.
  op2@bin_max = 40.
;  op2@bin_min = 1.0*minvali
;  op2@bin_max = 1.0*maxvali
;  nbn = maxvali-minvali
  nbn = 40


   diri = ""
   fili = systemfunc("ls ./pe*pm2_5_01.nc4")
   nfili= dimsizes(fili)
   print("nfili1="+nfili)
   f    = addfiles(diri+fili,"r")
;   lat    = f[:]->Latitude
;   lon    = f[:]->Longitude
   obs    = f[:]->Observation ;
   R      = f[:]->Obs_Minus_Forecast_unadjusted ;Observation
   aflag  = f[:]->Analysis_Use_Flag
   R = mask(R,(aflag.eq.1), True)
   a = R
   b = R
   delete(fili)
   delete(f)
   delete(obs)
   delete(R) 
   printMinMax (a, False)
   printMinMax (b, False)

  ap  = pdfx(a, nbn, op2)    ; specify 40 bins
  asciiwrite("ap.txt",ap)
  bp  = pdfx(b, nbn, op2)    ; rather than the default 25
  asciiwrite("bp.txt",bp)


   fili = systemfunc("ls  ./pe*pm2_5_02.nc4")
   nfili= dimsizes(fili)
   print("nfili2="+nfili)
   f    = addfiles(diri+fili,"r")
;   lat    = f[:]->Latitude
;   lon    = f[:]->Longitude
   obs    = f[:]->Observation ;
   R      = f[:]->Obs_Minus_Forecast_unadjusted ;Observation
   R = mask(R,(aflag.eq.1), True)
   c = R
   delete(fili)
   delete(f)
   delete(obs)
   delete(R)
   printMinMax (c, False)
  cp  = pdfx(c, nbn, op2)
  asciiwrite("cp.txt",cp)

   ;diri = "vrf_da_b/2019080412_i2019080400/"
   ;fili = systemfunc("cd "+diri+" ; ls pe*pm2_5_01.nc4")
   fili = systemfunc("ls ./pe*pm2_5_01.nc4")
   nfili= dimsizes(fili)
   print("nfili3="+nfili)
   f    = addfiles(diri+fili,"r")
;   lat    = f[:]->Latitude
;   lon    = f[:]->Longitude
   obs    = f[:]->Observation ;
   R      = f[:]->Obs_Minus_Forecast_unadjusted ;Observation
   R = mask(R,(aflag.eq.1), True)
   d = obs-R
   delete(fili)
   delete(f)
   delete(obs)
   delete(R)
  printMinMax (d, False)
  dp  = pdfx(d, nbn, op2)
  asciiwrite("dp.txt",dp)
;------------------------------------------------
; 
;------------------------------------------------
  nVar    = 4
  nBin    = ap@nbins          ; retrieve the number of bins

  xx      = new ( (/nVar, nBin/), typeof(ap))

  xx(0,:) = ap@bin_center     ; assign appropriate "x" axis values
  xx(1,:) = bp@bin_center
  xx(2,:) = cp@bin_center
  xx(3,:) = cp@bin_center

  yy      = new ( (/nVar, nBin/), typeof(ap))
  yy(0,:) = (/ ap /)
  yy(1,:) = (/ bp /)
  yy(2,:) = (/ cp /)
  yy(3,:) = (/ dp /)

  wopt = True
  wopt@fout  = "frq.txt"
  write_matrix (transpose(yy), "5f13.4" , wopt)

  wks  = gsn_open_wks ("x11","frq_pm2p5_00Z")                ; send graphics to PNG file
  res  = True
  res@xyLineThicknesses        = (/2.0,2.0,2.0,2.0/)        
  res@xyLineColors             = (/"black","red","green","blue"/)  
  res@xyMonoDashPattern        = True              ; all solid 
  res@tiYAxisString            = "PDF (%)"

  ;res@gsnXYBarChart            = True              ; Create bar plot
  ;res@gsnXYBarChartOutlineOnly = True

  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .25                 ; move units right
  res@pmLegendOrthogonalPosF = -0.4                ; move units down
  res@pmLegendWidthF         = 0.125               ; Change width and
  res@pmLegendHeightF        = 0.15                ; height of legend.
  res@lgPerimOn              = True                ; turn off/on box around
  res@lgLabelFontHeightF     = .015                ; label font height
  res@xyExplicitLegendLabels = (/"OBS","BKG","ANL"/)  ; create explicit labels
  res@tiMainString           = "PDF: "+nVar+" Distributions: " + nBin +" Bins"
;  res@tiMainString           = "PDF: Three Distributions: 40 Bins"
  plot = gsn_csm_xy (wks, xx(0:2,:), yy(0:2,:), res)
end
