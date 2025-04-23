; Verify WRF ARW forecasts against GFS met analysis
; Hongli Wang 02-01-2013
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external ex_dct "/glade/u/home/hlwang/util/dct2df90.so"

begin

  start_date = "2012082000"
  end_date   = "2012083000"

  VERIFY_HOUR     = (/00/)   ;  (/00,24,48,72/)   ; (/0, 12, 24, 36/)
  vNAM            = (/"U","V","T","Z","RH"/) ; (/"U"/)  ; (/"U","V","T","Z","RH"/) 
  UnitLabel       = (/"m s-1","m s-1","C","m","%"/)
  pressure_levels = (/500./) 

  interval   = 24 

  REG_DIR = "/glade/scratch/hlwang/afwa/tc_36/"
  EXP_DIR = "con_b_cq1/"
  RC_DIR  = REG_DIR + "rc/"
  FC_DIR  = REG_DIR + EXP_DIR + "fc/"

  SUFFIX  = "conb-gfs"

; read a WRF output to get dimensions
    ymdh = start_date
    yr=str_get_cols(ymdh,0,3)
    mn=str_get_cols(ymdh,4,5)
    dy=str_get_cols(ymdh,6,7)
    hr=str_get_cols(ymdh,8,9)
    WRF_DATE=yr + "-" + mn + "-" + dy + "_" + hr + ":00:00"
;    print("Verifying " + WRF_DATE )

  FC_FIL=FC_DIR + ymdh + "/" + "wrfout_d01_" + WRF_DATE
  f=addfile(FC_FIL + ".nc","r")

  attNam="WEST-EAST_GRID_DIMENSION"
  nLon=f@$attNam$ - 1
  attNam="SOUTH-NORTH_GRID_DIMENSION"
  nLat=f@$attNam$ - 1
  nTim=dimsizes(VERIFY_HOUR)
  nVar=dimsizes(vNAM)
  nLev=dimsizes(pressure_levels)
  attNam="DX"
  dx=f@$attNam$/1000.0  ; to KM 

  nspt  = nLat-1
  spt1d = new ((/nspt/),float)
  length1d = new ((/nspt/),float)
  spttot= new ((/nVar,nspt/),float)
  spttot= 0.0

  RMSE=new((/nTim, nVar, nLev, nLat, nLon/), "float")
  RMSE = 0 
  DIFF=new((/nTim, nVar, nLev, nLat, nLon/), "float")
  DIFF= 0

  ;printVarSummary(RMSE)

  nTim = nTim - 1
  nVar = nVar - 1
  nLev = nLev - 1

  quotes = integertochar(34)
  ymdh = start_date
  n = 0
  do while ( ymdh .le. end_date )

    yr=str_get_cols(ymdh,0,3)
    mn=str_get_cols(ymdh,4,5)
    dy=str_get_cols(ymdh,6,7)
    hr=str_get_cols(ymdh,8,9)

    INIT_DATE=yr + "-" + mn + "-" + dy + " " + hr + ":00"
    print("Verifying " + INIT_DATE )

    do t=0, nTim

      ;print("  For " + sprinti("%0.2i" ,VERIFY_HOUR(t)) + "H" )

;whl   if ( VERIFY_HOUR(t) .eq. 0 ) then
;        VERIFY_DATE = ymdh
;        FC_FIL= REG_DIR + EXP_DIR + ymdh + "/wrf/wrfinput_d01" 
;      else
        CMD="date -u -d " + \
             quotes + VERIFY_HOUR(t) + " hours " + INIT_DATE + quotes + " +" + \
             quotes + "%Y%m%d%H" + quotes
        VERIFY_DATE = systemfunc(CMD)

        yr=str_get_cols(VERIFY_DATE,0,3)
        mn=str_get_cols(VERIFY_DATE,4,5)
        dy=str_get_cols(VERIFY_DATE,6,7)
        hr=str_get_cols(VERIFY_DATE,8,9)

        WRF_DATE=yr + "-" + mn + "-" + dy + "_" + hr + ":00:00"

        FC_FIL  = FC_DIR + ymdh + "/" + "wrfout_d01_" + WRF_DATE
;      end if
      RC_FIL=RC_DIR + VERIFY_DATE + "/" + "wrfinput_d01"
      print("  RC : " + RC_FIL )
      print("  FC : " + FC_FIL )

      RC=addfile(RC_FIL + ".nc", "r")
      P_RC = wrf_user_getvar(RC, "pressure", -1)

      FC=addfile(FC_FIL + ".nc", "r")
      P_FC = wrf_user_getvar(FC, "pressure", -1)

      do i=0, nVar

        ;print("    processing for " + vNAM(i) )

        if ( vNAM(i) .eq. "U" ) then
          var="ua"
        else if ( vNAM(i) .eq. "V" ) then
          var="va"
        else if ( vNAM(i) .eq. "T" ) then
          var="tc"
        else if ( vNAM(i) .eq. "Z" ) then
          var="z"
        else if ( vNAM(i) .eq. "RH" ) then
          var="rh"
        else
        end if
        end if
        end if
        end if
        end if

        VAR_FC = wrf_user_getvar(FC,var,-1)
        FC_PLANE = wrf_user_intrp3d(VAR_FC,P_FC,"h",pressure_levels,0.,False)

        VAR_RC = wrf_user_getvar(RC,var,-1)
        RC_PLANE = wrf_user_intrp3d(VAR_RC,P_RC,"h",pressure_levels,0.,False)

        FC_PLANE=FC_PLANE-RC_PLANE

;        printVarSummary(FC_PLANE)

        ;if(any(ismissing(VAR_RC))) then
        ;  print("      " + vNAM(i) + " contains some missing values. Beware.")
        ;end if

        ;if(any(ismissing(VAR_FC))) then
        ;  print("      " + vNAM(i) + " contains some missing values. Beware.")
        ;end if

        ;if(any(ismissing(RC_PLANE))) then
        ;  print("      " + vNAM(i) + "_PLANE contains some missing values. Beware.")
        ;end if

        ;if(any(ismissing(FC_PLANE))) then
        ;  print("      " + vNAM(i) + "_PLANE contains some missing values. Beware.")
        ;end if
;  dct
  kd    = 1
;  nspt  = nLat-1
;  spt1d = new ((/nspt/),float)
;  length1d = new ((/nspt/),float)
  ex_dct :: dct2d(nLon,nLat,nspt,kd,dx,FC_PLANE(0,:,:),spt1d,length1d)
  spttot(i,:) = spttot(i,:) + spt1d

  ixx = ispan(1,nspt,1)

;************************************************
; plotting parameters
;************************************************
 ftype = "eps"
 pre_figname = SUFFIX+"_spt1d_"+vNAM(i)+"_"+nTim+"_"+ymdh
 figname = pre_figname+"."+ftype
; wks   = gsn_open_wks (ftype,figname)                            ; open workstation

 res                  = True                                 ; plot mods desired
 res@tiMainString     = "Variance Spetra: "+vNAM(i)          ; add title
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True

; plot  = gsn_csm_xy (wks,ixx,spt1d,res) ; create plot

;
;  spt1 = new ((/1,nspt/),float)
;  spt1 = spt1d
;  fmtf = nspt + "E11.3"
;  wopt = True
;  wopt@fout = pre_figname+"."+"txt" 
;  write_matrix(spt1,fmtf,wopt)
;  system ("mv spt1d.txt " + pre_figname+"."+"txt")

      end do  ; End of do loop for nVar
    end do    ; End of do loop for nTim

    n = n + 1

    CMD="date -u -d " + \
         quotes + interval + " hours " + INIT_DATE + quotes + " +" + \
         quotes + "%Y%m%d%H" + quotes
    ymdh = systemfunc(CMD)

  end do  ; End of do loop for ynmd 

  spttot = spttot/n
; print(spttot)
; print(length1d)
  spt = new ((/1,nspt/),float)
  spt(0,:) = length1d

  fmtf = nspt + "E11.3"
  wopt = True
  wopt@fout = "length"+"."+"txt"
  write_matrix(spt,fmtf,wopt)
;************************************************
; plotting parameters
;************************************************
 ftype = "eps"
 do i=0,nVar-1
 pre_figname = "spt1d_"+vNAM(i)+"_time_ave"
 figname = SUFFIX+pre_figname+"."+ftype
 wks   = gsn_open_wks (ftype,figname)                            ; open workstation

 res                  = True                                 ; plot mods desired
 res@tiMainString     = "Variance Spetra: "+vNAM(i)          ; add title
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True

 plot  = gsn_csm_xy (wks,ixx,spttot(i,:),res) ; create plot

 end do
; nVar variables together
 figname = SUFFIX+"_spt1d_vars_time_ave."+ftype
 wks   = gsn_open_wks (ftype,figname) 
 res                  = True                                 ; plot mods desired
 res@tiMainString     = "Variance Spetra: "+VERIFY_HOUR(0)+"Z"          ; add title
 res@tiXAxisString = "Length Scale (km)"  ; xaxis string
 res@tiYAxisString = "Variance Spectra"                  ; yaxis string 
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True
 res@xyLineThicknesses = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0/)   ; make 2nd lines thicker
 res@xyDashPatterns    = (/0,1,2,3,4,5,6,7/)
 res@trYMaxF           =   100.                    ; axis max
 res@trYMinF           =   0.000001                    ; axis min

 res@xyLabelMode            = "Custom"             ; label a line
 res@xyExplicitLabels       = vNAM                 ; text to use
 res@xyLineLabelFontHeightF = 0.020                ; font height
; res@xyLineLabelConstantSpacingF = 10.0
 res@xyLineLabelFontColors       = (/"black","gray","red","blue","purple","brown","orange","coral"/)     ; label color
 res@xyLineColors                = (/"black","gray","red","blue","purple","brown","orange","coral"/)

; xaxislbl = (/length1d(0),length1d(9),length1d(99)/)
 xaxislbl = (/"8352","5011","2277","1193","533","248","125"/)
 xaxisloc = (/1,2,5,10,23,50,100/)
; xaxislbl = length1d
; xaxisloc = ispan(1,nspt,1)
 res@tmXBOn      = True
 res@tmXBMode    = "Explicit"
 res@tmXBValues  = xaxisloc 
 res@tmXBLabels  = xaxislbl
; res@tmLabelAutoStride = True    ; nice stride on labels  

 plot  = gsn_csm_xy (wks,ixx,spttot,res) ; create plot

  fmtf = nspt + "E11.3"
  wopt = True
  wopt@fout = SUFFIX+"_spt1d_vars_time_ave."+"txt"
  write_matrix(spttot,fmtf,wopt)

end
