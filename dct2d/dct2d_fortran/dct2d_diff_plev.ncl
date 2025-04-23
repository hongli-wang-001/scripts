; Verify WRF ARW forecasts against GFS met analysis
; Hongli Wang and Jianyu Liu   01-30-2013
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external ex_dct "/glade/u/home/hlwang/util/dct2df90.so"

begin

  start_date = "2012082000"
  end_date   = "2012083000"

  VERIFY_HOUR     = (/48/)   ; (/0, 12, 24, 36/)
  vNAM            = (/"U","V","T","Z","RH"/) ; (/"U"/)  ; (/"U","V","T","Z","RH"/) 
  UnitLabel       = (/"m s-1","m s-1","C","m","%"/)
  pressure_levels = (/500./) 

  interval   = 24 

  REG_DIR = "/glade/scratch/hlwang/afwa/tc_36/"
  EXP_DIR = "con_b_cq1/"
  FC_DIR  = REG_DIR + EXP_DIR + "fc/"
  EXP_DIR2= "iasi_b_cyc_bias_cq1/"
  FC_DIR2 = REG_DIR + EXP_DIR2 + "fc/"


  SUFFIX  = "iasi-con"

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
        FC_FIL2 = FC_DIR2 + ymdh + "/" + "wrfout_d01_" + WRF_DATE
;      end if
     
      print("  FC : " + FC_FIL )
      print("  FC2: " + FC_FIL2 )
      FC=addfile(FC_FIL + ".nc", "r")
      FC2=addfile(FC_FIL2 + ".nc", "r")

      P_FC = wrf_user_getvar(FC, "pressure", -1)
      P_FC2= wrf_user_getvar(FC2, "pressure", -1)

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

        VAR_FC2 = wrf_user_getvar(FC2,var,-1)
        FC_PLANE2 = wrf_user_intrp3d(VAR_FC2,P_FC2,"h",pressure_levels,0.,False)

        FC_PLANE=FC_PLANE2-FC_PLANE

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
 wks   = gsn_open_wks (ftype,figname)                            ; open workstation

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
 res@tiMainString     = "Variance Spetra"          ; add title
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True
 res@xyLineThicknesses = (/1.0,1.0,1.0,1.0,1.0/)               ; make 2nd lines thicker

 res@xyLineColors      = (/"black","black","red","blue","gray"/)          ; change line color
 plot  = gsn_csm_xy (wks,ixx,spttot,res) ; create plot

  fmtf = nspt + "E11.3"
  wopt = True
  wopt@fout = SUFFIX+"_spt1d_vars_time_ave."+"txt"
  write_matrix(spttot,fmtf,wopt)

end
