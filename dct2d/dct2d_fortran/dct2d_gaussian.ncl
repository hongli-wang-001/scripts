;ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
; Author : Hongli Wang   MMM/NCAR   01/27/2013 
;=========================================================;
;
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;
external ex_dct "/glade/u/home/hlwang/util/dct2df90.so"
;
  dx=1.0
  kd=1
  av = 0.0
  sd = 1.0
  jymd=256
  ixmd=512
  nspt=255
  spt1d =  new ((/nspt/),float)
  length1d = new ((/nspt/),float)
  nor = random_normal(av,sd,(/jymd,ixmd/))
  ex_dct :: dct2d(ixmd,jymd,nspt,kd,dx,nor,spt1d,length1d)

  ixx = ispan(1,nspt,1)

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("eps","gaussian_spctra")                ; open workstation

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Gaussian Spectra"          ; add title
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True

 plot  = gsn_csm_xy (wks,ixx,spt1d,res) ; create plot

;
;  fmtf = nspt + "f9.2"
;  wopt = True
;  wopt@fout = "rmd.txt"
;  write_matrix(rmd(:,1:nend,fmtf,wopt)

