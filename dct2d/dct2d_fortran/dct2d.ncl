;ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
; Author : Hongli Wang   MMM/NCAR   01/27/2013 
; Input: 1. WRF file to provide lat and lon: fg
;           two WRF output files to obtain accumulated rainfall: fcst1.nc fcst2.nc; fg can be fcst1 or fcst2
;        2. St4 rainfall data
;        3. three external f90 subroutines
;        4. model precipitation
; Output: 1. rmd.txt  : st4 at model grid 
;         2. fss.txt and fbs.txt  
;         3. ts.txt
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
  dx=36.0
  fc1= addfile("dd-box.nc","r")
  fc2= addfile("fcst2.nc","r")
  theta1 = fc1->T
  theta2 = fc2->T
  mu  = fc1->MU
;  dtheta = theta2-theta1
  idimm=dimsizes(mu)
  print(idimm)
  ixmd=idimm(2)
  jymd=idimm(1)
  kzmd=idimm(0)
  nspt=jymd-1 ; min(ix,jy)-1

;  spt1d =  new ((/nspt/),float)
;  length1d = new ((/nspt/),float)
  kd=1 
  av = 0.0
  sd = 1.0
;  nor = random_normal(av,sd,(/jymd,ixmd/))
;  ex_dct :: dct2d(ixmd,jymd,nspt,kd,dx,mu(0,:,:),spt1d,length1d)
;  ex_dct :: dct2d(ixmd,jymd,nspt,kd,dx,nor,spt1d,length1d)
  jymd=512
  ixmd=512
  nspt=511
  spt1d =  new ((/nspt/),float)
  length1d = new ((/nspt/),float)
  nor = random_normal(av,sd,(/jymd,ixmd/))
  ex_dct :: dct2d(ixmd,jymd,nspt,kd,dx,nor,spt1d,length1d)

  ixx = ispan(1,nspt,1)

;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("ps","xy")                ; open workstation

 res                  = True                     ; plot mods desired
 res@tiMainString     = "Basic XY plot"          ; add title
 res@xyYStyle         = "Log"
 res@xyComputeYMin    = True
 res@xyXStyle         = "Log"
 res@xyComputeXMin    = True

 plot  = gsn_csm_xy (wks,ixx,spt1d,res) ; create plot

;
  fmtf = nspt + "f9.2"
  wopt = True
;  wopt@fout = "rmd.txt"
;  write_matrix(rmd(:,1:nend,fmtf,wopt)

