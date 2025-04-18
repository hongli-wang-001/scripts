
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  A = addfile("./file1_fcst.nc","rw")
  D = addfile("./pm25.nc","r")
vNames = (/ "aso4i","aso4j","aso4k","acli","ano3i","ano3j","ano3k","anh4i","anh4j","anh4k", "aclj", "aclk","anai", "anaj","aseacat","aeci","aecj", "amgj", "akj",  "acaj","afej", "aalj","asij","atij","amnj","asoil","acors","aothri","aothrj","axyl1j","axyl2j","axyl3j","atol1j","alvpo1i","alvpo1j","asvpo1i","asvpo1j","asvpo2i","asvpo2j","asvpo3j","aivpo1j","alvoo1i","alvoo2i","asvoo1i","asvoo2i","atol2j","atol3j","abnz1j","abnz2j","abnz3j","abnz3j","aiso1j","aiso2j","aiso3j","atrp1j","atrp2j","asqtj","aalk1j","aalk2j","apah1j","apah2j","apah3j","aorgcj","aolgbj","aolgaj","alvoo1j","alvoo2j","asvoo1j","asvoo2j","asvoo3j","apcsoj" /)

 scale = A->aso4i      
 scale = 1.0

 pm2d  = A->pm25_tot(0,63,:,:)
 
 pmobs = D->interpolated_pm25 
 ndxy = dimsizes(pmobs)
 ndx = ndxy(1)
 ndy = ndxy(0)
 print(ndx+" "+ndy)
 print(dimsizes(pmobs))
 print(dimsizes(pm2d))

 ; Set NaN to missing values
  if (any(isnan_ieee(pmobs))) then
    if(.not.isatt(pmobs,"_FillValue")) then
      pmobs@_FillValue = 1.e20
    end if
    replace_ieeenan (pmobs, pmobs@_FillValue, 0)
  end if

  opt = True
  opt@PrintStat = True
  statb = stat_dispersion(pmobs, opt )
 diff = pmobs-pm2d

 do j=0,ndy-1
 do i=0,ndx-1
    if (.not.ismissing(diff(j,i)) .and. diff(j,i).gt.0) then
       if (pm2d(j,i).gt.0) then
         if (pmobs(j,i)/pm2d(j,i).ge.10.0) then
           print(j+" "+i+" "+diff(j,i)+" "+pmobs(j,i)+" "+pm2d(j,i)+" "+pmobs(j,i)/pm2d(j,i))
         end if
         scale(0,5:63,j,i)=min((/pmobs(j,i)/pm2d(j,i)*2.0,100.0/))
       else
         print("BKG PM is 0! "+j+" "+i+" "+diff(j,i)+" "+pmobs(j,i)+" "+pm2d(j,i))
       end if
    end if
 end do
 end do

 nNames = dimsizes (vNames)   ; Number of variables on the file
 do n=0, nNames-1              ; Loop through each variable

    print(n+" "+vNames(n))
    va = A->$vNames(n)$        ; Read the variable to memory
    gName = vNames(n) 
    va = va*scale
    va = where(va.ge.10000.0, 10000.0, va)
    A->$gName$ = va
    delete(va)
 end do

end
