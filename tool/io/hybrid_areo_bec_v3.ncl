
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  A = addfile("./std1.nc","rw")
  B = addfile("./std2.nc","r")
  C = addfile("./std3.nc","r")
vNames = (/ "aso4i","aso4j","aso4k","acli","ano3i","ano3j","ano3k","anh4i","anh4j","anh4k", "aclj", "aclk","anai", "anaj","aseacat","aeci","aecj", "amgj", "akj",  "acaj","afej", "aalj","asij","atij","amnj","asoil","acors","aothri","aothrj","axyl1j","axyl2j","axyl3j","atol1j","alvpo1i","alvpo1j","asvpo1i","asvpo1j","asvpo2i","asvpo2j","asvpo3j","aivpo1j","alvoo1i","alvoo2i","asvoo1i","asvoo2i","atol2j","atol3j","abnz1j","abnz2j","abnz3j","abnz3j","aiso1j","aiso2j","aiso3j","atrp1j","atrp2j","asqtj","aalk1j","aalk2j","apah1j","apah2j","apah3j","aorgcj","aolgbj","aolgaj","alvoo1j","alvoo2j","asvoo1j","asvoo2j","asvoo3j","apcsoj" /)

 nNames = dimsizes (vNames)   ; Number of variables on the file

 do n=0, nNames-1              ; Loop through each variable

    print(n+" "+vNames(n))
    va = A->$vNames(n)$        ; Read the variable to memory
    vb = B->$vNames(n)$ 
    vc = C->$vNames(n)$
    gName = vNames(n) 
    va = 0.0*va + 0.8*vb+1.0*vc
    A->$gName$ = va
    delete(va)
    delete(vb)
 end do

end
