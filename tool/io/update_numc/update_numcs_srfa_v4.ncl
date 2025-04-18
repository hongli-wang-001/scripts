begin
;---Open files

   a        = addfile("bg_tracer.nc","r")
   c        = addfile("rst_tracer.nc","rw")

  numcNames = (/"numatkn","numacc","numcor"/)
  srfaNames = (/"srfatkn","srfacc","srfcor"/)

;   dstr  = asciiread("var.list", -1, "string")
 vNames = (/ "aso4i", "ano3i", "anh4i", "anai", "acli", "aeci", "aothri", "alvpo1i", "asvpo1i", "asvpo2i", "alvoo1i", "alvoo2i", "asvoo1i", "asvoo2i", "aso4j", "ano3j", "anh4j", "anaj", "aclj", "aecj", "aothrj", "afej", "asij", "atij", "acaj", "amgj", "amnj", "aalj", "akj", "alvpo1j", "asvpo1j", "asvpo2j", "asvpo3j", "aivpo1j", "axyl1j", "axyl2j", "axyl3j", "atol1j", "atol2j", "atol3j", "abnz1j", "abnz2j", "abnz3j", "aiso1j", "aiso2j", "aiso3j", "atrp1j", "atrp2j", "asqtj", "aalk1j", "aalk2j", "apah1j", "apah2j", "apah3j", "aorgcj", "aolgbj", "aolgaj", "alvoo1j", "alvoo2j", "asvoo1j", "asvoo2j", "asvoo3j", "apcsoj", "aso4k", "asoil", "acors", "aseacat", "aclk", "ano3k", "anh4k"/)
 nNames = dimsizes (vNames)   ; Number of variables on the file

; igroup = new(69,integer)
 igroup = (/1,1,1,1,  1,1,1, 1,1,1, 1,1,1,1,3,  1,2,  2,2,2,2,2,2,2,2,3,3,1,2,2,2,2,2,1,2,1,2,1,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/)
 igroup(0:13) = 1
 igroup(14:62) = 2
 igroup(63:) = 3
 gNames = (/"AT","AC", "CO"/)

  density = asciiread("density_v70.txt" , -1, "float")
  dstr = vNames
  printVarSummary(dstr)
  nvar = dimsizes(dstr)
  nden = dimsizes(density)
  print("nvar= "+nvar +" "+nden)
  do n=0,nvar-1
  print("VAR= " +dstr(n)+" "+igroup(n)+" "+gNames(igroup(n)-1) +" "+density(n))
  end do 

  nd_thin = 1
  it  = 0
  aso4j = a->aso4j(it,:,:,:)
  ndims = dimsizes(aso4j)
  ;print(ndims)

  ave = new((/3, ndims(0),ndims(1),ndims(2)/),float)
  ave_bg = new((/3, ndims(0),ndims(1),ndims(2)/),float)
  ave = 0.0
  ave_bg = 0.0
  scale = ave

  nl  = 63 
  do n=0,nvar-1
  print(dstr(n)+" "+igroup(n)+" "+gNames(igroup(n)-1))
  temp1  = a->$dstr(n)$(it,:,:,:)    
  temp2  = c->$dstr(n)$(it,:,:,:)  
  ;printVarSummary(temp1)
  ;printVarSummary(temp2)
  ;printVarSummary(ave)
  ;temp   = temp2-flt2dble(temp1)   
  ave_bg(igroup(n)-1,:,:,:)= ave_bg(igroup(n)-1,:,:,:)+temp1/density(n)
  ave(igroup(n)-1,:,:,:)= ave(igroup(n)-1,:,:,:)+temp2/density(n)
  end do 
  scale = 1.0
  scale = where(ave_bg>0.0,ave/ave_bg,1)

  do k=0,ndims(0)-1
  print(avg(scale(0,k,:,:))+" "+avg(scale(1,k,:,:))+" "+avg(scale(2,k,:,:)))
  end do

  do i=0,2
    print("Update_numc= "+numcNames(i))
    temp  = a->$numcNames(i)$
    printVarSummary(temp)
    temp(0,:,:,:)=scale(i,:,:,:)*temp(0,:,:,:)
    printVarSummary(temp)
    c->$numcNames(i)$ = temp 
    delete(temp)
  end do

  do i=0,2
    print("Update_srfa= "+srfaNames(i))
    temp  = a->$srfaNames(i)$
    printVarSummary(temp)
    temp(0,:,:,:)=scale(i,:,:,:)*temp(0,:,:,:)
    printVarSummary(temp)
    c->$srfaNames(i)$ = temp
    delete(temp)
  end do

end

