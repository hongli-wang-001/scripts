dset fsi.dat 
title "FSI"
undef -99999.9
options yrev
ydef 181 linear -90.000000 1
xdef 360 linear 0.000000 1.000000
zdef 6   levels 850 700 500 250 200 100 
tdef 1 linear 12Z15mar2015  1dy
vars 3
u      6    1,  1,  0,  0 Surface pressure [hPa]
v      6    1,  1,  0,  0 Surface pressure [hPa]
t      6    1,  1,  0,  0 Surface pressure [hPa]
endvars

