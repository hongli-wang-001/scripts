#!/bin/csh -x
#
#Total 70 aero species
#I: 14 J: K:7
touch lam_cmaq.yaml.cmaq.T70
#foreach avar (aso4i  ano3i   anh4i   anai    acli    aeci    aothri alvpo1i asvpo1i asvpo2i alvoo1i alvoo2i asvoo1i asvoo2i ) 
#foreach avar (aso4j  ano3j   anh4j   anaj    aclj    aecj    aothrj afej    asij    atij    acaj    amgj    amnj    aalj   ) 
#foreach avar (akj    alvpo1j asvpo1j asvpo2j asvpo3j aivpo1j axyl1j axyl2j  axyl3j  atol1j  atol2j  atol3j  abnz1j  abnz2j ) 
#foreach avar (abnz3j aiso1j  aiso2j  aiso3j  atrp1j  atrp2j  asqtj  aalk1j  aalk2j  apah1j  apah2j  apah3j  aorgcj  aolgbj ) 
foreach avar (aolgaj alvoo1j alvoo2j asvoo1j asvoo2j asvoo3j apcso  aso4k   asoil   acors   aseacat aclk    ano3k   anh4k )

cat <<EOF >> lam_cmaq.yaml.cmaq.T70  
- long name: ${avar}
  io name: ${avar}
  units: ugkg-1

EOF

end
