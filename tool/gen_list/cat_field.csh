#!/bin/csh -x
#
#Total 70 aero species
#I: 14 J: K:7
#aso4i,ano3i,anh4i,anai,acli,aeci,aothri,alvpo1i,asvpo1i,asvpo2i,alvoo1i,alvoo2i,asvoo1i,asvoo2i,
#                        aso4j,ano3j,anh4j,anaj,aclj,aecj,aothrj,afej,asij,atij,acaj,amgj,amnj,aalj,akj,
#                        alvpo1j,asvpo1j,asvpo2j,asvpo3j,aivpo1j,axyl1j,axyl2j,axyl3j,atol1j,atol2j,atol3j,
#                        abnz1j,abnz2j,abnz3j,aiso1j,aiso2j,aiso3j,atrp1j,atrp2j,asqtj,aalk1j,aalk2j,apah1j,
#                        apah2j,apah3j,aorgcj,aolgbj,aolgaj,alvoo1j,alvoo2j,asvoo1j,asvoo2j,asvoo3j,apcsoj,
#                        aso4k,asoil,acors,aseacat,aclk,ano3k,anh4k
#touch addpm.add
#foreach avar (aso4i  ano3i   anh4i   anai    acli    aeci    aothri alvpo1i asvpo1i asvpo2i alvoo1i alvoo2i asvoo1i asvoo2i ) 
#foreach avar (aso4j  ano3j   anh4j   anaj    aclj    aecj    aothrj afej    asij    atij    acaj    amgj    amnj    aalj   ) 
#foreach avar (akj    alvpo1j asvpo1j asvpo2j asvpo3j aivpo1j axyl1j axyl2j  axyl3j  atol1j  atol2j  atol3j  abnz1j  abnz2j abnz3j aiso1j  aiso2j  aiso3j  atrp1j  atrp2j  asqtj  aalk1j  aalk2j  apah1j  apah2j  apah3j  aorgcj  aolgbj aolgaj alvoo1j alvoo2j asvoo1j asvoo2j asvoo3j apcsoj )
foreach avar (aso4k   asoil   acors   aseacat aclk    ano3k   anh4k )
#touch setfield.add
cat <<EOF >> setfield.add  
   atlas::Field traj_${avar} = afieldsetTraj.field("${avar}");
   atlas::Field tl_${avar} = afieldsetTL.field("${avar}");
EOF

cat <<EOF >> autofield.add  
    auto tl_${avar}_view = atlas::array::make_view<double, 2>(tl_${avar});
    auto traj_${avar}_view = atlas::array::make_view<double, 2>(traj_${avar});
EOF

cat <<EOF >> TLfield.add  
        tl_fine_particulate_matter_view(jnode, level) +=
                      tl_${avar}_view(jnode, level)*traj_pm25co_view(jnode, level) ;
EOF

cat <<EOF >> setADfield.add 
    atlas::Field traj_${avar} = afieldsetTraj.field("${avar}");
    atlas::Field ad_${avar} = afieldsetAD.field("${avar}");
EOF

cat <<EOF >> autoADfield.add  
    auto ad_${avar}_view = atlas::array::make_view<double, 2>(ad_${avar});
    auto traj_${avar}_view = atlas::array::make_view<double, 2>(traj_${avar});
EOF

cat <<EOF >> ADfield.add 
        ad_${avar}_view(jnode, level) +=
                      ad_fine_particulate_matter_view(jnode, level)*traj_pm25co_view(jnode, level);
EOF
end
