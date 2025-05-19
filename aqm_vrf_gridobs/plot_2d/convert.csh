convert -delay 1 `ls -l obs*.png | grep -E '(00.png|03.png|06.png|09.png|12.png|15.png|18.png|21.png)' | awk '{print $9}'` OBS_PM2p5_every3h_delay1.gif
convert -delay 1 `ls -l obs*.png | grep -E '(00.png|06.png|12.png|18.png)' | awk '{print $9}'` OBS_PM2p5_every6h_delay1.gif
convert -delay 5 `ls -l obs*.png | grep -E '(00.png|06.png|12.png|18.png)' | awk '{print $9}'` OBS_PM2p5_every6h.gif
