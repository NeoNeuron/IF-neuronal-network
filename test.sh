#!/bin/bash
# temporal execution script
# ./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv 1 1000,20000 raster.csv
# cp ./data/raster/raster.csv ./data/raster/raster_x.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul14/t3/I.csv 1 1000,20000 lfp_x.csv
# ./bin/potential.out /media/kyle/Drive/ResearchData/Jul14/t3/V.csv 1 1000,20000 potential_x.csv

for ((i = 0; i < 4; i ++))
do
  # ./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv $i 1000,20000 raster_y.csv
  ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul21/t1/postI.csv $i 1000,20000 lfp.csv
  # ./bin/potential.out /media/kyle/Drive/ResearchData/Jul14/t3/V.csv $i 1000,20000 potential.csv
  cp ./data/lfp/lfp.csv ./data/lfp/lfp_y.csv
  # cp ./data/potential/potential.csv ./data/potential/potential_y.csv
  ./bin/mi.out 2 0.125 80,200
  ./bin/sta.out data/raster/raster.csv data/lfp/lfp.csv -10,25
  # ./bin/sta.out data/raster/raster.csv data/potential/potential.csv -10,25
  python pys/py.py ll$i.png
  # python py.py rr3.png
done
