#!/bin/bash
# temporal execution script
./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv 1 1000,20000 raster.csv
# cp ./data/raster/raster.csv ./data/raster/raster_x.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul14/t3/I.csv 1 0,20000 lfp_x.csv

for ((i = 0; i < 4; i ++))
do
  # ./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv $i 1000,20000 raster_y.csv
  ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul14/t3/I.csv $i 1000,20000 lfp.csv
  # cp ./data/lfp/lfp_y.csv ./data/lfp/lfp.csv
  ./bin/mi.out 1 0.25 40,100
  ./bin/sta.out data/raster/raster.csv data/lfp/lfp.csv -10,25
  python py.py
done
