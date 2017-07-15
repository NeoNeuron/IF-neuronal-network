#!/bin/bash
# temporal execution script
./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv 0 0,20000 raster.csv
./bin/lfp.out /media/kyle/Drive/ResearchData/Jul14/t3/I.csv 2 0,20000 lfp.csv
./bin/mi.out 1 0.25 100,100
