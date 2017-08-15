#!/bin/bash
# temporal execution script
# ./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv 1 1000,20000 raster.csv
# cp ./data/raster/raster.csv ./data/raster/raster_x.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul14/t3/I.csv 1 1000,20000 lfp_x.csv
# ./bin/potential.out /media/kyle/Drive/ResearchData/Jul14/t3/V.csv 1 1000,20000 potential_x.csv

# for ((i = 0; i < 4; i ++))
# do
#   # ./bin/raster.out /media/kyle/Drive/ResearchData/Jul14/t3/raster.csv $i 1000,20000 raster_y.csv
#   ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul21/t1/postI.csv $i 1000,20000 lfp.csv
#   # ./bin/potential.out /media/kyle/Drive/ResearchData/Jul14/t3/V.csv $i 1000,20000 potential.csv
#   cp ./data/lfp/lfp.csv ./data/lfp/lfp_y.csv
#   # cp ./data/potential/potential.csv ./data/potential/potential_y.csv
#   ./bin/mi.out 2 0.125 80,200
#   ./bin/sta.out data/raster/raster.csv data/lfp/lfp.csv -10,25
#   # ./bin/sta.out data/raster/raster.csv data/potential/potential.csv -10,25
#   python pys/py.py ll$i.png
#   # python py.py rr3.png
# done

# 5,32,37,43,45

# ./bin/raster.out /media/kyle/Drive/ResearchData/Jul23/t1/rasterPre.txt 0 1000,20000 raster1.csv
# ./bin/raster.out /media/kyle/Drive/ResearchData/Jul23/t2/rasterPre.txt 0 1000,20000 raster2.csv
# ./bin/raster.out /media/kyle/Drive/ResearchData/Jul23/t3/rasterPre.txt 0 1000,20000 raster3.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul23/t1/postI.txt 5,32,37 1000,20000 lfp1.csv
# # ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul23/t2/postI.txt 5,32,37 1000,20000 lfp2.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul23/t3/postI.txt 5,32,37 1000,20000 lfp3.csv
# cp ./data/raster/raster3.csv ./data/raster/raster.csv
# cp ./data/lfp/lfp3.csv ./data/lfp/lfp.csv
# ./bin/mi.out 1 0.25 40,100
# ./bin/sta.out data/raster/raster.csv data/lfp/lfp.csv -10,25
# python pys/py.py bb

./bin/raster.out ~/github/ifnet/data/tmp/raster.csv 0 1000,20000 raster.csv
# ./bin/lfp.out /media/kyle/Drive/ResearchData/Jul27/t0005/preI.txt 0 1000,20000 lfp_p.csv
for ((i = 0; i < 1; i ++))
do
  # ./bin/raster.out /media/kyle/Drive/ResearchData/Jul27/t0005/rasterPost.txt $i 1000,20000 raster_p.csv
  ./bin/lfp.out ~/github/ifnet/data/tmp/I.csv $i 1000,20000 lfp.csv
  # ./bin/potential.out /media/kyle/Drive/ResearchData/Jul27/t0005/postV.txt $i,1,2,3,4,5,6,7,8,9  1000,20000 potential.csv
  # cp ./data/raster/raster3.csv ./data/raster/raster.csv
  # cp ./data/lfp/lfp_y.csv ./data/lfp/lfp.csv
  ./bin/mi.out 1 0.25 100,200
  ./bin/sta.out data/raster/raster.csv data/lfp/lfp.csv -25,50
  # ./bin/lcc.out data/raster/raster.csv data/lfp/lfp.csv 0.25 100,200
  python pys/py.py fig65_single.png
done
