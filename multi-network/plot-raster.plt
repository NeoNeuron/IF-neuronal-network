set xrange [0:10000]
set yrange [-1:100]
set xlabel "Time(ms)"
set ylabel "Neuronal indices"

plot for [i=2:200] '/media/kyle/Drive/ResearchData/Mar09/rasterPre.txt' u i:1 t '' w p lc rgb 'blue' ps 0.5 pt 7

