set term pop
set output
set xlabel 'Time Delay (ms)'
set mxtics 10

set ylabel 'Mutual Information(Bits)'
set mytics 5
set format y '%.2te%+T'
set ytics add ('0' 0)

set grid

#set terminal postscript color enhanced eps
#set output './figure-eps/current_draw.eps'
#set size 1,1

plot './file-dat/tdmi_ordered.dat' u 1:2 t 'TDMI ordered' w lp, \
'./file-dat/tdmi_rand.dat' u 1:2 t 'TDMI rand' w lp
