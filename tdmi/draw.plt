set xlabel 'Delayed (ms)'
set mxtics 10

set ylabel 'Mutual Information(Bits)'
set mytics 5

set grid
set title 'TDMI of #65 neuron and LFP from its first-order connected post-neurons'

set terminal postscript color enhanced eps
set output './figure_eps/current_draw.eps'
set size 1,1

plot './file_dat/tdmi_ordered.dat' u 1:2 t 'TDMI ordered' w lp, \
'./file_dat/tdmi_rand.dat' u 1:2 t 'TDMI rand' w lp
