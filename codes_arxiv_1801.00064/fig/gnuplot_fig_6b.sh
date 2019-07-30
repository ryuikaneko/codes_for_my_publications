#!/bin/bash

#==== subroutine ====
function make_fig
{

name=$1

gnuplot << EOF

set style fill transparent solid 0.5 noborder

set style line  1 lt 1 pt 4  lc rgbcolor "#ff0200" lw 2
set style line  2 lt 1 pt 6  lc rgbcolor "#0021d2" lw 2
set style line 10 lt 1 pt 2  lc rgbcolor "#000000" lw 1

set term postscript eps enhanced color
#set size 0.6
#set size 0.4,0.35
set size 0.4,0.3
set xlabel "Energy" offset 0,0.5
set ylabel "DOS" offset 2,0
set key left bottom
#set key invert
set key samplen 2
#set key width -3

set xrange [-6.5:6.5]
set yrange [0:1.5]
set ytics 0.5
set tics front

set output "fig_${name}.eps"

#set label 1 "{/*.8 U/t=2, V/t=2, 001122CDW}" at graph 0.15,0.85
set label 1 "{/*.8 U/t=2, V/t=2, 001122CDW}" at graph 0.12,0.85

set label 2 "(b)" at graph -0.25,0.99

p \
"../dat/fig_6/dat_L120U2.00V2.00_001122CDW" i 2 u (\$2-1.68299083471737987/2):4 w filledcurve x1 ls 1 ti "", \
"../dat/fig_6/dat_L120U2.00V2.00_001122CDW" i 2 u (\$2-1.68299083471737987/2):4 w l ls 1 ti ""

EOF

#----

}

#==== run ====
make_fig 6b

rm -f *~
