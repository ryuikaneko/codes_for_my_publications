#!/bin/bash

#==== subroutine ====
function make_fig
{

name=$1

gnuplot << EOF

set style line  1 lt 1 pt 4  lc rgbcolor "#ff0200"
set style line  2 lt 1 pt 6  lc rgbcolor "#d67519"
set style line  3 lt 1 pt 8  lc rgbcolor "#267e17"
set style line  4 lt 1 pt 6  lc rgbcolor "#0021d2"
set style line  5 lt 1 pt 12 lc rgbcolor "#8b2bd3"
set style line  6 lt 1 pt 14 lc rgbcolor "#2b3230"
set style line  7 lt 1 pt 2  lc rgbcolor "#000000"
set style line  8 lt 2 pt 2  lc rgbcolor "#000000" dt (3,3)

do for [i=1:100] {
  set style line i linewidth 2
}

set pointsize 0.65
set bar 0.65

#set term pdf enhanced color
set term postscript eps enhanced color
#set size ratio 0.6
set size 0.4
set xlabel "V_2/t" offset 0.0,0.5
set ylabel "Energy" offset 1.5,0.0
set key samplen 1.5
set key width -2
#set key width -4
set key left top
#set key right top
#set key right bottom

#set xrange [0:10]
set xrange [0:3]
set yrange [:]
set ytics 2
set mytics 2

set label 1 "U/t=1" at graph 0.65,0.2

#set output "fig_${name}.pdf"
set output "fig_${name}.eps"

p \
"< awk '\$2>1.25{print \$0}' ../dat/fig_7/catdat_L24U1.00_TMI | grep -v '#'" u 2:4 w p ls 4 ti "TMI state", \
"< awk '\$2>=1.00{print \$0}' ../dat/fig_7/catdat_L24U1.00_220200CDW" u 2:4 w p ls 1 ti "Ground state", \
"< awk '\$2<1.00{print \$0}' ../dat/fig_7/catdat_L24U1.00_metal" u 2:4 w p ls 1 ti ""

EOF

#----

#ps2pdf -dEPSCrop fig_${name}.eps fig_${name}.pdf
#rm fig_${name}.eps

}

#==== run ====
make_fig 7

rm -f *~
