#!/bin/bash

#==== subroutine ====
function make_fig
{

name=$1
L=$2

gnuplot << EOF

set style line  1 lt 1 pt 4  lc rgbcolor "#ff0200"
set style line  2 lt 1 pt 6  lc rgbcolor "#ff8a1b"
set style line  3 lt 1 pt 8  lc rgbcolor "#267e17"
set style line  4 lt 1 pt 10 lc rgbcolor "#0021d2"
set style line  5 lt 1 pt 12 lc rgbcolor "#8b2bd3"
set style line  6 lt 1 pt 14 lc rgbcolor "#2b3230"

do for [i=1:20] {
  set style line i linewidth 2
}

set pointsize 1
set bar 1

set term postscript eps enhanced color
#set term pdf enhanced color
set size 0.5,0.4
set xlabel "U/t" offset 0,0.5
set ylabel "[E/L^2 - (U+6V)]/t" offset 1,-0.5
set key samplen 2
set key width -1.5
set key left top

set output "fig_${name}.eps"
#set output "fig_${name}.pdf"

set label 10 "(a)" at graph -0.22,0.95

p \
"< grep -v '#' ../dat/fig_4/dat_normal_metal/catdat_L${L}UpV0" u 1:(\$4-\$1-\$2*6) w lp ls 1 ti "normal metal", \
"< grep -v '#' ../dat/fig_4/dat_AF_insulator/catdat_L${L}UpV0" u 1:(\$4-\$1-\$2*6) w lp ls 2 ti "AF insulator", \
"< grep -v '#' ../dat/fig_4/dat_CO_AF_insulator/catdat_L${L}UpV0" u 1:(\$4-\$1-\$2*6) w lp ls 3 ti "AF+CO insulator", \
"< grep -v '#' ../dat/fig_4/dat_FM_metal/catdat_L${L}UpV0" u 1:(\$4-\$1-\$2*6) w lp ls 4 ti "FM metal"

EOF

#----

}

#==== run ====
#make_fig 4a 12
make_fig 4a 24
#make_fig 4a 120

rm -f *~
