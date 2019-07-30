#!/bin/bash

#==== subroutine ====
function make_fig
{

name=$1
V2=$2
Uc=$3
Uc2=$4

gnuplot << EOF

set style line  1 lt 1 pt 4  lc rgbcolor "#ff0200"
#set style line  2 lt 1 pt 6  lc rgbcolor "#ff8a1b"
set style line  2 lt 1 pt 6  lc rgbcolor "#d67519"
set style line  3 lt 1 pt 8  lc rgbcolor "#267e17"
set style line  4 lt 1 pt 10 lc rgbcolor "#0021d2"
set style line  5 lt 1 pt 12 lc rgbcolor "#8b2bd3"
set style line  6 lt 1 pt 14 lc rgbcolor "#2b3230"
set style line  7 lt 1 pt 2  lc rgbcolor "#000000"

# https://oku.edu.mie-u.ac.jp/~okumura/stat/colors.html
set style line 101 lt 1 pt 4  lc rgbcolor "#ffd1d1" # red
set style line 102 lt 1 pt 4  lc rgbcolor "#ffff99" # yellow
set style line 103 lt 1 pt 4  lc rgbcolor "#cbf266" # yellow-green
set style line 104 lt 1 pt 4  lc rgbcolor "#b4ebfa" # blue
set style line 105 lt 1 pt 4  lc rgbcolor "#edc58f" # orange
set style line 106 lt 1 pt 4  lc rgbcolor "#87e7b0" # green
set style line 107 lt 1 pt 4  lc rgbcolor "#c7b2de" # purple

do for [i=1:100] {
  set style line i linewidth 2
}

set pointsize 0.65
set bar 0.65

#set term pdf enhanced color
set term postscript eps enhanced color
#set size ratio 0.6
set size 0.4,0.35
set xlabel "U/t" offset 0.0,0.5
set ylabel "m^{/Symbol a}" offset 2.0,0.0
set key samplen 1.5
#set key width -3
#set key left top
#set key right top
#set key right bottom
set key right top at graph 1.3,0.9
set rmargin 6

#set xrange [0:10]
set xrange [0:8]
set yrange [-0.55:0.55]
set ytics 0.2
set mytics 2

#set label 1 "" at graph 0.03,0.9
set label 1 "V_2/t=${V2}" at graph 0.6,0.75

set label 2 "(b)" at graph -0.25,0.99

#set label 11 "220200" at graph 0.04,0.77
#set label 12 "CDW" at graph 0.06,0.65
#set label 13 "001122" at graph 0.04,0.32
#set label 14 "CDW" at graph 0.06,0.20
#set label 15 "001122" at graph 0.53,0.32
#set label 16 "CDW+AF" at graph 0.50,0.20

#set output "fig_${name}.pdf"
set output "fig_${name}.eps"

set tics front
set label 1 front
set label 2 front
#set label 11 front
#set label 12 front
#set label 13 front
#set label 14 front
#set label 15 front
#set label 16 front

#set style arrow 1 size character 1,20 filled linewidth 3
#set arrow 1 from first 1.0,1.25 to 0.12,1.1 arrowstyle 1
#set arrow 1 front

set samples 10000
x1=${Uc}
x2=${Uc2}
ymin=-0.55
ymax=0.55
p \
(x<x1)?ymax:ymin w filledcurves x1 ls 102 ti "", \
(x>x1&&x<x2)?ymax:ymin w filledcurves x1 ls 101 ti "", \
(x>x2)?ymax:ymin w filledcurves x1 ls 106 ti "", \
"< awk '\$1>${Uc}{print \$0}' ../dat/fig_5/dat_V4.00" u 1:(-\$8):(0) w e ls 1 ti "A", \
"" u 1:(-\$7):(0) w e ls 2 ti "B", \
"" u 1:(-\$6):(0) w e ls 3 ti "C", \
"" u 1:(-\$5):(0) w e ls 4 ti "D", \
"" u 1:(-\$4):(0) w e ls 5 ti "E", \
"" u 1:(-\$9):(0) w e ls 6 ti "F", \
\
"< awk '\$1<${Uc}{print \$0}' ../dat/fig_5/dat_V4.00" u 1:(-\$5):(0) w e ls 1 ti "", \
"" u 1:(-\$6):(0) w e ls 2 ti "", \
"" u 1:(-\$7):(0) w e ls 3 ti "", \
"" u 1:(-\$8):(0) w e ls 4 ti "", \
"" u 1:(-\$9):(0) w e ls 5 ti "", \
"" u 1:(-\$4):(0) w e ls 6 ti ""

EOF

#----

}

#==== run ====
make_fig 5b 4 0.125 2.375

rm -f *~
