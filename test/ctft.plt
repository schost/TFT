set xlabel "d"
set ylabel "time"
set term post eps enhanced color

set out "ctft.eps"
set style line 1 lt 1 lw 3 lc 3
set style line 2 lt 1 lw 3 lc 7
set style line 3 lt 1 lw 3 lc 1
plot "ctft.dat" using 1:2 with lines ls 1 title "ctft",\
     "ctft.dat" using 1:3 with lines ls 2 title "ntl",\
     "ctft_middle.dat" using 1:2 with lines ls 3 title "middle product"
