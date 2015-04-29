set xlabel "time"
set ylabel "Nonwetting source"
set key left center
#set grid
#set logscale x 10
#set yrange [0:1]
set xrange [0:1]
 set terminal postscript eps color solid lw 3
 set output "clin_Q.eps"

plot "clin-flux.txt" u 1:2   t "volume integral",\
     "clin-flux.txt" u 1:3   t "boundary integral",\
     "flux-anal.txt" u 1:2 w l t "analytic"

pause -1
