set xlabel "distance"
set ylabel "Saturation"
set key center top
#set grid
#set xrange [0:1]
set xrange [0:0.02]
set terminal postscript eps color solid lw 3
# set output "clin_sol.eps"
 set output "clin_sol-zoom.eps"



plot "clin10.txt" u 1:2 w l t "t=0.1",\
     "clin20.txt" u 1:2 w l t "t=0.2",\
     "clin40.txt" u 1:2 w l t "t=0.4",\
     "clin60.txt" u 1:2 w l t "t=0.6",\
     "clin80.txt" u 1:2 w l t "t=0.8",\
     "clin100.txt" u 1:2 w l t "t=1.0"

pause -1
