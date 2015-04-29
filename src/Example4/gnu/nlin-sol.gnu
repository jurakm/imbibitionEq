set xlabel "distance"
set ylabel "Saturation"
set key center top
#set grid
set xrange [0:1]
#set xrange [0:0.02]
#set terminal postscript eps color solid lw 3
# set output "clin_sol.eps"
# set output "clin_sol-zoom.eps"



plot "../del_0.05nlin10.txt" u 1:2 w l t "t=0.1",\
     "../del_0.05nlin20.txt" u 1:2 w l t "t=0.2",\
     "../del_0.05nlin40.txt" u 1:2 w l t "t=0.4",\
     "../del_0.05nlin60.txt" u 1:2 w l t "t=0.6",\
     "../del_0.05nlin80.txt" u 1:2 w l t "t=0.8",\
     "../del_0.05nlin100.txt" u 1:2 w l t "t=1.0"

pause -1
