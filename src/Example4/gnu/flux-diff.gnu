set xlabel "time"
set ylabel "Nonwetting source"
set key  right  top
#set title "Y(t) = 0.5 + 0.20* sin(2*PI*t), delta=0.05"
set title "imbibition"
set grid
set xrange [0:1]
set yrange [-1:1]
#set terminal postscript eps color solid lw 3
#set output "diff_Q_Ex4-0.05-3.0.eps"

plot "../del_0.05clin-flux.txt" u 1:2 w l  t "constant linearized",\
     "../del_0.05vlin-flux.txt" u 1:2 w l  t "variable linearized",\
     "../del_0.05nlin-flux.txt" u 1:2 w l  t "nonlinear"

#     "../del_0.05nlin-flux.txt" u 1:($4/100)  t "boundary value/100"

pause -1
