set style data lines
unset autoscale
set xrange [0:1]
set yrange [0:1]
do for [i=1:99] { plot sprintf('anLin-%d.txt', i) using 1:2; pause 0.5 }
pause -1
