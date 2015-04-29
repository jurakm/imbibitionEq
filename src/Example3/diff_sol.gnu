set xlabel "distance"
set ylabel "Saturation"
set key center bottom
set xrange [0:0.02]
set yrange [0.14:0.55]

do for [no = 1:100] {

plot sprintf("del_0.05clin%d.txt", no) u 1:2 w l t "c",  sprintf("del_0.05nlin%d.txt",no) u 1:2 w l t "n", sprintf("del_0.05vlin%d.txt", no) u 1:2 w l t "v"

pause 1
}
