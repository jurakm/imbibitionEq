set xrange [0:0.001]
set yrange [0.:1.01]
do for [ii=1:99] {
    plot sprintf('../del_0.05vlin%d.txt',ii) w l t "vlin" ;
    pause 0.1
}
pause -1

#, sprintf('vlin%d.txt',ii) w l t "vlin", sprintf('nlin%d.txt',ii) w l t "nlin"
