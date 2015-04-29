set xrange [0:1.0]
set yrange [0.:1.01]
do for [ii=0:99] {
    plot sprintf('../del_0.05nlin%d.txt',ii) w l t "nlin" ;
    pause 0.2
}
pause -1

#, sprintf('vlin%d.txt',ii) w l t "vlin", sprintf('nlin%d.txt',ii) w l t "nlin"
