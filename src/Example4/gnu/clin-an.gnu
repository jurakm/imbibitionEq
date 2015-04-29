set xrange [0:1]
set yrange [0.:1.01]
do for [ii=1:99] {
    plot sprintf('../del_0.05clin%d.txt',ii) w l t "clin" ;
    pause 0.1
}
pause -1

#, sprintf('vlin%d.txt',ii) w l t "vlin", sprintf('nlin%d.txt',ii) w l t "nlin"
