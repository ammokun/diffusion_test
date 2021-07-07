set y2tics
set logscale y2
set format y2 "10^{%L}"
set key outside

set y2label "error sqrt((q-qt)^2"
plot "output.dat" u 1:2 title"q0"
replot "output.dat" u 1:3 title"q"
replot "output.dat" u 1:4 title"qt" w l lw 2
replot "output.dat"  u 1:5 axis x1y2  title"error"