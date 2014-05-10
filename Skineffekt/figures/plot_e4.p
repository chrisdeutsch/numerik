set terminal epslatex color colortext size 5, 3

set output 'e4.tex'

set title ''

set xlabel 'Radius $\rho$ [cm]'
set xrange [*:*]

set ylabel 'Amplitude von $j$ [Fr/s/cm$^2$]'
set yrange [0:*]

set y2label 'Phase von $j$ [rad]' 
set y2range[*:*]
set y2tics
set ytics nomirror

set tics scale 0.7

set key right bottom

plot "../plots/e4.txt" using 1:2 smooth csplines lt 1 axes x1y1 t'\footnotesize Amplitude', "../plots/e4.txt" using 1:3 smooth csplines lt 2 axes x1y2 t'\footnotesize Phase' 

unset output