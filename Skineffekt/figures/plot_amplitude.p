set terminal epslatex color colortext size 3.5, 2.1

set output 'amplitude.tex'

set title ''

set xlabel 'Radius $\rho$ [cm]'
set xrange [*:*]
#set xtics 0, 5
#set mxtics 5

set ylabel 'Amplitude von $j$ [Fr/s/cm$^2$]'
set yrange [*:*]
#set ytics 1.61, 0.001

set tics scale 0.7

set key left

plot "../plots/e3.txt" using 1:2 smooth csplines lt 1 t'Amplitude'

unset output