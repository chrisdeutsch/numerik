set terminal epslatex color colortext size 8, 3

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

plot "../plots/e3.txt" using 1:2 smooth csplines lt 1 t'$\omega=10^3 1$/s', '../plots/e4.txt' using 1:2 smooth csplines lt 2 t'$\omega=10^4 1$/s', '../plots/e5.txt' using 1:2 smooth csplines lt 3 t'$\omega=10^5 1$/s', '../plots/e6.txt' using 1:2 smooth csplines lt 4 t'$\omega=10^6 1$/s'




unset output