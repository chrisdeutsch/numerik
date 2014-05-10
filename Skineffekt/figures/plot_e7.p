set terminal epslatex color colortext size 5, 3

set output 'e7.tex'

set title ''

set xlabel 'Radius $\rho$ [cm]'
set xrange [*:*]

set ylabel '$|j_0(\rho)|$ [Fr/s/cm$^2$]'
set yrange [0:*]

set y2label '$\phi(\rho)$ [rad]' 
set y2range[-2*pi:2*pi]
set y2tics ('$-\pi$' -pi, $0$ 0, '$\pi$' pi)
set ytics nomirror

set tics scale 0.7

set key left top

plot "../plots/e7.txt" using 1:2 smooth csplines lt 1 axes x1y1 t'\footnotesize Amplitude', "../plots/e7.txt" using 1:3 with lines lt 2 axes x1y2 t'\footnotesize Phase' 

unset output