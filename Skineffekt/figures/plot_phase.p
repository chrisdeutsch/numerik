set terminal epslatex color colortext size 3.5, 2.1

set output 'phase.tex'

set title ''

set xlabel 'Radius $\rho$'
set xrange [*:*]
#set xtics 0, 5
#set mxtics 5

set ylabel 'Phase der Stromdichte $j$'
set yrange [*:*]
#set ytics 1.61, 0.001

set tics scale 0.7

set key left top

plot "data.csv" using 1:3 lt 1 t'Phase ' 

unset output