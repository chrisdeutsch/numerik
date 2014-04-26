set terminal epslatex color colortext size 3.5, 2.1

set output 'ergebnis2d.tex'

set title ''

set xlabel 'Anzahl hinzugef\"{u}gter Schalen n'
set xrange [3:17]
set xtics 0, 5
set mxtics 5

set ylabel 'Madelung-Konstante $\tilde{\alpha}(n)$'
set yrange [1.610:1.616]
set ytics 1.61, 0.001

set tics scale 0.7

set key right bottom

a(x)=1.6155426267128247

plot "data2d.csv" using 1:2 lt 3 t'\footnotesize Madelung-Konstante $\tilde{\alpha}$', a(x) t'\footnotesize Literaturwert ${\alpha}$' lt 1

unset output