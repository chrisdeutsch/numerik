set terminal epslatex color colortext size 3.5, 2.1

set output 'ergebnis.tex'

set title ''

set xlabel 'Anzahl hinzugef\"{u}gter Schalen n'
set xrange [3:17]
set xtics 0, 5
set mxtics 5

set ylabel 'Madelungkonstante $\tilde{\alpha}(n)$'
set yrange [1.74745:1.74775]

set tics scale 0.7

set key right

a(x)=1.7475646946331822

plot "data.csv" using 1:2 lt 3 t'\footnotesize Madelung-Konstante', a(x) t'\footnotesize Literaturwert' lt 1

unset output