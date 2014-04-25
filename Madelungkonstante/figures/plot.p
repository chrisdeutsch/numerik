set terminal postscript eps color enhanced "Helvetica" 16 size 3, 2

set output 'ergebnis.eps'

set title ''

set xlabel 'Anzahl hinzugefuegter Schalen n'
set xrange [3:17]
set xtics 0, 5
set mxtics 5

set ylabel 'Madelungkonstante {/Symbol a}'
set yrange [1.74745:1.74775]

set key right

a(x)=1.7475646946331822

plot "data.csv" using 1:2 lt 3 t"Madelung-Konstante", a(x) t"Literaturwert" lt 1

unset output