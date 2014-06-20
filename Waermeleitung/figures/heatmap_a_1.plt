set terminal epslatex size 10cm,9cm color colortext

set output 'heatmap_a_1.tex'

set title 'Heat Map f{\"u}r Gitterabstand $a=0,01$'
unset key
set tic scale 0

# colors
set palette defined (0 "white", 0.01 "blue", 0.243 "red")
set cbrange [0.215:*]
set cblabel "Temperaturverteilung"
set cbtics ("$0.22$" 0.22, "$0.24$" 0.24)

set xrange [*:*]
set yrange [*:*]
set xtics ("$0$" 0, "$1$" 500/5, "$2$" 2*500/5, "$3$" 3*500/5, "$4$" 4*500/5, "$5$" 5*500/5)
set ytics ("$0$" 0, "$1$" 400/4, "$2$" 2*400/4, "$3$" 3*400/4, "$4$" 4*400/4)

set view map
splot '../Plotdaten/heatmap_a_1_gnu.txt' matrix with image

unset output