set terminal epslatex size 10cm,9cm color colortext

set output 'heatmap_a_10.tex'

set title 'Heat Map f{\"u}r Gitterabstand $a=0,1$'
unset key
set tic scale 0

# colors
set palette defined (0 "white", 0.01 "blue", 0.29 "red")
set cbrange [0.2:*]
set cblabel "Temperaturverteilung"
set cbtics ("$0.22$" 0.22, "$0.29$" 0.29)

set xrange [*:*]
set yrange [*:*]
set xtics ("$0$" 0, "$1$" 51/5, "$2$" 2*51/5, "$3$" 3*51/5, "$4$" 4*51/5, "$5$" 5*50/5)
set ytics ("$0$" 0, "$1$" 41/4, "$2$" 2*41/4, "$3$" 3*41/4, "$4$" 4*40/4)

set view map
splot '../Plotdaten/heatmap_a_10_gnu.txt' matrix with image

unset output