set terminal epslatex size 10cm,9cm color colortext

set output 'heatmap_a_20.tex'

set title 'Heat Map f{\"u}r Gitterabstand $a=0,2$'
unset key
set tic scale 0

# colors
set palette defined (0 "white", 0.01 "blue", 0.37 "red")
set cbrange [0.2:*]
set cblabel "Temperaturverteilung"
set cbtics ("$0$" 0.2, "$0.22$" 0.22, "$0.37$" 0.369)

set xrange [*:*]
set yrange [*:*]
set xtics ("$0$" 0, "$1$" 25/5, "$2$" 2*25/5, "$3$" 3*25/5, "$4$" 4*25/5, "$5$" 5*25/5)
set ytics ("$0$" 0, "$1$" 20/4, "$2$" 2*20/4, "$3$" 3*20/4, "$4$" 4*20/4)

set view map
splot '../Plotdaten/heatmap_a_20_gnu.txt' matrix with image

unset output