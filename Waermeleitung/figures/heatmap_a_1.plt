set terminal pngcairo

set output 'heatmap_a_1.png'

set title 'Heat Map f\''u r Gitterabstand $a=0,01$'
unset key
set tic scale 0

# colors
set palette defined (0 "white", 0.21 "white", 0.22 "blue", 0.245 "red")
set cbrange [0.21:*]
set cblabel "Temperaturverteilung"
set cbtics

set xrange [*:*]
set yrange [*:*]
set xtics ("0" 0, "1" 201/5, "2" 2*201/5, "3" 3*201/5, "4" 4*201/5, "5" 5*201/5)
set ytics ("0" 0, "1" 151/4, "2" 2*151/4, "3" 3*151/4, "4" 4*151/4)

set view map
splot 'heatmap_a_1_gnu.txt' matrix with image

unset output