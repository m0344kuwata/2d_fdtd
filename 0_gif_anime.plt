set term gif animate optimize delay 1 size 1920, 1080
set output 'movie.gif'


set xrange [0:1]
set xlabel "x [-]"
set xtics 0.5

set yrange [0:1]
set ylabel "y [-]"
set ytics 0.5

set zrange [-1:1]
set zlabel "E_Field [-]"
set ztics 0.5

do for[i=10001:10500]{
splot sprintf("./e_field_%d.txt", i)  
}

set out
set terminal wxt enhanced