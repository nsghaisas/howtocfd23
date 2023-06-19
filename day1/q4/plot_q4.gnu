 
set term post color
set output 'burger_2d.ps'
set dgrid3d 50,50,50
#set samples 20, 20
#set isosamples 21, 21
#set contour base
#set cntrparam levels 15
set style data lines
set title "2D Burger solution plot" 
set xlabel "X axis" 
set ylabel "Y axis" 
set zlabel "U-Magnitude " 
splot 'burger_2d.dat' using 1:2:5
