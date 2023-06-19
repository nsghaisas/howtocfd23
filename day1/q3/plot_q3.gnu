set term post color
set title "1D Burger Equation Solution"
set output 'Heaviside_q3.ps'
p './heaviside.dat' u 1:2 w l t 'Heaviside at t=0',\
  './burger_1d.dat' u 1:2 w l t 'MacCormak',
