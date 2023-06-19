set term post color
set title "Non-Linear Convection"
set output 'Heaviside_q2.ps'
p './heaviside.dat' u 1:2 w l t 'Heaviside at t=0',\
  './scheme_NLC_1.dat' u 1:2 w l t 'FTCS',\
  './scheme_NLC_2.dat' u 1:2 w l t 'FTBS',\
  './scheme_NLC_3.dat' u 1:2 w l t 'Lax',\
  './scheme_NLC_4.dat' u 1:2 w l t 'Lax-Wendroff',\
  './scheme_NLC_5.dat' u 1:2 w l t 'Beam Warming',\
  './scheme_NLC_6.dat' u 1:2 w l t 'MacCormak'
