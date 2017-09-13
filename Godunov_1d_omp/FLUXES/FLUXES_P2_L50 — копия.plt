reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1.0]
set yrange[0.0:0.3]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [0:0.3]

set output 'FLUXES_P2_L_N900_color40pt1v2_copy.png'
splot 'FLUXES_L_1_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P2_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
