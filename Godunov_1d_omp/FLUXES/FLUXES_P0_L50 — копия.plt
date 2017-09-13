reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1]
set yrange[0.0:0.6]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [0:0.3]

set output 'FLUXES_P0_L_N900_color40pt1v2_copy.png'
splot 'FLUXES_L_1_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
