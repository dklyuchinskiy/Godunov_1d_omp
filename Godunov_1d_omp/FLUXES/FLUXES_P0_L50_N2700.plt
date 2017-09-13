reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1]
set yrange[0.0:0.6]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [0:0.3]

set output 'FLUXES_P0_L_N900_color40pt1v2_lagrange.png'
splot 'FLUXES_L_1_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_2_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_3_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_4_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_6_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_7_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_8_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_10_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_11_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_12_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_14_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_15_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_16_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_18_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_19_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_20_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_22_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_23_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_24_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_26_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_27_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_28_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_30_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_31_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_32_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_34_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_35_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_36_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_38_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_39_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_40_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_42_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_43_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_44_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_46_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_47_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_49_P0_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
