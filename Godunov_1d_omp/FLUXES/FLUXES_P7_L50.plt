reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1]
set yrange[0.0:0.2]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [-1.0:1.0]

set output 'FLUXES_P7_L_N900_color40pt1v2_lagrange.png'
splot 'FLUXES_L_1_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_2_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_3_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_4_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_6_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_7_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_8_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_10_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_11_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_12_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_14_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_15_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_16_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_18_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_19_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_20_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_22_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_23_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_24_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_26_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_27_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_28_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_30_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_31_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_32_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_34_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_35_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_36_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_38_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_39_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_40_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_42_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_43_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_44_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_46_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_47_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_49_P7_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
