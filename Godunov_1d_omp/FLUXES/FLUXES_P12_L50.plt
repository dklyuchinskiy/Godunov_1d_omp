reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1]
set yrange[0.0:0.3]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [0:0.3]

set output 'FLUXES_P12_L_N900_color40pt1v2.png'
splot 'FLUXES_L_1_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_2_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_3_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_4_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_6_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_7_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_8_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_10_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_11_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_12_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_14_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_15_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_16_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_18_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_19_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_20_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_22_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_23_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_24_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_26_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_27_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_28_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_30_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_31_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_32_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_34_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_35_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_36_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_38_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_39_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_40_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_42_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_43_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_44_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_46_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_47_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_49_P12_N900.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
