reset
clear

set term png font "Times-Roman, 16"
set view map

set xrange[0:1]
set yrange[0.0:0.3]

#set xrange[0:20.0]
#set yrange[0:4]

set cbrange [0:0.3]

set output 'FLUXES_P2_L_N300_color40pt1v2_lagrang.png'
splot 'FLUXES_L_1_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_2_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_3_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_4_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_5_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_6_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_7_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_8_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_9_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_10_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_11_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_12_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_13_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_14_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_15_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_16_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_17_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_18_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_19_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_20_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_21_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_22_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_23_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_24_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_25_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_26_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_27_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_28_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_29_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_30_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_31_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_32_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_33_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_34_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_35_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_36_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_37_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_38_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_39_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_40_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_41_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_42_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_43_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_44_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_45_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_46_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_47_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_48_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle,\
'FLUXES_L_49_P2_N300.dat' u 1:2:3 w points pt 7 palette pointsize 1 notitle
