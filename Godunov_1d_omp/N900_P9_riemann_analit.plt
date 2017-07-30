reset
					  				  clear

									  				   ###################################
													   				   set term png font "Times-Roman, 16"

																	   				    ##################################

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0050.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0050.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0050.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0100.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0100.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0100.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0150.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0150.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0150.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0200.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0200.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0200.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0250.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0250.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0250.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0300.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0300.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0300.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0350.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0350.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0350.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0400.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0400.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0400.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0450.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0450.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0450.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0500.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0500.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0500.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0550.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0550.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0550.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0600.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0600.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0600.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0650.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0650.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0650.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0700.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0700.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0700.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0750.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0750.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0750.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0800.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0800.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0800.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0850.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0850.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0850.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0900.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0900.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0900.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.0950.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.0950.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.0950.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1000.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1000.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1000.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1050.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1050.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1050.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1100.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1100.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1100.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1150.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1150.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1150.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1200.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1200.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1200.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1250.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1250.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1250.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1300.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1300.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1300.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1350.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1350.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1350.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1400.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1400.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1400.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1450.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1450.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1450.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1500.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1500.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1500.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1550.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1550.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1550.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1600.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1600.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1600.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1650.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1650.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1650.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1700.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1700.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1700.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1750.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1750.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1750.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1800.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1800.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1800.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1850.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1850.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1850.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1900.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1900.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1900.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.1950.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.1950.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.1950.dat' using 1:4 w linespoints pt 7 title "numerical"

set xrange[0:1]

set yrange[0.00:1.10]
set ylabel "pressure"
set output 'workspace/900/P/P_9/N900_P9_P_analit_riemann_0.2000.png'
plot 'workspace/900/N900_P9_SLV0_TERM2_analit_0.2000.dat' using 1:4 w l lw 2 title "exact", 'workspace/900/N900_P9_SLV0_TERM2_L_0.2000.dat' using 1:4 w linespoints pt 7 title "numerical"


clear
reset


