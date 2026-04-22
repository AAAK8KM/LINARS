# Set the input file format to CSV
set datafile separator ","

# Labels and Title
set xlabel "time, ms"
set ylabel "Precision"
set title "Error rate vs time"

# Visual styling
set grid
set key top right
set style data points
set pointsize 1.2
set autoscale
set terminal png
set logscale y
set output 'Impl_TimeVsErr.png' 

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'CGD.csv'  using 5:4 smooth unique title "CGD" with points pt 7 lc rgb "red", \
     'sor.csv'  using 5:4 smooth unique title "SOR" with points pt 7 lc rgb "yellow", \
     'SteeperGD.csv'  using 5:4 smooth unique title "SteeperGD" with points pt 6 lc rgb "blue", \
     'gauszeidel.csv'  using 5:4 smooth unique title "GausZeidel" with points pt 6 lc rgb "green", \
     'ssorchb.csv'  using 5:4 smooth unique title "SSOR+Cheb" with points pt 6 lc rgb "purple", \
     'gmres.csv'  using 5:4 smooth unique title "GMRES(4)" with points pt 2 lc rgb "black"
