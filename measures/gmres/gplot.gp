# Set the input file format to CSV
set datafile separator ","

# Labels and Title
set xlabel "n"
set ylabel "Time"
set title "Num vs time"

# Visual styling
set grid
set key top right
set style data points
set pointsize 1.2
set autoscale
set terminal png
set logscale y
set output 'NumVsTime.png' 

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'gmresm.csv'  using 1:5 smooth unique title "GMRES(4)" with points pt 6 lc rgb "black"