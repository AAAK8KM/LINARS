# Set the input file format to CSV
set datafile separator " "

# Labels and Title
set xlabel "tau"
set ylabel "N"
set title "Error rate vs time"

# Visual styling
set grid
set key top right
set style data points
set pointsize 1.2
set autoscale
set terminal png
set logscale y
set output 'NoverTau4.png' 

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'tau4.csv'  using 1:2 smooth unique title "Simple" with points pt 7 lc rgb "red"
