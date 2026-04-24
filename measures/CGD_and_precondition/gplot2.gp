# Set the input file format to CSV
set datafile separator ","

# Labels and Title
set xlabel "n"
set ylabel "Precision"
set title "Error rate vs steps"

# Visual styling
set grid
set key top right
set style data points
set pointsize 1.2
set autoscale
set terminal png
set logscale y
set output 'Impl_NumVsErr.png' 

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'CGD.csv'  using 1:4 smooth unique title "no precondition" with points pt 7 lc rgb "red", \
     'PCCGD.csv'  using 1:4 smooth unique title "precondition" with points pt 6 lc rgb "blue"
