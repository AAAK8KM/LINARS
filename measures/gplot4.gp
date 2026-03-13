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
set output 'Impl_NumVsRealErr.png' 

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'simple.csv'  using 1:2 smooth unique title "Simple" with points pt 7 lc rgb "red", \
     'jacobi.csv'  using 1:2 smooth unique title "Jacobi" with points pt 6 lc rgb "blue", \
     'gauszeidel.csv'  using 1:2 smooth unique title "GausZeidel" with points pt 6 lc rgb "green"
