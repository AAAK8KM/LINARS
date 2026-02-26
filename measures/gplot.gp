# Set the input file format to CSV
set datafile separator ","

# Labels and Title
set xlabel "n"
set ylabel "Time, ms"
set title "Average Time for matrix 3000x3000 on vector multiplication"

# Visual styling
set grid
set key top right
set style data points
set pointsize 1.2

# Plot command
# 'every ::1' skips the header line
# 'smooth unique' groups by column 1 (n) and calculates the average of the Y column
plot 'res.csv' every ::1 using 1:2 smooth unique title "Avg Dense matrix" with points pt 7 lc rgb "blue", \
     'res.csv' every ::1 using 1:3 smooth unique title "Avg CSR matrix" with points pt 6 lc rgb "red"