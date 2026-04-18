set datafile separator ","
set terminal pdfcairo enhanced color
set output "ex7_strong_scaling.pdf"
set grid
set key outside
set xlabel "MPI processes"
set ylabel "Execution time (s)"
set title "Exercise 7 strong scaling"
plot for [ex=1:6] "ex7_results.csv" using (strcol(1) eq "strong" && int(strcol(2)) == ex ? $3 : 1/0):6 with linespoints title sprintf("ex%d", ex)

set output "ex7_speedup.pdf"
set ylabel "Speedup vs 1 process"
set title "Exercise 7 speedup"
time1(ex) = real(system(sprintf("awk -F, 'BEGIN{v=0} $1==\"strong\" && $2==\"%d\" && $3==\"1\" {v=$6; exit} END{print v}' ex7_results.csv", ex)))
plot for [ex=1:6] "ex7_results.csv" using (strcol(1) eq "strong" && int(strcol(2)) == ex ? $3 : 1/0):(time1(ex)/$6) with linespoints title sprintf("ex%d", ex)
