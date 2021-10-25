set terminal postscript color enhanced
set output "op_coeff_bf.ps"
set title "Bound-free opacity coefficient for HI" font ", 14"
set xrange[0:14]
set key left
set xlabel "1/{/Symbol l} [{/Symbol m}m^{-1}]" font ", 12"
set ylabel "log_{10}(k_{bf}({/Symbol l}))" font ", 12"
plot  "dataB5V.dat" u 1:(log10($2)) with lines title "B5V",\
"dataG2V.dat" u 1:(log10($2)) with lines title "G2V"
