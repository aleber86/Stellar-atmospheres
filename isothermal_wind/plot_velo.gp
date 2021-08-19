#set terminal postscript color enhanced
#set output "vel_curve.ps"
set key off
set xlabel "r/r_{c}" font ", 14"
set ylabel "v/v_{c}" font ", 14"
set xrange[0:10]
set yrange[0:5]
plot "v-vs-r.dat" index 1 u ($2<=1&&$1<=1?$1:1/0):2 with lines lw 2 lt rgb "blue",\
"v-vs-r.dat" index 0 u ($2>=1&&$1>=1?$1:1/0):2 with lines lw 2 lt rgb "blue",\
"v-vs-r.dat" index 0 u ($2>=1&&$1<=1?$1:1/0):2 with lines lw 1 lt rgb "red",\
"v-vs-r.dat" index 1 u ($2<=1&&$1>=1?$1:1/0):2 with lines lw 1 lt rgb "red",\
for [i=4:99] "v-vs-r.dat" index i u 1:2 with lines lw 1 lt rgb "black"
