rm benchmarks.dat
rm freezeout_max_radius.dat
make clean
make

for i in {1..15}
do
	echo "\n\n\nRunning event $i\n"
	./cpu-vah
done
