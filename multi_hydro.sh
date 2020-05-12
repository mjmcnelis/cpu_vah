rm benchmarks.dat

for i in {1..50}
do
	echo "\n\n\nRunning event $i\n"
	sh hydro.sh
done