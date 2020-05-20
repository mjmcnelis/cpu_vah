cd python
sh generate_samples.sh $1
cd ..

rm -r output
mkdir output
rm -r semi_analytic
mkdir semi_analytic

make clean
make

for ((n = 1; n <= $1; n++))
do
    echo
    echo "Running event $n"
    ./cpu_vah $n
    echo
done

# run a single hydro event per model parameter sample
# $1 = total number of model parameter samples
