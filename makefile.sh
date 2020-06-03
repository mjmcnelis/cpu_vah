
if [ $1 = "icpc" ]; then
	echo "Using icpc compiler"
    cd makefiles
    cp Makefile_icpc ../Makefile
else
    echo "Using gcc compiler"
    cd makefiles
    cp Makefile_gcc ../Makefile
fi


