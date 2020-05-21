
if [ $1 = "icpc" ]; then
	echo "Using icpc compiler"
    cp Makefile_icpc Makefile
else
    echo "Using gcc compiler"
    cp Makefile_gcc Makefile
fi

