#!/bin/bash

for i in `seq 1 18`
do
	echo "mpiexec -np $i ./test ./testdata/sample"
	mpiexec -np $i ./mpi/test ./testdata/sample
done

for i in `seq 1 30`
do
	echo "mpiexec -np $i ./test ./testdata/Inputfile"
	mpiexec -np $i ./mpi/test ./testdata/Inputfile
done
