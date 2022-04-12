#!/bin/bash
rm -f *.vtk l2.out
date>diagnostics.txt
mpirun -np 24 ./mpiflow >> diagnostics.txt
date>>diagnostics.txt
