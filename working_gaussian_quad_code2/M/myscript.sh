#!/usr/bin/env bash
rm a.out
ifort gaussquad.f90
./a.out
rm -r running
mkdir running

ifort -c Tulli13_1990.f90
ifort Tulli13_1990.f90 progTulli12_1990.f90 -mkl

cp fort.23 running
cp a.out running
cp parallel_script_v4.py running
cp raw_x.txt running
cp raw_w.txt running
cp job.sh running
cd running

python3 parallel_script_v4.py 



