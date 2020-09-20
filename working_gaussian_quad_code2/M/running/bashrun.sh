#!/bin/bash
mkdir 20
cp a.out 20
cp job.sh 20
cp fort.23 20
cp raw_x.txt 20
cp raw_w.txt 20
mv input.txt 20
cd 20
qsub job.sh	
cd .. 
