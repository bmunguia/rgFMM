#!/bin/bash -l
#PBS -l nodes=1:ppn=24
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q gpu
#PBS -d .

regent.py ../FMM.rg -i ../test7.dat -o output7.dat -p 2 -ll:cpu 2 -ll:csize 32768 -hl:prof 1 -logfile log_%
