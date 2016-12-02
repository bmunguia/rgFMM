#!/bin/bash -l
#PBS -l nodes=1:ppn=24
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -q gpu
#PBS -d .

regent.py ../FMM.rg -ll:cpu 8 -ll:csize 32768 -i ../test7.dat -o output7.dat -p 8 -hl:prof 1 -logfile log_%
