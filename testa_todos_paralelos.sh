#!/bin/bash

for TAM in 512 1024 2048 4096; do
    for num_threads in 1 2 4 6 8 10 12 14 16; do
        ./paralelo_omp_inteligente $TAM 8 $num_threads
    done;
    for num_threads in 1 2 4 6 8 10 12 14 16; do
        ./paralelo_omp_inteligente_strassen $TAM 8 $num_threads
    done;
    ./paralelo_acc_inteligente $TAM 
    ./paralelo_acc_inteligente_strassen $TAM 
done;
