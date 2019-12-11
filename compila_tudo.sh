#!/bin/bash

echo "Compilando sequencial_inteligente.cpp" &&
pgc++ sequencial_inteligente.cpp -o sequencial_inteligente &&
echo "Compilando sequencial_ingenuo.cpp" &&
pgc++ sequencial_ingenuo.cpp -o sequencial_ingenuo &&
echo "Compilando paralelo_omp_inteligente.cpp" &&
pgc++ -mp paralelo_omp_inteligente.cpp -o paralelo_omp_inteligente &&
echo "Compilando paralelo_omp_inteligente_strassen.cpp" &&
pgc++ -mp paralelo_omp_inteligente_strassen.cpp -o paralelo_omp_inteligente_strassen &&
echo "Compilando paralelo_acc_inteligente.cpp" &&
pgc++ -acc -ta=nvidia paralelo_acc_inteligente.cpp -o paralelo_acc_inteligente &&
echo "Compilando paralelo_acc_inteligente_strassen.cpp" &&
pgc++ -acc -ta=nvidia paralelo_acc_inteligente_strassen.cpp -o paralelo_acc_inteligente_strassen &&
echo "Tudo compilado"