#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double** mult( size_t matrix_tam , double** A , double** B ){
    int linha,coluna,k;
    
    // Alocando a matriz de retorno
    double** C = malloc( sizeof(double*) * matrix_tam );
    for( linha = 0 ; linha < matrix_tam ; linha++ ){
        C[linha] = malloc( sizeof(double) * matrix_tam );
    }

    for( linha = 0 ; linha < matrix_tam ; linha++ ){
        for( coluna = 0 ; coluna < matrix_tam ; coluna++ ){
            C[linha][coluna] = 0;
            for( k = 0 ; k < matrix_tam ; k++ ){
                C[linha][coluna] += A[linha][k] * B[k][coluna] ;
            }
        }
    }
    return C;
}

void printa_matriz( size_t matrix_tam , double** A ){
    int lin,col;
    for( lin = 0 ; lin < matrix_tam ; lin++ ){
        for( col = 0 ; col < matrix_tam ; col++ ){
            printf("%.2f " , A[lin][col]);
        }
        printf("\n");
    }
}

int main( int argc , char** argv ){

    int matrix_tam = 1024;
    if( argc >= 2 ){
        matrix_tam = atoi( argv[1] );
    }
    
    int lin,col;

    // Declarando as matrizes A e B separadamente para terem espaços separados na memória
    double** A = malloc( sizeof(double*)*matrix_tam );
    for( lin = 0 ; lin < matrix_tam ; lin++ ){
        double* nova_linha = malloc( sizeof(double)*matrix_tam );
        for( col = 0 ; col < matrix_tam ; col++ ){
            nova_linha[col] = ( (float)rand() )/RAND_MAX;
            // nova_linha[col] = lin; // Facil de checar corretude
        }
        A[lin] = nova_linha;
    }
    double** B = malloc( sizeof(double*)*matrix_tam );
    for( lin = 0 ; lin < matrix_tam ; lin++ ){
        double* nova_linha = malloc( sizeof(double)*matrix_tam );
        for( col = 0 ; col < matrix_tam ; col++ ){
            nova_linha[col] = ( (float)rand() )/RAND_MAX;
            // nova_linha[col] = col; // Facil de checar corretude
        }
        B[lin] = nova_linha;
    }

    double initial_t = omp_get_wtime();
    double** C = mult( matrix_tam , A , B );
    double final_t = omp_get_wtime();

    printf("Duracao da multiplicacao(TAM=%d): %.3fs\n" , matrix_tam , final_t - initial_t );
    
    if( matrix_tam <= 16 ){
        printf("A:\n");
        printa_matriz(matrix_tam , A);
        printf("B:\n");
        printa_matriz(matrix_tam , B);
        printf("C:\n");
        printa_matriz(matrix_tam , C);
    }

}
