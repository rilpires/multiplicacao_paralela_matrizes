#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double* gera_transposta( size_t matrix_tam , double* A ){
    double initial_t = omp_get_wtime();
    double* ret = malloc( sizeof(double)*matrix_tam*matrix_tam );
    int i , j;
    for( j = 0 ; j < matrix_tam ; j++ ){
        for( i = 0 ; i < matrix_tam ; i++ ){
            ret[i+j*matrix_tam] = A[j+i*matrix_tam];
        }
    }
    double final_t = omp_get_wtime();
    printf("Duracao do gera_transposta: %.3fs\n" , final_t-initial_t);
    return ret;
}

double* mult( size_t matrix_tam , double* A , double* B ){
    double* C = malloc( sizeof(double) * matrix_tam * matrix_tam );
    int linha,coluna,k;
    double temp;
    double *C_linha_pointer , *A_linha_pointer , *B_coluna_pointer; 
    double* BT = gera_transposta(matrix_tam,B);

    for( linha = 0 ; linha < matrix_tam ; linha++ ){
        C_linha_pointer = C + linha*matrix_tam;
        A_linha_pointer = A + linha*matrix_tam;
        for( coluna = 0 ; coluna < matrix_tam ; coluna++ ){
            temp = 0;
            B_coluna_pointer = BT + coluna*matrix_tam;
            for( k = 0 ; k < matrix_tam ; k ++ ){
                temp += A_linha_pointer[k] * B_coluna_pointer[k] ;
            }
            C_linha_pointer[coluna] = temp;
        }
    }
    free( BT );
    return C;
}

void printa_matriz( size_t matrix_tam , double* A ){
    int lin,col;
    for( lin = 0 ; lin < matrix_tam ; lin++ ){
        for( col = 0 ; col < matrix_tam ; col++ ){
            printf("%.2f " , A[col + lin*matrix_tam]);
        }
        printf("\n");
    }
}

int main( int argc , char** argv ){

    int matrix_tam = 1024;
    if( argc >= 2 ){
        matrix_tam = atoi( argv[1] );
    }
    
    double* A = malloc( sizeof(double)*matrix_tam*matrix_tam );
    double* B = malloc( sizeof(double)*matrix_tam*matrix_tam );
    for( int i = 0 ; i < matrix_tam*matrix_tam ; i++ ){
        A[i] = ((double)rand()) / RAND_MAX;
        B[i] = ((double)rand()) / RAND_MAX;
        // A[i] = i % matrix_tam; // Facil de checar corretude
        // B[i] = i / matrix_tam; // Facil de checar corretude
    }


    double initial_t = omp_get_wtime();
    double* C = mult( matrix_tam , A , B );
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
