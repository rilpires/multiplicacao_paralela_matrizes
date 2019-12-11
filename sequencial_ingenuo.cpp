#include "matriz_io.cpp"
#include <omp.h>

typedef double FLOAT;

FLOAT** mult( size_t matrix_tam , FLOAT** A , FLOAT** B ){
    int linha,coluna,k;
    
    // Alocando a matriz de retorno
    FLOAT** C = (FLOAT**)malloc( sizeof(FLOAT*) * matrix_tam );
    for( linha = 0 ; linha < matrix_tam ; linha++ ){
        C[linha] = (FLOAT*)malloc( sizeof(FLOAT) * matrix_tam );
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

int main( int argc , char** argv ){

    int matrix_tam = 1024;
    if( argc >= 2 )
        matrix_tam = atoi( argv[1] );

    char nome_matriz_a[64] , nome_matriz_b[64] , nome_resultado[64];
    sprintf( nome_matriz_a , "A_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_matriz_b , "B_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_resultado , "resultado_%d_%d.txt" , matrix_tam , matrix_tam );

    FLOAT** A = le_matriz<FLOAT>( nome_matriz_a , matrix_tam );
    FLOAT** B = le_matriz<FLOAT>( nome_matriz_b , matrix_tam );

    FLOAT initial_t = omp_get_wtime();
    FLOAT** C = mult( matrix_tam , A , B );
    FLOAT final_t = omp_get_wtime();

    printf("Duracao da multiplicacao(TAM=%d): %.3fs\n" , matrix_tam , final_t - initial_t );
    
    if( matrix_tam <= 1024 )
        escreve_matriz( nome_resultado , matrix_tam , C );

    printf("[Corretude]\nOs 3 primeiros valores de C:\n" );
    printf("C: %f %f %f \n" , C[0] , C[1] , C[2] );
}
