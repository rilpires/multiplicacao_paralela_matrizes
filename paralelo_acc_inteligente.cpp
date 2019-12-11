#include "matriz_io.cpp"
#include <omp.h>

/**
 * Compilar com:
 *      pg++ -acc paralelo_acc_inteligente.cpp -o paralelo_acc_inteligente
 * 
 * Executar com:
 *      ./paralelo_acc_inteligente [TAMANHO_MATRIZ]
**/

typedef double FLOAT;

FLOAT* gera_transposta( size_t matrix_tam , FLOAT* A ){
    FLOAT initial_t = omp_get_wtime();
    FLOAT* ret = (FLOAT*)malloc( sizeof(FLOAT)*matrix_tam*matrix_tam );
    size_t matrix_vetor_tam = matrix_tam*matrix_tam;
    int i , j;

    #pragma acc data copyin(A[:matrix_vetor_tam]) copyout(ret[:matrix_vetor_tam])
    #pragma acc parallel loop collapse(2)
    for( j = 0 ; j < matrix_tam ; j++ ){
        for( i = 0 ; i < matrix_tam ; i++ ){
            ret[i+j*matrix_tam] = A[j+i*matrix_tam];
        }
    }

    FLOAT final_t = omp_get_wtime();
    return ret;
}

FLOAT* mult( size_t matrix_tam , FLOAT* A , FLOAT* B ){
    FLOAT*  C = (FLOAT*)malloc( sizeof(FLOAT) * matrix_tam * matrix_tam );
    int     linha,coluna,i;
    FLOAT   temp_sum;
    FLOAT*  BT = gera_transposta(matrix_tam,B);
    size_t  matrix_vetor_tam = matrix_tam*matrix_tam;

    #pragma acc data \
        copyin(A[:matrix_vetor_tam],BT[:matrix_vetor_tam]) \
        copyout(C[:matrix_vetor_tam])
    {
        #pragma acc parallel loop independent tile(8,8)
        for( linha = 0 ; linha < matrix_tam ; linha++ )
        for( coluna = 0 ; coluna < matrix_tam ; coluna ++ )
        {
            for( temp_sum = 0 , i = 0 ; i < matrix_tam ; i ++ )
                temp_sum += A[i + linha*matrix_tam] * BT[i + coluna*matrix_tam];
            C[coluna + linha*matrix_tam] = temp_sum;
        }
    }
    free( BT );
    return C;
}

int main( int argc , char** argv ){

    int matrix_tam = 1024 ;
    if( argc >= 2 ) matrix_tam = atoi( argv[1] );

    char nome_matriz_a[64] , nome_matriz_b[64] , nome_resultado[64];
    sprintf( nome_matriz_a , "A_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_matriz_b , "B_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_resultado , "resultado_%d_%d.txt" , matrix_tam , matrix_tam );

    FLOAT* A = le_matriz_vetor<FLOAT>( nome_matriz_a , matrix_tam );
    FLOAT* B = le_matriz_vetor<FLOAT>( nome_matriz_b , matrix_tam );

    FLOAT initial_t = omp_get_wtime();
    FLOAT* C = mult( matrix_tam , A , B );
    FLOAT final_t = omp_get_wtime();

    printf("[ACC] - Duracao da multiplicacao (TAM=%d): %.3fs\n" , matrix_tam , final_t - initial_t );

    if( matrix_tam <= 1024 )
        escreve_matriz_sequencial( nome_resultado , matrix_tam , C );

    
    printf("[Corretude]\nOs 3 primeiros valores de C:\n" );
    printf("C: %f %f %f \n" , C[0] , C[1] , C[2] );
    
}
