#include "matriz_io.cpp"
#include <omp.h>

typedef double FLOAT;

#ifndef OPERATIONS_PER_LOOP
    #define OPERATIONS_PER_LOOP 4
#endif

FLOAT* gera_transposta( size_t matrix_tam , FLOAT* A ){
    FLOAT initial_t = omp_get_wtime();
    FLOAT* ret = (FLOAT*)malloc( sizeof(FLOAT)*matrix_tam*matrix_tam );
    int i , j;

    #pragma acc data copyin(A[:matrix_tam*matrix_tam]) copyout(ret[:matrix_tam*matrix_tam])
    #pragma acc parallel loop collapse(2)
    for( j = 0 ; j < matrix_tam ; j++ ){
        for( i = 0 ; i < matrix_tam ; i++ ){
            ret[i+j*matrix_tam] = A[j+i*matrix_tam];
        }
    }

    FLOAT final_t = omp_get_wtime();
    // printf("Duracao do gera_transposta: %.3fs\n" , final_t-initial_t);
    return ret;
}

FLOAT* mult( size_t matrix_tam , FLOAT* A , FLOAT* B , size_t block_tam ){
    FLOAT* C = (FLOAT*)malloc( sizeof(FLOAT) * matrix_tam * matrix_tam );
    int block_x,block_y,linha,coluna,k, num_blocos = matrix_tam/block_tam , lim1 , lim2 ;
    FLOAT temp_sum[8],temp_A;
    FLOAT *C_linha_pointer , *A_linha_pointer;
    FLOAT* B_coluna_pointers[8];
    FLOAT* BT = gera_transposta(matrix_tam,B);

    #pragma acc data copyin(A[:matrix_tam*matrix_tam],B[:matrix_tam*matrix_tam]) copyout(C[:matrix_tam*matrix_tam])
    #pragma acc parallel loop collapse(2)
    {
        for( block_x = 0 ; block_x < num_blocos ; block_x++ )
        for( block_y = 0 ; block_y < num_blocos ; block_y++ )
        {
            lim1 = (block_y+1)*block_tam;
            for( linha = block_y*block_tam ; linha < lim1 ; linha++ ){
                C_linha_pointer = C + linha*matrix_tam;
                A_linha_pointer = A + linha*matrix_tam;
                lim2 = (block_x+1)*block_tam;
                for( coluna = block_x*block_tam ; coluna < lim2 ; coluna += OPERATIONS_PER_LOOP ){
                    temp_sum[0] = 0;
                    B_coluna_pointers[0] = BT + (coluna+0)*matrix_tam;
                    #if OPERATIONS_PER_LOOP >= 2
                    temp_sum[1] = 0;
                    B_coluna_pointers[1] = BT + (coluna+1)*matrix_tam;
                    #if OPERATIONS_PER_LOOP >= 4
                    temp_sum[2] = 0;
                    temp_sum[3] = 0;
                    B_coluna_pointers[2] = BT + (coluna+2)*matrix_tam;
                    B_coluna_pointers[3] = BT + (coluna+3)*matrix_tam;
                    #if OPERATIONS_PER_LOOP == 8
                    temp_sum[4] = 0;
                    temp_sum[5] = 0;
                    temp_sum[6] = 0;
                    temp_sum[7] = 0;
                    B_coluna_pointers[4] = BT + (coluna+4)*matrix_tam;
                    B_coluna_pointers[5] = BT + (coluna+5)*matrix_tam;
                    B_coluna_pointers[6] = BT + (coluna+6)*matrix_tam;
                    B_coluna_pointers[7] = BT + (coluna+7)*matrix_tam;
                    #endif
                    #endif
                    #endif
                    for( k = 0 ; k < matrix_tam ; k ++ ){
                        temp_A = A_linha_pointer[k];
                        temp_sum[0] += temp_A * B_coluna_pointers[0][k] ;
                        #if OPERATIONS_PER_LOOP >= 2
                        temp_sum[1] += temp_A * B_coluna_pointers[1][k] ;
                        #if OPERATIONS_PER_LOOP >= 4
                        temp_sum[2] += temp_A * B_coluna_pointers[2][k] ;
                        temp_sum[3] += temp_A * B_coluna_pointers[3][k] ;
                        #if OPERATIONS_PER_LOOP == 8
                        temp_sum[4] += temp_A * B_coluna_pointers[4][k] ;
                        temp_sum[5] += temp_A * B_coluna_pointers[5][k] ;
                        temp_sum[6] += temp_A * B_coluna_pointers[6][k] ;
                        temp_sum[7] += temp_A * B_coluna_pointers[7][k] ;
                        #endif
                        #endif
                        #endif
                    }
                    C_linha_pointer[coluna+0] = temp_sum[0];
                    #if OPERATIONS_PER_LOOP >= 2
                    C_linha_pointer[coluna+1] = temp_sum[1];
                    #if OPERATIONS_PER_LOOP >= 4
                    C_linha_pointer[coluna+2] = temp_sum[2];
                    C_linha_pointer[coluna+3] = temp_sum[3];
                    #if OPERATIONS_PER_LOOP == 8
                    C_linha_pointer[coluna+4] = temp_sum[4];
                    C_linha_pointer[coluna+5] = temp_sum[5];
                    C_linha_pointer[coluna+6] = temp_sum[6];
                    C_linha_pointer[coluna+7] = temp_sum[7];
                    #endif
                    #endif
                    #endif
                }
            }
        }
    }

    free( BT );
    return C;
}

int main( int argc , char** argv ){

    if( (OPERATIONS_PER_LOOP != 1 ) 
    && (OPERATIONS_PER_LOOP != 2 )
    && (OPERATIONS_PER_LOOP != 4 )
    &&( OPERATIONS_PER_LOOP != 8 ) ){
        printf("OPERATIONS_PER_LOOP deve ser 1 , 2 , 4 , 8 ou 16");
        return 2;
    }

    int matrix_tam = 1024 , block_tam , num_threads = 1;
    if( argc >= 2 ) matrix_tam = atoi( argv[1] );
    if( argc >= 3 ) block_tam = atoi( argv[2] ); else block_tam = matrix_tam;
    if( argc >= 4 ) num_threads = atoi( argv[3] ); else num_threads = 1;
    omp_set_num_threads(num_threads);
    if( block_tam < OPERATIONS_PER_LOOP ){
        printf("block_tam deve ser menor ou igual a OPERATIONS_PER_LOOP \n");
        return 2;
    }

    char nome_matriz_a[64] , nome_matriz_b[64] , nome_resultado[64];
    sprintf( nome_matriz_a , "A_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_matriz_b , "B_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_resultado , "resultado_%d_%d.txt" , matrix_tam , matrix_tam );

    FLOAT* A = le_matriz_vetor<FLOAT>( nome_matriz_a , matrix_tam );
    FLOAT* B = le_matriz_vetor<FLOAT>( nome_matriz_b , matrix_tam );

    FLOAT initial_t = omp_get_wtime();
    FLOAT* C = mult( matrix_tam , A , B , block_tam);
    FLOAT final_t = omp_get_wtime();

    printf("Duracao da multiplicacao(TAM=%d,block_tam=%d,OPERATIONS_PER_LOOP=%d,num_threads=%d): %.3fs\n" 
        , matrix_tam , block_tam , OPERATIONS_PER_LOOP , num_threads , final_t - initial_t );

    if( matrix_tam <= 1024 )
        escreve_matriz_sequencial( nome_resultado , matrix_tam , C );

}
