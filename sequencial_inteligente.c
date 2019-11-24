#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef OPERATIONS_PER_LOOP
    #define OPERATIONS_PER_LOOP 4
#endif

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
    // printf("Duracao do gera_transposta: %.3fs\n" , final_t-initial_t);
    return ret;
}

double* mult( size_t matrix_tam , double* A , double* B , size_t block_tam ){
    double* C = malloc( sizeof(double) * matrix_tam * matrix_tam );
    int block_x,block_y,linha,coluna,k, num_blocos = matrix_tam/block_tam , lim1 , lim2 ;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp_A;
    double *C_linha_pointer , *A_linha_pointer;
    double *B_coluna_pointer1 , *B_coluna_pointer2 , *B_coluna_pointer3 , *B_coluna_pointer4; 
    double *B_coluna_pointer5 , *B_coluna_pointer6 , *B_coluna_pointer7 , *B_coluna_pointer8; 
    double* BT = gera_transposta(matrix_tam,B);

    for( block_x = 0 ; block_x < num_blocos ; block_x++ )
    for( block_y = 0 ; block_y < num_blocos ; block_y++ )
    {

        lim1 = (block_y+1)*block_tam;
        for( linha = block_y*block_tam ; linha < lim1 ; linha++ ){
            C_linha_pointer = C + linha*matrix_tam;
            A_linha_pointer = A + linha*matrix_tam;
            lim2 = (block_x+1)*block_tam;
            for( coluna = block_x*block_tam ; coluna < lim2 ; coluna += OPERATIONS_PER_LOOP ){
                temp1 = 0;
                B_coluna_pointer1 = BT + (coluna+0)*matrix_tam;
                #if OPERATIONS_PER_LOOP >= 2
                temp2 = 0;
                B_coluna_pointer2 = BT + (coluna+1)*matrix_tam;
                #if OPERATIONS_PER_LOOP >= 4
                temp3 = 0;
                temp4 = 0;
                B_coluna_pointer3 = BT + (coluna+2)*matrix_tam;
                B_coluna_pointer4 = BT + (coluna+3)*matrix_tam;
                #if OPERATIONS_PER_LOOP == 8
                temp5 = 0;
                temp6 = 0;
                temp7 = 0;
                temp8 = 0;
                B_coluna_pointer5 = BT + (coluna+4)*matrix_tam;
                B_coluna_pointer6 = BT + (coluna+5)*matrix_tam;
                B_coluna_pointer7 = BT + (coluna+6)*matrix_tam;
                B_coluna_pointer8 = BT + (coluna+7)*matrix_tam;
                #endif
                #endif
                #endif
                for( k = 0 ; k < matrix_tam ; k ++ ){
                    temp_A = A_linha_pointer[k];
                    temp1 += temp_A * B_coluna_pointer1[k] ;
                    #if OPERATIONS_PER_LOOP >= 2
                    temp2 += temp_A * B_coluna_pointer2[k] ;
                    #if OPERATIONS_PER_LOOP >= 4
                    temp3 += temp_A * B_coluna_pointer3[k] ;
                    temp4 += temp_A * B_coluna_pointer4[k] ;
                    #if OPERATIONS_PER_LOOP == 8
                    temp5 += temp_A * B_coluna_pointer5[k] ;
                    temp6 += temp_A * B_coluna_pointer6[k] ;
                    temp7 += temp_A * B_coluna_pointer7[k] ;
                    temp8 += temp_A * B_coluna_pointer8[k] ;
                    #endif
                    #endif
                    #endif
                }
                C_linha_pointer[coluna+0] = temp1;
                #if OPERATIONS_PER_LOOP >= 2
                C_linha_pointer[coluna+1] = temp2;
                #if OPERATIONS_PER_LOOP >= 4
                C_linha_pointer[coluna+2] = temp3;
                C_linha_pointer[coluna+3] = temp4;
                #if OPERATIONS_PER_LOOP == 8
                C_linha_pointer[coluna+4] = temp5;
                C_linha_pointer[coluna+5] = temp6;
                C_linha_pointer[coluna+6] = temp7;
                C_linha_pointer[coluna+7] = temp8;
                #endif
                #endif
                #endif
            }
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

    #if (OPERATIONS_PER_LOOP != 1 )&&(OPERATIONS_PER_LOOP != 2 )&&(OPERATIONS_PER_LOOP != 4 )&&(OPERATIONS_PER_LOOP != 8 )
        printf("OPERATIONS_PER_LOOP deve ser 1 , 2 , 4 ou 8");
        return 2;
    #endif

    int matrix_tam = 1024;
    int block_tam;
    if( argc >= 2 )
        matrix_tam = atoi( argv[1] );
    if( argc >= 3 )
        block_tam = atoi( argv[2] );
    else
        block_tam = matrix_tam;

    double* A = malloc( sizeof(double)*matrix_tam*matrix_tam );
    double* B = malloc( sizeof(double)*matrix_tam*matrix_tam );
    int i;
    for( i = 0 ; i < matrix_tam*matrix_tam ; i++ ){
        A[i] = ((double)rand()) / RAND_MAX;
        B[i] = ((double)rand()) / RAND_MAX;
        // A[i] = i % matrix_tam; // Facil de checar corretude
        // B[i] = i / matrix_tam; // Facil de checar corretude
    }


    double initial_t = omp_get_wtime();
    double* C = mult( matrix_tam , A , B , block_tam);
    double final_t = omp_get_wtime();

    printf("Duracao da multiplicacao(TAM=%d,block_tam=%d): %.3fs\n" , matrix_tam , block_tam , final_t - initial_t );
    
    if( matrix_tam <= 16 ){
        printf("A:\n");
        printa_matriz(matrix_tam , A);
        printf("B:\n");
        printa_matriz(matrix_tam , B);
        printf("C:\n");
        printa_matriz(matrix_tam , C);
    }

}
