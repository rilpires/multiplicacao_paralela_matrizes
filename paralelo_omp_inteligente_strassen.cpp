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

    #pragma omp parallel for private(i,j) shared(ret,A,matrix_tam)
    for( j = 0 ; j < matrix_tam ; j++ ){
        for( i = 0 ; i < matrix_tam ; i++ ){
            ret[i+j*matrix_tam] = A[j+i*matrix_tam];
        }
    }

    FLOAT final_t = omp_get_wtime();
    // printf("Duracao do gera_transposta: %.3fs\n" , final_t-initial_t);
    return ret;
}

FLOAT* mult( size_t matrix_tam , FLOAT* A , FLOAT* B , size_t block_tam=8 ){
    FLOAT* C = (FLOAT*)malloc( sizeof(FLOAT) * matrix_tam * matrix_tam );
    int block_x,block_y,linha,coluna,k, num_blocos = matrix_tam/block_tam , lim1 , lim2 ;
    FLOAT temp_sum[8],temp_A;
    FLOAT *C_linha_pointer , *A_linha_pointer;
    FLOAT* B_coluna_pointers[8];
    FLOAT* BT = gera_transposta(matrix_tam,B);

    #pragma omp parallel
    #pragma omp master
    {
        for( block_x = 0 ; block_x < num_blocos ; block_x++ )
        for( block_y = 0 ; block_y < num_blocos ; block_y++ )
        {
            lim1 = (block_y+1)*block_tam;
            for( linha = block_y*block_tam ; linha < lim1 ; linha++ ){
                C_linha_pointer = C + linha*matrix_tam;
                A_linha_pointer = A + linha*matrix_tam;
                lim2 = (block_x+1)*block_tam;
                #pragma omp task private(k,B_coluna_pointers,temp_sum,coluna,temp_A) firstprivate(A_linha_pointer,C_linha_pointer,block_x,lim2)
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


FLOAT* strassen( size_t matrix_tam , FLOAT* M1 , FLOAT* M2 ){
    
    FLOAT* Z = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam);

    FLOAT* A =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* B =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* C =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* D =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* E =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* F =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* G =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* H =  (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P0 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P1 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P2 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P3 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P4 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P5 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P6 = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* FH = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AB = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* CD = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* GE = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AD = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* EH = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* BD = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* GH = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AC = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* EF = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);    
    
    size_t linha,coluna;
    int temp1 , temp2 , temp3 , temp4;
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=0;linha<matrix_tam/2;linha++){
        temp1 = linha*matrix_tam , temp2 = linha*matrix_tam/2;
        for(coluna=0;coluna<matrix_tam/2;coluna++){
            temp3 = coluna+temp1;
            temp4 = coluna+temp2;
            A[temp4] = M1[temp3];
            E[temp4] = M2[temp3];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=0;linha<matrix_tam/2;linha++){
        temp1 = linha*matrix_tam , temp2 = linha*matrix_tam/2;
        for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
            temp3 = coluna+temp1;
            temp4 = coluna+temp2-matrix_tam/2;
            B[temp4] = M1[temp3];
            F[temp4] = M2[temp3];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=matrix_tam/2;linha<matrix_tam;linha++){
        temp1 = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2;
        for(coluna=0;coluna<matrix_tam/2;coluna++){
            temp3 = coluna+temp1;
            temp4 = coluna+temp2;
            C[temp4] = M1[temp3];
            G[temp4] = M2[temp3];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=matrix_tam/2;linha<matrix_tam;linha++){
        temp1 = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2;
        for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
            temp3 = coluna+temp1;
            temp4 = coluna+temp2-matrix_tam/2;
            D[temp4] = M1[temp3];
            H[temp4] = M2[temp3];
        }
    }

    size_t temp_tam2 = matrix_tam/2;
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2) shared(matrix_tam)
    for( linha=0 ; linha<temp_tam2 ; linha++ ){
        temp1=linha*matrix_tam/2 , temp2;
        for( coluna=0 ; coluna<temp_tam2 ; coluna++ ){
            temp2 = coluna+temp1;
            FH[temp2] = F[temp2] - H[temp2];
            AB[temp2] = A[temp2] + B[temp2];
            CD[temp2] = C[temp2] + D[temp2];
            GE[temp2] = G[temp2] - E[temp2];
            AD[temp2] = A[temp2] + D[temp2];
            EH[temp2] = E[temp2] + H[temp2];
            BD[temp2] = B[temp2] - D[temp2];
            GH[temp2] = G[temp2] + H[temp2];
            AC[temp2] = A[temp2] - C[temp2];
            EF[temp2] = E[temp2] + F[temp2];
        }
    }
    // Ate que dimensao continuar a recursao com o Strassen
    if(matrix_tam <= 1024 ){
        P0 = mult( matrix_tam/2 , A , FH  , 8 );
        P1 = mult( matrix_tam/2 , AB , H  , 8 );
        P2 = mult( matrix_tam/2 , CD , E  , 8 );
        P3 = mult( matrix_tam/2 , D , GE  , 8 );
        P4 = mult( matrix_tam/2 , AD , EH , 8 );
        P5 = mult( matrix_tam/2 , BD , GH , 8 );
        P6 = mult( matrix_tam/2 , AC , EF , 8 );
    } else {
        P0 = strassen( matrix_tam/2 , A , FH  );
        P1 = strassen( matrix_tam/2 , AB , H  );
        P2 = strassen( matrix_tam/2 , CD , E  );
        P3 = strassen( matrix_tam/2 , D , GE  );
        P4 = strassen( matrix_tam/2 , AD , EH );
        P5 = strassen( matrix_tam/2 , BD , GH );
        P6 = strassen( matrix_tam/2 , AC , EF );
    }
    
    free(A);free(C);free(D);free(B);
    free(E);free(F);free(G);free(H);
    free(FH);free(AB);free(CD);free(GE);free(AD);
    free(EH);free(BD);free(GH);free(AC);free(EF);
    
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=0;linha<matrix_tam/2;linha++){
        temp1 = linha*matrix_tam , temp2 = linha*matrix_tam/2;
        for(coluna=0;coluna<matrix_tam/2;coluna++){
            temp3 = temp1+coluna;
            temp4 = temp2+coluna;
            Z[temp3] = P3[temp4] + P4[temp4] + P5[temp4] - P1[temp4];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=0;linha<matrix_tam/2;linha++){
        temp1 = linha*matrix_tam , temp2 = linha*matrix_tam/2;
        for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
            temp3 = temp1 + coluna;
            temp4 = temp2 + coluna - matrix_tam/2;
            Z[temp3] = P0[temp4] + P1[temp4];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=matrix_tam/2;linha<matrix_tam;linha++){
        temp1 = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2;
        for(coluna=0;coluna<matrix_tam/2;coluna++){
            temp3 = temp1 + coluna;
            temp4 = temp2 + coluna;
            Z[temp3] = P2[temp4] + P3[temp4];
        }
    }
    #pragma omp parallel for private(linha,coluna) private(temp1,temp2,temp3,temp4) shared(matrix_tam)
    for(linha=matrix_tam/2;linha<matrix_tam;linha++){
        temp1 = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2;
        for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
            temp3 = temp1 + coluna;
            temp4 = temp2 + coluna - matrix_tam/2;
            Z[temp3] = P0[temp4] + P4[temp4] - P2[temp4] - P6[temp4];
        }
    }

    free(P0);free(P1);free(P2);free(P3);
    free(P4);free(P5);free(P6);

    return Z;
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
    if( block_tam < OPERATIONS_PER_LOOP ){
        printf("block_tam deve ser menor ou igual a OPERATIONS_PER_LOOP \n");
        return 2;
    }
    omp_set_num_threads(num_threads);

    char nome_matriz_a[64] , nome_matriz_b[64] , nome_resultado[64];
    sprintf( nome_matriz_a , "A_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_matriz_b , "B_%d_%d.txt" , matrix_tam , matrix_tam );
    sprintf( nome_resultado , "resultado_%d_%d.txt" , matrix_tam , matrix_tam );

    FLOAT* A = le_matriz_vetor<FLOAT>( nome_matriz_a , matrix_tam );
    FLOAT* B = le_matriz_vetor<FLOAT>( nome_matriz_b , matrix_tam );

    FLOAT initial_t = omp_get_wtime();
    FLOAT* C = strassen( matrix_tam , A , B );
    FLOAT final_t = omp_get_wtime();

    printf("[OMP STRASSEN](TAM=%d,block_tam=%d,OPERATIONS_PER_LOOP=%d,num_threads=%d): %.3fs\n" 
        , matrix_tam , block_tam , OPERATIONS_PER_LOOP , num_threads , final_t - initial_t );

    if( matrix_tam <= 1024 )
        escreve_matriz_sequencial( nome_resultado , matrix_tam , C );

    printf("[Corretude]\nOs 3 primeiros valores de C:\n" );
    printf("C: %f %f %f \n" , C[0] , C[1] , C[2] );

}
