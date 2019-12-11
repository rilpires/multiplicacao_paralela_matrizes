#include "matriz_io.cpp"
#include <omp.h>
#include <openacc.h>

typedef double FLOAT;

void gera_transposta( size_t matrix_tam , FLOAT* MAT , FLOAT* MAT_T ){
    int i , j;
    FLOAT initial_t = omp_get_wtime();
    #pragma acc data    deviceptr(MAT,MAT_T)
    #pragma acc parallel loop tile(32,32)
    {
        for( i = 0 ; i < matrix_tam ; i++ )
        for( j = 0 ; j < matrix_tam ; j++ )
            MAT_T[i+j*matrix_tam] = MAT[j+i*matrix_tam];
    }
    FLOAT final_t = omp_get_wtime();
    //printf("Duracao do gera transposta: %fs\n" , final_t - initial_t );
}

void mult( size_t matrix_tam , FLOAT* MAT1 , FLOAT* MAT2 , FLOAT* PROD ){

    int     linha,coluna,i,j;
    FLOAT   temp_sum; 
    FLOAT*  MAT2_T = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam);
    
    // Funciona, mas demora muito mais. 
    // Fazemos a transposta na regiao data a seguir
    //gera_transposta( matrix_tam , MAT2 , MAT2_T );

    #pragma acc data \
        deviceptr(MAT1 , MAT2 , MAT2_T , PROD)
    {
        #pragma acc parallel loop tile(32,32)
        for( i = 0 ; i < matrix_tam ; i++ )
        for( j = 0 ; j < matrix_tam ; j++ )
            MAT2_T[i+j*matrix_tam] = MAT2[j+i*matrix_tam];
        
        #pragma acc parallel loop independent tile(8,8)
        for( linha = 0 ; linha < matrix_tam ; linha++ )
        for( coluna = 0 ; coluna < matrix_tam ; coluna ++ )
        {
            for( temp_sum = 0 , i = 0 ; i < matrix_tam ; i ++ )
                temp_sum += MAT1[i + linha*matrix_tam] * MAT2_T[i + coluna*matrix_tam];
            PROD[coluna + linha*matrix_tam] = temp_sum;
        }
    }
    
}

FLOAT* strassen( size_t matrix_tam , FLOAT* M1 , FLOAT* M2 ){

    FLOAT* Z = (FLOAT*)malloc(sizeof(FLOAT)*matrix_tam*matrix_tam);
    FLOAT *A =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* B =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* C =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* D =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* E =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* F =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* G =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* H =  (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P0 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P1 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P2 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P3 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P4 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P5 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* P6 = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* FH = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AB = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* CD = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* GE = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AD = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* EH = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* BD = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* GH = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* AC = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);
    FLOAT* EF = (FLOAT*)acc_malloc(sizeof(FLOAT)*matrix_tam*matrix_tam/4);  

    #pragma acc data copyin( M1[:matrix_tam*matrix_tam] , M2[:matrix_tam*matrix_tam] ) \
                     deviceptr( A , B , C , D , \
                                E , F , G , H ,  \
                                FH , AB , CD , GE , AD , \
                                EH , BD , GH , AC , EF , \
                                P0 , P1 , P2 , P3 ,  P4 , P5 , P6 ) \
                     copyout( Z[:matrix_tam*matrix_tam] )
    {
        int linha,coluna;

        #pragma acc parallel loop
        for(linha=0;linha<matrix_tam/2;linha++){
            int temp = linha*matrix_tam , temp2 = linha*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=0;coluna<matrix_tam/2;coluna++){
                temp3 = coluna+temp;
                temp4 = coluna+temp2;
                A[temp4] = M1[temp3];
                E[temp4] = M2[temp3];
            }
        }
        #pragma acc parallel loop
        for(linha=0;linha<matrix_tam/2;linha++){
            int temp = linha*matrix_tam , temp2 = linha*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
                temp3 = coluna+temp;
                temp4 = coluna+temp2-matrix_tam/2;
                B[temp4] = M1[temp3];
                F[temp4] = M2[temp3];
            }
        }
        #pragma acc parallel loop
        for(linha=matrix_tam/2;linha<matrix_tam;linha++){
            int temp = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=0;coluna<matrix_tam/2;coluna++){
                temp3 = coluna+temp;
                temp4 = coluna+temp2;
                C[temp4] = M1[temp3];
                G[temp4] = M2[temp3];
            }
        }
        #pragma acc parallel loop
        for(linha=matrix_tam/2;linha<matrix_tam;linha++){
            int temp = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
                temp3 = coluna+temp;
                temp4 = coluna+temp2-matrix_tam/2;
                D[temp4] = M1[temp3];
                H[temp4] = M2[temp3];
            }
        }


        size_t temp_tam2 = matrix_tam/2;
        #pragma acc parallel loop
        for( linha=0 ; linha<temp_tam2 ; linha++ ){
            size_t temp=linha*matrix_tam/2 , temp2;
            for( coluna=0 ; coluna<temp_tam2 ; coluna++ ){
                temp2 = coluna+temp;
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

        mult( matrix_tam/2 , A , FH, P0 );
        mult( matrix_tam/2 , AB , H, P1 );
        mult( matrix_tam/2 , CD , E, P2 );
        mult( matrix_tam/2 , D , GE, P3 );
        mult( matrix_tam/2 , AD , EH, P4 );
        mult( matrix_tam/2 , BD , GH, P5 );
        mult( matrix_tam/2 , AC , EF, P6 );

        #pragma acc parallel loop tile(8,8)
        for(linha=0;linha<matrix_tam/2;linha++){
            int temp = linha*matrix_tam , temp2 = linha*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=0;coluna<matrix_tam/2;coluna++){
                temp3 = temp+coluna;
                temp4 = temp2+coluna;
                Z[temp3] = P3[temp4] + P4[temp4] + P5[temp4] - P1[temp4];
            }
        }
        #pragma acc parallel loop tile(8,8)
        for(linha=0;linha<matrix_tam/2;linha++){
            int temp = linha*matrix_tam , temp2 = linha*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
                temp3 = temp + coluna;
                temp4 = temp2 + coluna - matrix_tam/2;
                Z[temp3] = P0[temp4] + P1[temp4];
            }
        }
        #pragma acc parallel loop tile(8,8)
        for(linha=matrix_tam/2;linha<matrix_tam;linha++){
            int temp = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=0;coluna<matrix_tam/2;coluna++){
                temp3 = temp + coluna;
                temp4 = temp2 + coluna;
                Z[temp3] = P2[temp4] + P3[temp4];
            }
        }
        #pragma acc parallel loop tile(8,8)
        for(linha=matrix_tam/2;linha<matrix_tam;linha++){
            int temp = linha*matrix_tam , temp2 = (linha-matrix_tam/2)*matrix_tam/2 , temp3 , temp4 ;
            for(coluna=matrix_tam/2;coluna<matrix_tam;coluna++){
                temp3 = temp + coluna;
                temp4 = temp2 + coluna - matrix_tam/2;
                Z[temp3] = P0[temp4] + P4[temp4] - P2[temp4] - P6[temp4];
            }
        }

    }

    acc_free(A ); 
    acc_free(B );
    acc_free(C );
    acc_free(D );
    acc_free(E );
    acc_free(F );
    acc_free(G );
    acc_free(H );
    acc_free(P0);
    acc_free(P1);
    acc_free(P2);
    acc_free(P3);
    acc_free(P4);
    acc_free(P5);
    acc_free(P6);
    acc_free(FH);
    acc_free(AB);
    acc_free(CD);
    acc_free(GE);
    acc_free(AD);
    acc_free(EH);
    acc_free(BD);
    acc_free(GH);
    acc_free(AC);
    acc_free(EF);

    return Z;
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
    FLOAT* C = strassen( matrix_tam , A , B );
    FLOAT final_t = omp_get_wtime();

    printf("[ACC STRASSEN] - Duracao da multiplicacao (TAM=%d): %.3fs\n" , matrix_tam , final_t - initial_t );

    if( matrix_tam <= 1024 )
        escreve_matriz_sequencial( nome_resultado , matrix_tam , C );

    printf("[Corretude]\nOs 3 primeiros valores de C:\n" );
    printf("C: %f %f %f \n" , C[0] , C[1] , C[2] );
    
}
