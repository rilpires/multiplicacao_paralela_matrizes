#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

template <typename T>
T* le_matriz_vetor( const char* nome_arquivo , size_t tam_matrix ){
    T* ret = new T[tam_matrix*tam_matrix]();
    ifstream f(nome_arquivo);
    if( f.is_open() )
        for( int i = 0 ; i < tam_matrix*tam_matrix ; i++ )
            f >> ret[i];
    else
        for( int i = 0 ; i < tam_matrix*tam_matrix ; i++ )
            ret[i] = ((T)rand())/RAND_MAX ;
    
    f.close();
    return ret;
}

template <typename T>
T** le_matriz( const char* nome_arquivo , size_t tam_matrix ){
    T** ret = new T*[tam_matrix]();
    for( int i = 0 ; i < tam_matrix ; i++ )
        ret[i] = new T[tam_matrix]();
    ifstream f(nome_arquivo);
    
    if( f.is_open() )
        for( int i = 0 ; i < tam_matrix ; i++ )
            for( int j = 0 ; j < tam_matrix ; j++ )
                f >> ret[i][j];
    else
        for( int i = 0 ; i < tam_matrix ; i++ )
            for( int j = 0 ; j < tam_matrix ; j++ )
                ret[i][j] = ((T)rand())/RAND_MAX ;
    
    f.close();
    return ret;
}

template <typename T>
void escreve_matriz_sequencial( const char* nome_arquivo , size_t tam_matrix , T* matriz){
    ofstream f(nome_arquivo);
    for( int i = 0 ; i < tam_matrix ; i++ ){
        for( int j = 0 ; j < tam_matrix ; j++ )
            f << matriz[ i*tam_matrix + j ] << "\t";
        f << endl;
    }
    f.close();
}

template <typename T>
void escreve_matriz( const char* nome_arquivo , size_t tam_matrix , T** matriz){
    ofstream f(nome_arquivo);
    for( int i = 0 ; i < tam_matrix ; i++ ){
        for( int j = 0 ; j < tam_matrix ; j++ )
            f << matriz[i][j] << "\t";
        f << endl;
    }
    f.close();
}

template <typename T>
void escreve_matriz_aleatoria( const char* nome_arquivo , size_t tam_matrix ){
    ofstream f(nome_arquivo);
    for( int i = 0 ; i < tam_matrix ; i++ ){
        for( int j = 0 ; j < tam_matrix ; j++ )
            f << setprecision(2) << setw(4) << setfill('0') << ((T)rand())/RAND_MAX << "\t";
        f << endl;
    }
    f.close();
}



