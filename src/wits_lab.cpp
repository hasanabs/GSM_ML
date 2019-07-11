#include<vector>
#include <iostream>
#include <math.h>
#include <complex>
#include <time.h>
#include <stdlib.h>
#include <chrono>
#include<iomanip>
#include <limits>
#include <algorithm>
#include "wits_lab.h"
#define PI 3.1415926536
using namespace std;
chrono::time_point<chrono::system_clock> start, end_time;

long WITS::permutasi(int n){
    if (n<=1) {
        return 1;
    }
    else{
        return n*WITS::permutasi(n-1);
    }
}

long WITS::kombinasi(int x, int y){
    return WITS::permutasi(x)/WITS::permutasi(y)/WITS::permutasi(x-y);
}

vector<vector<int>> WITS::matrix_bitgenerator(int row,int column){
    vector<vector<int>> matrik;   matrik.clear();
    for (int x=0;x<row;x++){
        vector<int> temp;
        for (int y=0;y<column;y++){
            temp.push_back(rand()%2);
        }
        matrik.push_back(temp);
    }
    return matrik;
}

void WITS::makecombi(vector<vector<int>> &matrik, vector<int> &vect, int x,int y, int left){
    if (y==0){
        matrik.push_back(vect);
        return;
    }
    for (int i=left;i<=x;i++){
        vect.push_back(i);
        WITS::makecombi(matrik, vect, x, y-1, i+1);
        vect.pop_back();
    }
}

vector<vector<int>> WITS::nfirst_nck(vector<vector<int>> &matrik, int n){
    vector<vector<int>> matrik_temp; matrik_temp.clear();
    for (int i=0;i<n;i++){
        matrik_temp.push_back(matrik[i]);
    }
    return matrik_temp;
}

vector<vector<int>> WITS::nchoosek(int x,int y){
    vector<vector<int>> matrik; matrik.clear();
    vector<int> vect;
    makecombi(matrik, vect,x,y,1);
    return matrik;
}

void WITS::show_matrix(vector<vector<int>> &matrik){
    for (long unsigned int x=0;x<matrik.size();x++){
        for (long unsigned int y=0;y<matrik[x].size();y++){
            cout <<noshowpos << matrik[x][y]<< " ";
        }
        cout << endl;
    }
}

void WITS::init_rand(){
    srand(time(0));
}

void WITS::tic(){
    start = chrono::system_clock::now();
}

void WITS::toc(){
    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start; 
    std::time_t finish_at = std::chrono::system_clock::to_time_t(end_time);
    cout<<"Elapsed time: " <<setprecision(6) <<elapsed_seconds.count() <<" Second\n";
    cout<<"Finished computation at: " <<std::ctime(&finish_at)<<endl;
}

double WITS::randn()
{
    double temp1;
    double temp2;
    double result;
    int p=1;

    while( p > 0 )
    {
        temp2 = ( rand() / ( (double)RAND_MAX ) ); 
        if ( temp2 == 0 ){
        p = 1;
        }
        else{
        p = -1;
        }
    }
    temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
    result = sqrt( -2.0 * log( temp2 ) ) * temp1;
    return result;  
}

void WITS::show_complex(complex<double> &komplek, int precise){
    cout.setf( ios::fixed );
    cout <<setprecision(precise) <<noshowpos <<setw(precise+3) <<real(komplek);
    cout <<setprecision(precise) <<showpos <<setw(precise+2) <<imag(komplek) <<"i ";   
}

vector<vector< complex<double> >> WITS::complex_randn(int row, int column, double normalization){
    vector<vector< complex<double>>> matrik;
    vector< complex<double>> vect;
    matrik.clear();
    for (int i=0;i<row;i++){
        for (int j=1;j<=column;j++){
            WITS a01;
            complex<double> mycomplex(a01.randn()*normalization,a01.randn()*normalization);
            vect.push_back(mycomplex);
        }
        matrik.push_back(vect);
        vect.clear();
    }
    return matrik;
}

void WITS::show_complex_matrix(vector<vector< complex<double>>> &matrik, int precise){
    for (long unsigned int x=0;x<matrik.size();x++){
        for (long unsigned int y=0;y<matrik[x].size();y++){
            WITS a02;
            a02.show_complex(matrik[x][y], precise);
        }
        cout << endl;
    }
}

long int WITS::bin2dec(vector<vector<int>> &matrik,int firstbit, int lastbit){ //first bit start from 1
    int length=lastbit-firstbit+1;
    int matrix_size=matrik.size()*matrik[0].size();
    try{
        if (length>matrix_size){
            throw matrix_size;
        }
        else{
            int matrix_width=matrik[0].size();
            //initialization pointer position
            int Row = lastbit/(matrix_width+1);
            int Column=(lastbit-1)%matrix_width;
            int power=0;
            long int result=0;
            for(int i=length;i>0;i--){
                result=result+matrik[Row][Column]*pow(2,power);
                power++;
                if (Column%matrix_width==0){
                    Row--;
                    Column=matrix_width;
                }
                Column--;
            }
            return result; result=0;
        }
    }
    catch (int n){
        cout<<"ERROR: Out of range conversion bit, your lenght should less then " <<n <<endl;
        return 0;
    }
}

vector<vector< complex<double>>> WITS::mul_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB){
    int row=matrikA.size();
    int column=matrikB[0].size();
    int sum=matrikA[0].size();
    vector<vector<complex<double>>> result(row,vector<complex<double>>(column,0)); // zeros(row,column) in complex
    try{
        if (matrikA[0].size()==matrikB.size()){
            complex<double> temp;  complex <double> scalar_mul;
            for (int i=0;i<row;i++){
                for (int j=0;j<column;j++){
                    for (int k=0;k<sum;k++){
                        scalar_mul = matrikA[i][k]*matrikB[k][j];
                        temp=temp+scalar_mul;
                    }
                    result[i][j] = temp;
                    temp={0,0};
                }
            }
            return result;
        }
        else{
            throw 'e';
        }
    }
    catch(char x){
         cout << "ERROR: Matrix dimensions are not possible to do multiplication operations\n";
         return result; 
    }
}

vector<vector< complex<double>>> WITS::sum_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB){
    vector<vector<complex<double>>> result(matrikA.size(),vector<complex<double>>(matrikA[0].size(),0)); // zeros(row,column) in complex
    try{
        if(matrikA.size()==matrikB.size() && matrikA[0].size()==matrikB[0].size()){
            for (unsigned int i=0;i<matrikA.size();i++){
                for (unsigned int j=0;j<matrikA[0].size();j++){
                    result[i][j]=matrikA[i][j]+matrikB[i][j];
                }
            }
        return result;
        }
        else{
            throw 'e';
        }
    }
    catch(char x){
        cout << "ERROR: Matrix dimensions are not possible to do summation operations\n";
        return result; 
    }
    result.clear();
}

vector<vector< complex<double>>> WITS::sub_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB){
    vector<vector<complex<double>>> result(matrikA.size(),vector<complex<double>>(matrikA[0].size(),0)); // zeros(row,column) in complex
    try{
        if(matrikA.size()==matrikB.size() && matrikA[0].size()==matrikB[0].size()){
            for (unsigned int i=0;i<matrikA.size();i++){
                for (unsigned int j=0;j<matrikA[0].size();j++){
                    result[i][j]=matrikA[i][j]-matrikB[i][j];
                }
            }
        return result;
        }
        else{
            throw 'e';
        }
    }
    catch(char x){
        cout << "ERROR: Matrix dimensions are not possible to do subtraction operations\n";
        return result;
    }
    result.clear();
}

double WITS::frobenius_norm(vector<vector< complex<double>>> &matrik){
    double result=0;
    for (unsigned int i=0;i<matrik.size();i++){
        for (unsigned int j=0;j<matrik[0].size();j++){
            result=result+norm(matrik[i][j]);//pow(real(matrik[i][j]),2)+pow(imag(matrik[i][j]),2);
        }
    }
    return sqrt(result);
}

complex<double> GSM::BPSK(int n){
    vector<complex<double>> conselation={complex<double>(-1,0), complex<double>(1,0)};
    return conselation[n];
}

complex<double> GSM::QPSK(int n){
    vector<complex<double>> conselation={complex<double>(-1,-1)/sqrt(2), complex<double>(-1,1)/sqrt(2), complex<double>(1,-1)/sqrt(2), complex<double>(1,1)/sqrt(2)};
    return conselation[n];
}

complex<double> GSM::QAM16(int n){
    vector<complex<double>> conselation=
    {
        complex<double>(-3,-3)/sqrt(10), complex<double>(-3,-1)/sqrt(10), complex<double>(-3,3)/sqrt(10),  complex<double>(-3,1)/sqrt(10),
        complex<double>(-1,-3)/sqrt(10), complex<double>(-1,-1)/sqrt(10), complex<double>(-1,3)/sqrt(10),  complex<double>(-1,1)/sqrt(10),
        complex<double>(3,-3)/sqrt(10), complex<double>(3,-1)/sqrt(10), complex<double>(3,3)/sqrt(10),  complex<double>(3,1)/sqrt(10),
        complex<double>(1,-3)/sqrt(10), complex<double>(1,-1)/sqrt(10), complex<double>(1,3)/sqrt(10),  complex<double>(1,1)/sqrt(10)
    };
    return conselation[n];
}

vector<vector<int>> GSM::generate_order_bin(int bit_lenght){
    int N[bit_lenght];
    fill_n(N,bit_lenght,1); //Inizialitation 1
    int number_bit=pow(2,bit_lenght);
    vector<vector<int>> result(number_bit,vector<int>(bit_lenght,0)); // x=zeros(2^bit,bit_lenght) in integer
    for (int row=0;row<pow(2,bit_lenght);row++){
        for (int column=0;column<bit_lenght;column++){
            int x=pow(2,(bit_lenght-1-column));
            if (row%x==0){
                N[column]^=1;
            }
            result[row][column]=N[column];
        }
    }
    return result;
}

vector<vector<complex<double>>> GSM::encoding(vector<vector<int>> &bin_data, int L1, int M, int Nt, int Np, vector<vector<int>> &TAC){
    vector<vector<complex<double>>> x_vect(Nt,vector<complex<double>>(1,0)); // x=zeros(Nt,1) in complex
    try{
        if (M==2 || M==4 || M==16){
            int bit_each_symbol=log(M)/log(2);
            WITS a03;
            int selected_TAC=a03.bin2dec(bin_data,1, L1);
            for (int i=0;i<Np;i++){
                if (M==2){
                    x_vect[TAC[selected_TAC][i]-1][0]=GSM::BPSK(a03.bin2dec(bin_data,L1+1+bit_each_symbol*i, L1+bit_each_symbol+bit_each_symbol*i));
                }
                else if (M==4){
                    x_vect[TAC[selected_TAC][i]-1][0]=GSM::QPSK(a03.bin2dec(bin_data,L1+1+bit_each_symbol*i, L1+bit_each_symbol+bit_each_symbol*i));
                }
                else if(M==16){
                     x_vect[TAC[selected_TAC][i]-1][0]=GSM::QAM16(a03.bin2dec(bin_data,L1+1+bit_each_symbol*i, L1+bit_each_symbol+bit_each_symbol*i));
                }
            }
        }
        else{
            throw M;
        }
    }
    catch (int n){
        cout << "ERROR: Your M conselation is = " <<n <<" its not suitable for this library\n";
    }
    return x_vect;
}