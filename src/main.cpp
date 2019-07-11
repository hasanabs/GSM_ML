#include <iostream>
#include <math.h>
#include<vector>
#include <complex>
#include "wits_lab.h"
using namespace std;
WITS gsm; GSM ant;
int Nt=16; int Np=2; int Nr=4; int M=4;
int SNR_min=0, SNR_max=20, step=2;
int Es=Np;
unsigned long transmit_repeat=500;

long combination = gsm.kombinasi(Nt,Np);
long L1=floor(log(combination)/log(2));
long L2= Np*log(M)/log(2);
long L=L1+L2;
long N=pow(2,L1);

int main() {
    gsm.tic();
    gsm.init_rand();

    /*---------TAC---------*/
    vector<vector<int>> prob_TAC = gsm.nchoosek(Nt,Np);
    vector<vector<int>> TAC=gsm.nfirst_nck(prob_TAC,N);
    /*---------Generate compare vector x---------*/
    vector<vector<int>> compare_bit= ant.generate_order_bin(L);
    vector<vector<complex<double>>> x_compare(Nt,vector<complex<double>>(pow(2,L),0)); // x=zeros(Nt,2^L) in complex
    for(int i=0;i<pow(2,L);i++){
        vector<vector<int>> temp01={compare_bit[i]};
        vector<vector<complex<double>>> temp02= ant.encoding(temp01,L1,M,Nt,Np,TAC);
        for (int j=0;j<Nt;j++){
            x_compare[j][i]=temp02[j][0];
        }
        temp02.clear();
    }

    /*---------SNR Loop---------*/
    for (int SNR=SNR_min;SNR<=SNR_max;SNR+=step){
        unsigned long err_acumulation=0;
        for (unsigned long transmit_time=0;transmit_time<transmit_repeat;transmit_time++){
            vector<vector<int>> bit = gsm.matrix_bitgenerator(1,L);
            vector<vector<complex<double>>> H = gsm.complex_randn(Nr,Nt,sqrt(0.5));
            vector<vector<complex<double>>> x  = ant.encoding(bit,L1,M,Nt,Np,TAC);
            vector<vector<complex<double>>> Hx = gsm.mul_matrix(H,x);
            vector<vector<complex<double>>> n = gsm.complex_randn(Nr,1,sqrt(Es/2/pow(10,(double)SNR/(double)10)));
            vector<vector<complex<double>>> y = gsm.sum_matrix(Hx,n);
            // /*---------ML Detection---------*/
            vector<vector<complex<double>>> x_temp(Nt,vector<complex<double>>(1,0)); // vector Nt x 1
            double min_dist=18; int min_dist_index=0;
            for (int i=0; i<pow(2,L);i++){
                for (int j=0;j<Nt;j++){
                    x_temp[j][0]=x_compare[j][i];
                }
                vector<vector<complex<double>>> Hx_compare = gsm.mul_matrix(H,x_temp);
                vector<vector<complex<double>>>y_Hx=gsm.sub_matrix(y,Hx_compare);
                double distance=gsm.frobenius_norm(y_Hx);
                if (distance<min_dist){
                    min_dist=distance;
                    min_dist_index=i;
                }
            }
            vector<vector<int>> bit_detected= {compare_bit[min_dist_index]};
            for (int i=0;i<L;i++){
                err_acumulation+=bit[0][i]^bit_detected[0][i];
            }
        }
        cout.precision(5);
        cout <<"SNR = " << SNR <<", number of Error = " <<scientific <<(double) err_acumulation/(double) (L*transmit_repeat) <<endl;
    }
    cout<<fixed<<endl;
    gsm.toc();
}