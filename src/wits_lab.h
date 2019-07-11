// This function is to emulate frequently mathematical operation in that usualy utilize in MATLAB (R)
// Copyright (C) 2019 Wireless Information and Technology System LAB, NSYSU, Taiwan.
// Written by Hasan Albinsaid <hasanalbinsaid@hotmail.com>

#ifndef WITS_H
#define WITS_H
#include <vector>
#include <complex>
using namespace std;

class WITS{
    public:
        long permutasi(int n);
        long kombinasi (int x, int y);
        vector<vector<int>> matrix_bitgenerator(int row,int column);
        vector<vector<int>> nchoosek(int x,int y);
        void makecombi(vector<vector<int>> &matrik, vector<int> &vect, int x,int y, int left);
        vector<vector<int>> nfirst_nck(vector<vector<int>> &matrik, int n);
        void show_matrix(vector<vector<int>> &matrik);
        void init_rand();
        void tic();
        void toc();
        double randn();
        void show_complex(complex<double> &komplek, int precise);
        vector<vector< complex<double>>> complex_randn(int row, int column, double normalization);
        void show_complex_matrix(vector<vector< complex<double>>> &matrik, int precise);
        long int bin2dec(vector<vector<int>> &matrik,int firstbit, int lastbit); //order 1-n based on row first
        vector<vector< complex<double>>> mul_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB);
        vector<vector< complex<double>>> sum_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB);
        vector<vector< complex<double>>> sub_matrix(vector<vector< complex<double>>> &matrikA, vector<vector< complex<double>>> &matrikB);
        double frobenius_norm(vector<vector< complex<double>>> &matrik);
};
class GSM{
    public:
    complex<double> BPSK(int n);
    complex<double> QPSK(int n);
    complex<double> QAM16(int n);
    vector<vector<int>> generate_order_bin(int bit_lenght);
    vector<vector<complex<double>>> encoding(vector<vector<int>> &bin_data, int L1, int M, int Nt, int Np, vector<vector<int>> &TAC);
};

#endif