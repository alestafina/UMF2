#pragma once

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Fem {
private:
   vector<vector<double>> A, M, G;        // ������� � ��������� ������� 
   vector<vector<double>> loc_M, loc_G;   // ��������� ������� ����� � ���������
   vector<double> b;                      // ������ ������ ����� 
   vector<double> d;                      // ������ ������ ����� �������������� ������
   vector<double> q;                      // ������ ������� 
   vector<double> q_real;                 // ������ ��������� �����
   vector<double> grid;                   // ����� �� x
   vector<double> times;                  // ����� �� �������
   vector<double> q_init;                 // ��������� �������
   vector<double> loc_b;                  // ��������� ������ ������ ����� 
   vector<int> first_bc, second_bc, thrid_bc; // ����, � ������� ��������� ��

   int nx;                                // ���-�� �������� ���������, ���-�� �����
   int n_times;                           // ���������� ��������� �����
   double hx, ht, kx, kt;                 // ���� �����
   double betta, sigma;                   // ������������
   double w;                              // �������� ����������
public:
   Fem();
   ~Fem();

   void read_data();                      // ���� ������
  
   void making_grid();                    // ����������� �����
   void init_cond();
   
   void glob_M();                         // ���������� ���������� ������� ����
   void glob_G();                         // ���������� ���������� ������� ���������
   void glob_b(double t);                 // ���������� ����������� ������� ������ �����
   
   void time_scheme(int idx);             // ���������� ������� �����
   
   vector<double> matr_vec_mult(vector<double> &x, vector<vector<double>> &A); // ��������� ������� �� ������

   void boundary(double t);               // ���� ������� �������
   void LU();                             // LU-���������� � ������� ������� ��-��

   void errLine(ofstream &fout, int idx); // ������� ������ 
   void iterLine(double resid, int num, ofstream &fout);

   double residual();                     // ������� �������
   void relax(vector<double> &q_prev);    // ���������� ��������� ����������
   void FPI();                            // ����� ������� ��������
};

class Functions {
public:
   double u_real(double x, double t);
   double u_g(double x, double t);
   double f(double x, double t);
   double theta(double x, double t);
   double u_betta(double x, double t);
   double lambda(double u);
};
