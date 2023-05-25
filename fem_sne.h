#pragma once

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Fem {
private:
   vector<vector<double>> A, M, G;        // матрицы в ленточном формате 
   vector<vector<double>> loc_M, loc_G;   // локальные матрицы массы и жесткости
   vector<double> b;                      // вектор правой части 
   vector<double> d;                      // вектор правой части нестационарной задачи
   vector<double> q;                      // вектор решения 
   vector<double> q_real;                 // вектор ожидаемых весов
   vector<double> grid;                   // сетка по x
   vector<double> times;                  // сетка по времени
   vector<double> q_init;                 // начальное условие
   vector<double> loc_b;                  // локальный вектор правой части 
   vector<int> first_bc, second_bc, thrid_bc; // узлы, в которых выполнены ку

   int nx;                                // кол-во конечных элементов, кол-во узлов
   int n_times;                           // количество временных слоев
   double hx, ht, kx, kt;                 // шаги сетки
   double betta, sigma;                   // коэффициенты
   double w;                              // параметр релаксации
public:
   Fem();
   ~Fem();

   void read_data();                      // ввод данных
  
   void making_grid();                    // составление стеки
   void init_cond();
   
   void glob_M();                         // вычисление глобальной матрицы масс
   void glob_G();                         // вычисление глобальной матрицы жесткости
   void glob_b(double t);                 // вычисление глобального вектора правой части
   
   void time_scheme(int idx);             // двуслойная неявная схема
   
   vector<double> matr_vec_mult(vector<double> &x, vector<vector<double>> &A); // умножение матрицы на вектор

   void boundary(double t);               // учет краевых условий
   void LU();                             // LU-разложение и решение системы ур-ий

   void errLine(ofstream &fout, int idx); // функции вывода 
   void iterLine(double resid, int num, ofstream &fout);

   double residual();                     // подсчет невзяки
   void relax(vector<double> &q_prev);    // применение параметра релаксации
   void FPI();                            // метод простой итерации
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
