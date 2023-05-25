#include "fem_sne.h"

double Functions::u_real(double x, double t) {
   return x * t;
}

double Functions::u_g(double x, double t) {
   return x * t;
}

double Functions::f(double x, double t) {
   return x * t;
}

double Functions::theta(double x, double t) {
   return x * t;
}

double Functions::u_betta(double x, double t) {
   return x * t;
}

double Functions::lambda(double u) {
   return u;
}

Fem::Fem() {
   A = vector<vector<double>>();
   M = vector<vector<double>>();
   G = vector<vector<double>>();
   b = vector<double>();                   
   d = vector<double>();                   
   q = vector<double>();                   
   q_real = vector<double>();
   grid = vector<double>();
   times = vector<double>();
   q_init = vector<double>();
   loc_M = vector<vector<double>>();
   loc_G = vector<vector<double>>();
   loc_b = vector<double>();
   first_bc = vector<int>();
   second_bc = vector<int>();
   thrid_bc = vector<int>();

   nx = 0;                               
   n_times = 0;                          
   hx = 1, ht = 1, kx = 1, kt = 1;                
   betta = 0, sigma = 0;                  
   w = 1.6;                             
}

Fem::~Fem() {
   A.~vector();
   M.~vector();
   G.~vector();
   b.~vector();
   d.~vector();
   q.~vector();
   q_real.~vector();
   grid.~vector();
   times.~vector();
   q_init.~vector();
   loc_M.~vector();
   loc_G.~vector();
   loc_b.~vector();
   first_bc.~vector();
   second_bc.~vector();
   thrid_bc.~vector();
}

void Fem::read_data() {
   ifstream fgrid("grid.txt");

   fgrid >> betta >> sigma;

   fgrid >> nx >> hx >> kx;
   grid.resize(nx);
   fgrid >> grid[0];

   fgrid >> n_times >> ht >> kt;
   times.resize(n_times);
   fgrid >> times[0];

   fgrid.close();
}

void Fem::making_grid() {
   for (int i = 1; i < nx; i++)
      grid[i] = grid[0] + hx * pow(kx, i - 1);

   for (int i = 1; i < n_times; i++)
      times[i] = times[0] + ht * pow(kt, i - 1);

   loc_M.resize(2, vector<double>(2));
   loc_G.resize(2, vector<double>(2));

   loc_b.resize(2);

   d.resize(nx);
   q.resize(nx);
}

void Fem::init_cond() {
   Functions init;

   q_init.resize(nx);
   for (int i = 0; i < nx; i++)
      q_init[i] = init.u_real(grid[i], times[0]);
}

void Fem::glob_M() {
   double el1, el2;
   
   M.clear();
   M.resize(3, vector<double>(nx));

   for (int i = 0; i < nx - 1; i++) {
      el1 = (sigma * (grid[i + 1] - grid[i]) * 2) / 6.0;
      el2 = (sigma * (grid[i + 1] - grid[i])) / 6.0;

      loc_M[0][0] = loc_M[1][1] = el1;
      loc_M[0][1] = loc_M[1][0] = el2;

      M[0][i + 1] = loc_M[1][0];
      M[2][i] = loc_M[0][1];
      M[1][i] += loc_M[0][0];
      M[1][i + 1] += loc_M[1][1];
   }
}

void Fem::glob_G() {
   Functions L;
   double el1, el2;

   G.clear();
   G.resize(3, vector<double>(nx));

   for (int i = 0; i < nx - 1; i++) {
      el1 = (L.lambda(q[i]) + L.lambda(q[i + 1])) / (2.0 * (grid[i + 1] - grid[i]));
      el2 = - (L.lambda(q[i]) + L.lambda(q[i + 1])) / (2.0 * (grid[i + 1] - grid[i]));

      loc_G[0][0] = loc_G[1][1] = el1;
      loc_G[0][1] = loc_G[1][0] = el2;

      G[0][i + 1] = loc_G[1][0];
      G[2][i] = loc_G[0][1];
      G[1][i] += loc_G[0][0];
      G[1][i + 1] += loc_G[1][1];
   }
}

void Fem::glob_b(double t) { 
   Functions f;
   double el1, el2;

   b.clear();
   b.resize(nx);

   for (int i = 0; i < nx - 1; i++) {
      el1 = f.f(grid[i], t);
      el2 = f.f(grid[i + 1], t);

      loc_b[0] = (grid[i + 1] - grid[i]) * (2.0 * el1 + el2) / 6.0;
      loc_b[1] = (grid[i + 1] - grid[i]) * (el1 + 2.0 * el2) / 6.0;

      b[i] += loc_b[0];
      b[i + 1] += loc_b[1];
   }
}

void Fem::time_scheme(int idx) {
   double dt = times[idx] - times[idx - 1];
   vector<double> tmp(nx);
   vector<vector<double>> tmpM(3, vector<double>(nx));

   A.clear();
   A.resize(3, vector<double>(nx));

   glob_G();
   glob_M();

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < nx; j++) {
         A[i][j] = M[i][j] / dt + G[i][j];
      }
   }

   d.clear();
   d.resize(nx);

   glob_b(times[idx]);
   d = b;

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < nx; j++) {
         tmpM[i][j] = M[i][j] / dt;
      }
   }
   tmp = matr_vec_mult(q_init, tmpM);

   for (int i = 0; i < nx; i++) d[i] += tmp[i];

   boundary(times[idx]);
}

vector<double> Fem::matr_vec_mult(vector<double> &x, vector<vector<double>> &A) {
   vector<double> result(x.size());
   for (int k = 0; k < x.size(); k++) {
      for (int l = 0; l < 3; l++) {
         if (A[l][k] != 0)
            result[k] += A[l][k] * x[k + l - 1];
      }
   }
   return result;
}

void Fem::boundary(double t) {
   int type_left, type_right;
   Functions f;

   ifstream bc("boundary.txt");
   
   bc >> type_left >> type_right;
   
   switch (type_left)
   {
   case 1:
      A[1][0] = 1.0; A[2][0] = 0.0;
      d[0] = f.u_g(grid[0], t);
      break;
   case 2:
      d[0] += - f.theta(grid[0], t);
      break;
   case 3:
      A[1][0] += betta;
      d[0] += - betta * f.u_betta(grid[0], t);
      break;
   }
   switch (type_right)
   {
   case 1:
      A[1][nx - 1] = 1.0; A[0][nx - 1] = 0.0;
      d[nx - 1] = f.u_g(grid[nx - 1], t);
      break;
   case 2:
      d[nx - 1] += f.theta(grid[nx - 1], t);
      break;
   case 3:
      A[1][nx - 1] += betta;
      d[nx - 1] += betta * f.u_betta(grid[nx - 1], t);
      break;
   }

}

double Fem::residual() {
   vector<double> vec(nx);
   double res = 0.0; double nb = 0.0;
   for (int i = 0; i < nx; i++) {
      vec[i] = A[0][i] * q[i - 1] + A[1][i] * q[i] + A[2][i] * q[i + 1] - d[i];
      res += vec[i] * vec[i]; nb += d[i] * d[i];
   }
   return sqrt(res / nb);
}

void Fem::relax(vector<double> &q_prev) {
   for (int i = 0; i < nx; i++) q[i] = w * q[i] + (1.0 - w) * q_prev[i];
}


void Fem::errLine(ofstream &fout, int idx) {
   Functions real;
   fout << "Error:;;;";
   for (int i = 0; i < nx; i++) fout << abs(q[i] - real.u_real(grid[i], times[i])) << ';';
}

void Fem::iterLine(double resid, int num, ofstream &fout) {
   fout << num << ");" << resid << ";" << log(resid) << ";";
   for (int i = 0; i < nx; i++)
      fout << q[i] << ";";
   fout << endl;
}

void Fem::FPI() {
   ofstream fout("answer.txt");
   int num = 0;
   int maxiter = 1000;
   double eps = 1e-14;
   double resid = 10.0;

   vector<double> q_prev(nx);
   for (int i = 0; i < nx; i++) q[i] = 1.0;
   for (int i = 0; i < nx; i++) q_prev[i] = q[i];

   resid = residual();

   for (int i = 1; i < n_times; i++) {
      while (resid > eps && num < maxiter) {
         time_scheme(i);
         LU();
         relax(q_prev);
         num++;
         resid = residual();
         q_prev = q;
         iterLine(resid, num, fout);
      }
      errLine(fout, i);
   }
}

void Fem::LU() {
   for (int i = 1; i < nx; i++)
   {
      A[0][i] = A[0][i] / A[1][i - 1];
      A[1][i] = A[1][i] - A[0][i] * A[2][i - 1];
   }
   vector<double> y(nx);
   y[0] = d[0];
   for (int i = 1; i < nx; i++)
      y[i] = d[i] - A[0][i] * y[i - 1];
   q[nx - 1] = y[nx - 1] / A[1][nx - 1];
   for (int i = nx - 2; i >= 0; i--)
      q[i] = (y[i] - A[2][i] * q[i + 1]) / A[1][i];

   q_init = q;
}
