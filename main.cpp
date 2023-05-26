#include "fem_sne.h"

int main() {
   Fem matrix;
   matrix.read_data();
   matrix.making_grid();
   matrix.init_cond();
   matrix.FPI();
   //matrix.LU();
   return 0;
}
