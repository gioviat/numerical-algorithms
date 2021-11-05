#include "quadrature_header.h"

int main(){
  double a, x;
  double sum_trapezoidal = 0., sum_simp = 0., sum_gauss = 0.;
  double exact = 1.618194443708, err_gauss, err_simp, err_trap;
  int N = 3, m = 1;
  
  ofstream fdata;
  fdata.open("Si.dat");
  
  cout << scientific << setprecision(4);
  
  a = 0.;
  
  cout << "x1          " << "Gauss       " << "Simpson     " << "Trapez." << endl;
  cout << "----------------------------------------------" << endl;
  
  for(x = 1; x <= 15; x++){
    sum_gauss += QuadratureRuleGauss(F, a, x, m, N);
    sum_simp += QuadratureRuleSimpson(F, a, x, N);
    sum_trapezoidal += QuadratureRuleTrapezoidal(F, a, x, N);
    cout << x << "  " << sum_gauss << "  " << sum_simp << "  " << sum_trapezoidal << endl;
    fdata << x << " " << sum_gauss << endl;
    a = x;
  }
  
  err_gauss = fabs(sum_gauss - exact);
  err_simp = fabs(sum_simp - exact);
  err_trap = fabs(sum_trapezoidal - exact);
  
  cout << "----------------------------------------------" << endl;
  cout << "Error     " << "  " << err_gauss << "  " << err_simp << "  " << err_trap << endl;
  fdata.close();
  
  return 0;
}
