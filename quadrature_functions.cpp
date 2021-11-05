#include "quadrature_header.h"

double QuadratureRuleTrapezoidal(double (*F)(double), double a, double b, int N){
  
  /*
  Estimate the value of an integral over the interval [a, b] using the Midpoint method
  
  F     function to be integrated
  a     left side of the interval
  b     right side of the interval
  N     number of sub-intervals
  */
  
  double h, t, sum=0.;
  int i;
  
  h = (b - a)/N;
  
  for(i = 1; i < N; i++){
    t = a + i*h;
    sum += h*F(t);
  }
  
  sum += 0.5*h*F(a) + 0.5*h*F(b);
  
  return sum;
}

double QuadratureRuleSimpson(double (*F)(double), double a, double b, int N){
  
  /*
  Estimate the value of an integral over the interval [a, b] using the Simpson method
  
  F     function to be integrated
  a     left side of the interval
  b     right side of the interval
  N     number of sub-intervals
  */
  
  double h, t, sum=0., hred;
  int i;
  
  h = (b - a)/N;
  hred = h/3;
  
  for(i = 1; i < N; i++){
    t = a + i*h;
    sum += (2 + 2*(i%2))*F(t);
  }
 
  sum += F(a) + F(b);
  
  return sum*hred;
}

double QuadratureRuleGauss(double (*F)(double), double a, double b, int N, int Ng){
  double w[Ng], x[Ng];
  double sum, sumk;
  double x0, x1;
  int i, k;
  
  w[0] = 5./9.; w[1] = 8./9.; w[2] = 5./9.;
  x[0] = -sqrt(3./5.); x[1] = 0; x[2] = sqrt(3./5.);
  
  sum = 0.;
  
  for (i = 0; i < N; i++){
    x0 = (b - a)*.5;
    x1 = (a + b)*.5;
    sumk = 0.0;
    for (k = 0; k < Ng; k++){
      sumk += w[k]*F(x0*x[k] + x1);
    }
    sum += sumk;
  }
  
  sum *= x0;
  
  return sum;
}

double F(double t){
  if(t != 0) return sin(t)/t;
  return 1;
}