#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double QuadratureRuleSimpson(double (*F)(double), double, double, int);
double QuadratureRuleTrapezoidal(double (*F)(double), double, double, int);
double QuadratureRuleGauss(double (*F)(double), double, double, int, int);
double F(double);

