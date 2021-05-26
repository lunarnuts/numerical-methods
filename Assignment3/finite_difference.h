#ifndef montecarlo
#define montecarlo

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <vector>
#include <math.h>

using namespace std;

class Double_Barrier_knock_out_call_option
{
public:
    Double_Barrier_knock_out_call_option(const double &L,  // Lower Boundary
                                         const double &U,  // Upper Boundary
                                         const int &N,     // Temporal discretization steps
                                         const int &M,     // Spatial discretization steps
                                         const double &q,  // dividend yield
                                         const double &K,  // Strike price
                                         const double &r,  // Risk-free rate
                                         const double &v,  // volatility
                                         const double &T); // number of timesteps

    void explicit_euler_scheme();
    void implicit_euler_scheme();
    void crank_nicolson();
    void extrapolation();
    double call_price = 0.0;
    double SE = 0.0;
    double option_price_call_black_scholes(double S);

private:
    double L = 0;   // Lower Boundary
    double U = 0;   // Upper Boundary
    int N = 0;      // Temporal discretization steps
    int M = 0;      // Spatial discretization steps
    double q = 0.0; // dividend yield
    double K = 0.0; // initial Option price
    double r = 0.0; // Risk-free rate
    double v = 0.0; // volatility
    double T = 0.0; // time to maturity
    void thomas_algorithm(const vector<double> &a,
                          const vector<double> &b,
                          const vector<double> &c,
                          const vector<double> &d,
                          vector<double> &f);
    double Norm(const double &z);
};
#endif