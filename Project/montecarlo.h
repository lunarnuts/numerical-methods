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
#include "sobol.h"

using namespace std;

class Asian_vanilla_call_option
{
public:
    Asian_vanilla_call_option(const int &n,     // Number of divisions
                              const double &q,  // dividend yield
                              const double &S,  // initial Option price
                              const double &K,  // Strike price
                              const double &r,  // Risk-free rate
                              const double &v,  //volatility
                              const double &T,  //expiration time;
                              const double &L,  // number of batches
                              const double &m); // number of timesteps

    double option_price_call_black_scholes(); // time to maturity
    void quasi_monte_carlo_call_price();
    double call_price = 0.0;
    double SE = 0.0;
    double SD = 0.0;

private:
    int L;
    int m;
    int n;    // Number of trials
    double q; // dividend yield
    double S; // initial Option price
    double K; // Strike price
    double r; // Risk-free rate
    double v; //volatility
    double T; //expiration time
    double max(double a, double b);
    double N(const double &z);
    double get_gbm(); //gaussian RV generator
    double get_gbm(const double &stddev);
    double get_gbm(const double &mean, const double &stddev);
    double get_sobol(int index, int dimension, double number);
};
#endif