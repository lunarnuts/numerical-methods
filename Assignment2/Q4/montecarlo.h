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

using namespace std;

class European_vanilla_call_option
{
public:
    European_vanilla_call_option(const int &n,     // Number of divisions
                                 const double &q,  // dividend yield
                                 const double &S,  // initial Option price
                                 const double &K,  // Strike price
                                 const double &r,  // Risk-free rate
                                 const double &v,  //volatility
                                 const double &T); //expiration time;

    double option_price_call_black_scholes(); // time to maturity
    pair<double, double> monte_carlo_initial(int n);
    void monte_carlo_with_importance_sampling();
    double call_price = 0.0;
    double SD = 0;
    double SE = 0.0;
    double rho = 0;

private:
    int n;        // Number of trials
    double q;     // dividend yield
    double S;     // initial Option price
    double K;     // Strike price
    double r;     // Risk-free rate
    double v;     //volatility
    double T;     //expiration time
    double muHat; //selected mu
    double bsm;
    double max(double a, double b);
    double N(const double &z);
    double get_gbm(); //gaussian RV generator
    double get_gbm(const double &stddev);
    double get_gbm(const double &moydev, const double &stddev);
    double get_St(double mu, double stddev);
    pair<double, double> calculate_mu();
};
#endif