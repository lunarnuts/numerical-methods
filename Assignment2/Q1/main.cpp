#include "montecarlo.h"

int main(int argc, char *argv[])
{
    int n = 1000;    // Number of divisions
    double q = 0.0;  // dividend yield
    double S = 100;  // initial Option price
    double K = 20;   // Strike price
    double r = 0.02; // Risk-free rate
    double v = 0.2;  //volatility
    double T = 5.0;  //(double)1 / (double)12; //expiration time;

    if (argc > 1)
    {
        sscanf(argv[1], "%d", &n);
    }
    European_vanilla_call_option *option = new European_vanilla_call_option(n, q, S, K, r, v, T);
    double bsm = option->option_price_call_black_scholes();
    auto start = chrono::steady_clock::now();
    option->monte_carlo_call_price();
    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Estimated Price by BSM: " << bsm << endl;
}