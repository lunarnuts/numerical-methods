#include "montecarlo.h"

int main(int argc, char *argv[])
{
    int n = 1000;                      // Number of divisions
    double q = 0.0232;                 // dividend yield
    double S = 1868.99;                // initial Option price
    double K = 1870.0;                 // Strike price
    double r = 0.003866;               // Risk-free rate
    double v = 0.2979;                 //volatility
    double T = (double)1 / (double)52; //expiration time;
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
    cout << "Sample size: " << n << endl;
    cout << "Estimated Price by BSM: " << bsm << endl;
    cout << "Estimated Price by Monte-Carlo: " << option->call_price << endl;
    cout << "Estimated Standard Error: " << option->SE << endl;
    cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
         << option->call_price + option->SE * 1.96 << "]" << endl;
    cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
    cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;
}