#include "montecarlo.h"

int main(int argc, char *argv[])
{
    int n = 1000;                      // Number of divisions
    double q = 0.02;                   // dividend yield
    double S = 2000;                   // initial Option price
    double K = 2200;                   // Strike price
    double r = 0.005;                  // Risk-free rate
    double v = 0.3;                    //volatility
    double T = (double)1 / (double)12; //expiration time;
    if (argc > 1)
    {
        sscanf(argv[1], "%d", &n);
    }
    for (n = 1000; n < 10E7; n *= 10)
    {
        European_vanilla_call_option *option = new European_vanilla_call_option(n, q, S, K, r, v, T);
        double bsm = option->option_price_call_black_scholes();
        auto start = chrono::steady_clock::now();
        option->monte_carlo_call_price();
        auto end = chrono::steady_clock::now();
        auto diff = end - start;
        cout << "Sample size: " << n << endl;
        cout << "Estimated Price by BSM: " << bsm << endl;
        cout << "Estimated Price by Naive Monte-Carlo: " << option->call_price << endl;
        cout << "Estimated Standard Error: " << option->SE << endl;
        cout << "Estimated Standard Deviation: " << option->SD << endl;
        cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
             << option->call_price + option->SE * 1.96 << "]" << endl;
        cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
        cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;
        cout << "***********************************" << endl;
    }
}