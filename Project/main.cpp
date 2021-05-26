#include "montecarlo.h"

int main(int argc, char *argv[])
{
    int L = 10;        //number of batches
    int m = 50;        //number of dimension
    int n = 100;       // Number of trials
    double q = 0;      // dividend yield
    double S = 100.00; // initial Option price
    double K = 100.00; // Strike price
    double r = 0.10;   // Risk-free rate
    double v = 0.20;   //volatility
    double T = 1.00;   //expiration time;
    if (argc > 1)
    {
        sscanf(argv[1], "%d", &n);
    }
    for (int i = n; i < 10E5; i *= 10) //up until 100000 samples, as 10E6
                                       //and higher needs significantly more time
    {
        Asian_vanilla_call_option *option = new Asian_vanilla_call_option(i, q, S, K, r, v, T, L, m);
        auto start = chrono::steady_clock::now();
        option->quasi_monte_carlo_call_price();
        auto end = chrono::steady_clock::now();
        auto diff = end - start;
        cout << "Sample size: " << i << endl;
        cout << "Total no of trials: " << i * L << endl;
        cout << "Estimated Price by Quasi Monte-Carlo: " << option->call_price << endl;
        cout << "Estimated Standard Error: " << option->SE << endl;
        cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
             << option->call_price + option->SE * 1.96 << "]" << endl;
        cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
        cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;
        cout << "**************** end ********************" << endl;
    }
}