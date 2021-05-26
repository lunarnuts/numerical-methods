#include "montecarlo.h"
unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count(); //setup seed for RV generator
default_random_engine generator;

European_vanilla_call_option::European_vanilla_call_option(const int &n,    // Number of divisions
                                                           const double &q, // dividend yield
                                                           const double &S, // initial Option price
                                                           const double &K, // Strike price
                                                           const double &r, // Risk-free rate
                                                           const double &v, //volatility
                                                           const double &T) //expiration time;
{
    this->n = n;
    this->q = q;
    this->S = S;
    this->K = K;
    this->r = r;
    this->v = v;
    this->T = T;
};

double European_vanilla_call_option::max(double a, double b)
{
    return (b < a) ? a : b;
}

double European_vanilla_call_option::N(const double &z)
{
    if (z > 6.0)
    {
        return 1.0;
    }; // this guards against overflow
    if (z < -6.0)
    {
        return 0.0;
    };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0)
        n = 1.0 - n;
    return n;
};

double European_vanilla_call_option::get_gbm()
{
    double x = 0.0;
    double y = 0.0;
    double euclid_sq = 0.0;
    do
    {
        x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        euclid_sq = x * x + y * y;
    } while (euclid_sq >= 1.0);

    return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}

double European_vanilla_call_option::option_price_call_black_scholes()
{
    double time_sqrt = sqrt(T);
    double d1 = (log(S / K) + r * T) / (v * time_sqrt) + 0.5 * v * time_sqrt;
    double d2 = d1 - (v * time_sqrt);
    return S * N(d1) - K * exp(-r * T) * N(d2);
};

void European_vanilla_call_option::monte_carlo_call_price()
{
    double delta_R = (T * (r - q - 0.5 * v * v));
    double delta_SD = v * sqrt(T);
    double RV = get_gbm();
    double S_cur = S * exp(delta_R + delta_SD * RV);
    double temp = max(S_cur - K, 0.0);
    double sum = temp;
    double sum_sqr = temp * temp;
    double SE = 0;
    for (int i = 2; i <= n; i++)
    {
        RV = get_gbm();
        S_cur = S * exp(delta_R + delta_SD * RV);
        S_cur = max(S_cur - K, 0.0);
        sum += S_cur;
        sum_sqr += S_cur * S_cur;
    }

    this->call_price = sum / n * exp(-r * T);
    this->SE = sqrt((sum_sqr / n - sum * sum / (n * n)) / (n - 1));
}
