#include "montecarlo.h"
unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count(); //setup seed for RV generator
std::mt19937 generator;
std::normal_distribution<double> normal(0, 1.0);

Asian_vanilla_call_option::Asian_vanilla_call_option(const int &n,    // n of divisions
                                                     const double &q, // dividend yield
                                                     const double &S, // initial Option price
                                                     const double &K, // Strike price
                                                     const double &r, // Risk-free rate
                                                     const double &v, //volatility
                                                     const double &T, //expiration time;
                                                     const double &L, // n of Les
                                                     const double &m) // n of timesteps)
{
    this->n = n;
    this->L = L;
    this->m = m;
    this->q = q;
    this->S = S;
    this->K = K;
    this->r = r;
    this->v = v;
    this->T = T;
};

double Asian_vanilla_call_option::max(double a, double b)
{
    return (b < a) ? a : b;
}

double Asian_vanilla_call_option::N(const double &z)
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

double Asian_vanilla_call_option::get_gbm()
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
double Asian_vanilla_call_option::get_sobol(int index, int dimension, double number)
{
    const double s = sobol::sample(index + 1, dimension);
    //double inter = s + number - floor(s + number);
    return get_gbm(0, s) / s;
}
double Asian_vanilla_call_option::get_gbm(const double &stddev)
{

    return get_gbm() * stddev;
}
double Asian_vanilla_call_option::get_gbm(const double &mean, const double &stddev)
{
    return get_gbm(stddev) + mean;
}

double Asian_vanilla_call_option::option_price_call_black_scholes()
{
    double time_sqrt = sqrt(T);
    double d1 = (log(S / K) + (r - q) * T) / (v * time_sqrt) + 0.5 * v * time_sqrt;
    double d2 = d1 - (v * time_sqrt);
    return S * exp(-q * T) * N(d1) - K * exp(-r * T) * N(d2);
};

void Asian_vanilla_call_option::quasi_monte_carlo_call_price()
{
    double y = 0;  // Yi
    double y2 = 0; //Yi^2

    double deltaT = T / (double)m; //timestep
    double drift = (r - q - v * v * 0.5) * deltaT;

    for (int j = 0; j < L; j++) //generate L Les
    {
        double x = 0;
        double U[m]; //uniform dist vector
        for (int i = 0; i < m; i++)
        {
            srand(((int)time(0) * 10000 * (j + 1)) * L);
            U[i] = normal(generator);
        }
        for (int i = 0; i < n; i++) //generate n samples
        {
            double Price[m];
            double RV = get_sobol(i, 0, U[0]); //generate quasi random variable using  sobol sequence
            Price[0] = S * exp(drift + v * sqrt(deltaT) * RV);

            for (int d = 1; d < m; d++) //generate path
            {
                RV = get_sobol(i, d, U[d]);
                Price[d] = Price[d - 1] * exp(drift + v * sqrt(deltaT) * RV);
            }
            double sum = 0;

            for (int k = 0; k < m; k++)
            {
                sum += Price[k];
            }
            double Price_avg = sum / (double)(m);
            Price_avg = max(Price_avg - K, 0);
            x += Price_avg;
        }
        x = x / (double)(n);
        y += x;
        y2 += x * x;
    }
    double payoff = y / (double)(L);
    this->call_price = payoff * exp(-r * T);
    this->SE = sqrt(((y2 / L) - payoff * payoff) / (L - 1)); //
};
