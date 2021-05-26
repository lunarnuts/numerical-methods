#include "montecarlo.h"
unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count(); //setup seed for RV generator
std::mt19937 generator;
std::normal_distribution<double> normal(0, 1.0);
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
double European_vanilla_call_option::get_gbm(const double &stddev)
{

    return normal(generator) * stddev;
}
double European_vanilla_call_option::get_gbm(const double &mean, const double &stddev)
{
    return get_gbm(stddev) + mean;
}

double European_vanilla_call_option::option_price_call_black_scholes()
{
    double time_sqrt = sqrt(T);
    double d1 = (log(S / K) + (r - q) * T) / (v * time_sqrt) + 0.5 * v * time_sqrt;
    double d2 = d1 - (v * time_sqrt);
    double result = S * N(d1) * exp(-r * T);
    this->bsm = result;
    return result;
};
double European_vanilla_call_option::get_St(double mu, double stddev)
{
    return S * exp(mu + v * sqrt(T) * get_gbm(0, 1));
}
pair<double, double> European_vanilla_call_option::calculate_mu()
{
    double mu = (T * (r - q - 0.5 * v * v));
    cout << " Mu is " << mu << endl;
    double delta_SD = v * sqrt(T);
    double RV = get_gbm();

    double W = 0; //likelyhood ratio or Weight

    vector<pair<double, double>> SE;
    vector<pair<double, double>> MuandWeight;
    int m = 10000;
    for (int j = 0; j < 100; j++)
    {
        double Ysum = 0;
        double YYsum = 0;
        double muHat = mu + (1 / (double)1000 * j);
        double A = exp((mu * mu - muHat * muHat) * (-1 / (2 * v * v * T)));
        double C = 2 * (muHat - mu) * (-1 / (2 * v * v * T));
        for (int i = 0; i < m; i++)
        {
            double Xt = get_St(muHat, v * v * T);
            W = A * exp(C * log(Xt / S));
            double temp = (Xt > K) ? Xt * W * exp(-r * T) : 0;
            Ysum += temp;
            YYsum += temp * temp;
        }
        double Ym = Ysum / m;
        double YYm = YYsum / m;
        double sigma = sqrt((YYm - Ym * Ym));
        double se = sigma / sqrt(m - 1);
        MuandWeight.push_back(pair<double, double>{muHat, Ym});
        SE.push_back(pair<double, double>{se, muHat});
    }
    double min_se = 0;
    for (int i = 0; i < 100; i++)
    {
        if (SE[i].first < SE[min_se].first)
        {
            min_se = i;
        }
    }
    cout << "mu is " << SE[min_se].second << endl;
    return MuandWeight[min_se];
}

void European_vanilla_call_option::monte_carlo_with_importance_sampling()
{
    pair<double, double> temp = calculate_mu();
    this->muHat = temp.first;
    double mu = (T * (r - q - 0.5 * v * v));
    double delta_SD = v * sqrt(T);
    double RV = get_gbm();
    double X = 0;
    double W = 0;
    double FX = 0;
    double FXsum = 0;
    double FXXsum = 0;
    double A = exp((mu * mu - muHat * muHat) * (-1 / (2 * v * v * T)));
    double C = 2 * (muHat - mu) * (-1 / (2 * v * v * T));
    for (int i = 0; i < n; i++)
    {
        double X = get_St(muHat, 1.0);
        W = A * exp(C * log(X / S));
        FX = X > K ? X * W : 0;
        FXsum += FX;
        FXXsum += FX * FX;
    }
    double FXm = FXsum / n;
    double FXXm = FXXsum / n;
    double varFX = FXXm - FXm * FXm;
    double sd = sqrt(varFX);

    this->SD = sd;
    this->SE = sd / sqrt(n - 1);
    this->call_price = FXm * exp(-r * T);
}
