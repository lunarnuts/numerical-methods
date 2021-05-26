#include "finite_difference.h"
// Double Barrier knock-out option pricing by Finite Differene Method.
// Created by lunarnuts(C) 2021.
// Please do not copy to submit as your own work.

// Used https://www.quantstart.com/articles/C-Explicit-Euler-Finite-Difference-Method-for-Black-Scholes/
// as a reference

// Class constructor
Double_Barrier_knock_out_call_option::Double_Barrier_knock_out_call_option(const double &L, // Lower Boundary
                                                                           const double &U, // Upper Boundary
                                                                           const int &N,    // Spatial discretization steps
                                                                           const int &M,    // Temporal discretization steps
                                                                           const double &q, // Dividend yield
                                                                           const double &K, // Strike price
                                                                           const double &r, // Risk-free rate
                                                                           const double &v, // Volatility
                                                                           const double &T)
{
    this->L = L;
    this->U = U;
    this->N = N;
    this->M = M;
    this->q = q;
    this->K = K;
    this->r = r;
    this->v = v;
    this->T = T;
};

// Explicit Euler Scheme FDM
void Double_Barrier_knock_out_call_option::explicit_euler_scheme()
{
    vector<vector<double>> result_table(N, vector<double>(M, 0.0));
    vector<double> Ai(M - 1, 0.0);
    vector<double> Bi(M, 0.0);
    vector<double> Ci(M - 1, 0.0);
    vector<double> Xi(M, 0.0);
    vector<double> Fi(M, 0.0); //Boundary condition vector
    vector<vector<double>> Imatrix(M, vector<double>(M, 0.0));

    //fill in identity matrix
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (j == i)
                Imatrix[i][j] = 1;
        }
    }
    double Xmin = log(L);
    double Xmax = log(U);
    double deltaT = T / N;
    double deltaX = (Xmax - Xmin) / (M + 1);
    double mu = (r - q - 0.5 * v * v);
    double A = v * v / (2 * deltaX * deltaX) - mu / (2 * deltaX);
    double B = r + v * v / (deltaX * deltaX);
    double C = v * v / (2 * deltaX * deltaX) + mu / (2 * deltaX);

    //boundary conditions
    Fi[0] = -A * Xmin;
    Fi[M - 1] = -C * Xmax;

    // calculate A,B,C coeficients and vector X
    for (int i = 0; i < M; i++)
    {
        Xi[i] = Xmin + i * deltaX;
        if (i < M - 1)
        {
            Ai[i] = -A * deltaT;
            Ci[i] = -C * deltaT;
        }
        Bi[i] = B * deltaT;
    }
    result_table[0] = Xi;
    for (int i = 0; i < N - 1; i++)
    {
        vector<double> Di(M, 0.0);
        vector<double> Fn(M, 0.0);
        Fn = Xi; //result_table[i];
        Di[0] = Fn[0] - (-B * Fn[0] + C * Fn[1]) * i * deltaT - deltaT * i * Fi[0];
        Di[M - 1] = Fn[M - 1] - (A * Fn[M - 2] - B * Fn[M - 1]) * (i)*deltaT - deltaT * (i)*Fi[M - 1];
        for (int j = 1; j < M - 1; j++)
        {
            Di[j] = Fn[j] - (A * Fn[j - 1] - B * Fn[j] + C * Fn[j + 1]) * (i)*deltaT - deltaT * i * Fi[j];
        }

        result_table[i + 1] = Di;
    }

    for (int i = 0; i < N; i++)
    {
        cout << "T = " << i * deltaT << " : " << max(exp(result_table[i][484]) * exp(-q * T) - K * exp(-r * T), 0.0) << ", " << endl;
    }
    this->call_price = max(exp(result_table[N - 1][485]) * exp(-q * T) - K * exp(-r * T), 0.0);
};

// Explicit Euler Scheme FDM
void Double_Barrier_knock_out_call_option::implicit_euler_scheme()
{
    vector<vector<double>> result_table(N, vector<double>(M, 0.0));
    vector<double> Ai(M - 1, 0.0);
    vector<double> Bi(M, 0.0);
    vector<double> Ci(M - 1, 0.0);
    vector<double> Xi(M, 0.0);
    vector<double> Fi(M, 0.0); //Boundary condition vector
    vector<vector<double>> Imatrix(M, vector<double>(M, 0.0));

    //fill in identity matrix
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (j == i)
                Imatrix[i][j] = 1;
        }
    }
    double Xmin = log(L);
    double Xmax = log(U);
    double deltaT = T / N;
    double deltaX = (Xmax - Xmin) / (M + 1);
    double mu = (r - q - 0.5 * v * v);
    double A = v * v / (2 * deltaX * deltaX) - mu / (2 * deltaX);
    double B = r + v * v / (deltaX * deltaX);
    double C = v * v / (2 * deltaX * deltaX) + mu / (2 * deltaX);

    //boundary conditions
    Fi[0] = -A * Xmin;
    Fi[M - 1] = -C * Xmax;

    // calculate A,B,C coeficients and vector X
    for (int i = 0; i < M; i++)
    {
        Xi[i] = Xmin + i * deltaX;
        if (i < M - 1)
        {
            Ai[i] = -A * deltaT;
            Ci[i] = -C * deltaT;
        }
        Bi[i] = 1 + B * deltaT;
    }
    result_table[0] = Xi;
    for (int i = 0; i < N - 1; i++)
    {
        vector<double> Di(M, 0.0);
        vector<double> Fn(M, 0.0);
        Fn = Xi; //result_table[i];
        Di[0] = Fn[0] - (-B * Fn[0] + C * Fn[1]) * i * deltaT - deltaT * i * Fi[0];
        Di[M - 1] = Fn[M - 1] - (A * Fn[M - 2] - B * Fn[M - 1]) * (i)*deltaT - deltaT * (i)*Fi[M - 1];
        for (int j = 1; j < M - 1; j++)
        {
            Di[j] = (Fn[j] - deltaT * i * Fi[j]);
        }
        thomas_algorithm(Ai, Bi, Ci, Di, Fn);
        result_table[i + 1] = Fn;
    }

    for (int i = 0; i < N; i++)
    {
        cout << "T = " << i * deltaT << " : " << max(exp(result_table[i][512]) * exp(-q * T) - K * exp(-r * T), 0.0) << ", " << endl;
    }
    this->call_price = max(exp(result_table[N - 1][512]) * exp(-q * T) - K * exp(-r * T), 0.0);
}

// Crank-Nicolson Scheme FDM
void Double_Barrier_knock_out_call_option::crank_nicolson()
{

    vector<vector<double>> result_table(N, vector<double>(M, 0.0));
    vector<double> Ai(M - 1, 0.0);
    vector<double> Bi(M, 0.0);
    vector<double> Ci(M - 1, 0.0);
    vector<double> Xi(M, 0.0);
    vector<double> Fi(M, 0.0); //Boundary condition vector
    vector<vector<double>> Imatrix(M, vector<double>(M, 0.0));

    //fill in identity matrix
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (j == i)
                Imatrix[i][j] = 1;
        }
    }
    double Xmin = log(L);
    double Xmax = log(U);
    double deltaT = T / N;
    double deltaX = (Xmax - Xmin) / (M + 1);
    double mu = (r - q - 0.5 * v * v);
    double A = v * v / (2 * deltaX * deltaX) - mu / (2 * deltaX);
    double B = r + v * v / (deltaX * deltaX);
    double C = v * v / (2 * deltaX * deltaX) + mu / (2 * deltaX);

    //boundary conditions
    Fi[0] = -A * Xmin;
    Fi[M - 1] = -C * Xmax;

    // calculate A,B,C coeficients and vector X
    for (int i = 0; i < M; i++)
    {
        Xi[i] = Xmin + i * deltaX;
        if (i < M - 1)
        {
            Ai[i] = -A * 0.5 * deltaT;
            Ci[i] = -C * 0.5 * deltaT;
        }
        Bi[i] = 1 + B * 0.5 * deltaT;
    }
    result_table[0] = Xi;
    for (int i = 0; i < N - 1; i++)
    {
        vector<double> Di(M, 0.0);
        vector<double> Fn(M, 0.0);
        Fn = Xi; //result_table[i];
        Di[0] = Fn[0] - (-B * Fn[0] + C * Fn[1]) * i * deltaT - deltaT * i * Fi[0];
        Di[M - 1] = Fn[M - 1] - (A * Fn[M - 2] - B * Fn[M - 1]) * (i)*deltaT - deltaT * (i)*Fi[M - 1];
        for (int j = 1; j < M - 1; j++)
        {
            Di[j] = Fn[j] - 0.5 * (A * Fn[j - 1] - B * Fn[j] + C * Fn[j + 1]) * (i)*deltaT - deltaT * i * Fi[j];
        }
        thomas_algorithm(Ai, Bi, Ci, Di, Fn);
        result_table[i + 1] = Fn;
    }

    for (int i = 0; i < N; i++)
    {
        cout << "T = " << i * deltaT << " : " << max(exp(result_table[i][498]) * exp(-q * T) - K * exp(-r * T), 0.0) << ", " << endl;
    }
    this->call_price = max(exp(result_table[N - 1][498]) * exp(-q * T) - K * exp(-r * T), 0.0);
}

// Extrapolation
void Double_Barrier_knock_out_call_option::extrapolation()
{
    int N = 1;
    {
        vector<vector<double>> result_table(N, vector<double>(M, 0.0));
        vector<double> Ai(M - 1, 0.0);
        vector<double> Bi(M, 0.0);
        vector<double> Ci(M - 1, 0.0);
        vector<double> Xi(M, 0.0);
        vector<double> Fi(M, 0.0); //Boundary condition vector
        vector<vector<double>> Imatrix(M, vector<double>(M, 0.0));

        //fill in identity matrix
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < M; j++)
            {
                if (j == i)
                    Imatrix[i][j] = 1;
            }
        }
        double Xmin = log(L);
        double Xmax = log(U);
        double deltaT = T / N;
        double deltaX = (Xmax - Xmin) / (M + 1);
        double mu = (r - q - 0.5 * v * v);
        double A = v * v / (2 * deltaX * deltaX) - mu / (2 * deltaX);
        double B = r + v * v / (deltaX * deltaX);
        double C = v * v / (2 * deltaX * deltaX) + mu / (2 * deltaX);

        //boundary conditions
        Fi[0] = -A * Xmin;
        Fi[M - 1] = -C * Xmax;

        // calculate A,B,C coeficients and vector X
        for (int i = 0; i < M; i++)
        {
            Xi[i] = Xmin + i * deltaX;
            if (i < M - 1)
            {
                Ai[i] = -A * deltaT;
                Ci[i] = -C * deltaT;
            }
            Bi[i] = 1 + B * deltaT;
        }
        result_table[0] = Xi;
        for (int i = 0; i < N - 1; i++)
        {
            vector<double> Di(M, 0.0);
            vector<double> Fn(M, 0.0);
            Fn = Xi; //result_table[i];
            Di[0] = Fn[0] - (-B * Fn[0] + C * Fn[1]) * i * deltaT - deltaT * i * Fi[0];
            Di[M - 1] = Fn[M - 1] - (A * Fn[M - 2] - B * Fn[M - 1]) * (i)*deltaT - deltaT * (i)*Fi[M - 1];
            for (int j = 1; j < M - 1; j++)
            {
                Di[j] = (Fn[j] - deltaT * i * Fi[j]);
            }
            thomas_algorithm(Ai, Bi, Ci, Di, Fn);
            result_table[i + 1] = Fn;
        }
    }
    /*
    for (int i = 0; i < N; i++)
    {
        cout << "T = " << i * deltaT << " : " << max(exp(result_table[i][512]) * exp(-q * T) - K * exp(-r * T), 0.0) << ", " << endl;
    }
    this->call_price = max(exp(result_table[N - 1][512]) * exp(-q * T) - K * exp(-r * T), 0.0);
    */
}
// Thomas algorithm to multiply tridiag matrix
void Double_Barrier_knock_out_call_option::thomas_algorithm(const vector<double> &a,
                                                            const vector<double> &b,
                                                            const vector<double> &c,
                                                            const vector<double> &d,
                                                            vector<double> &f)
{
    size_t N = d.size();
    std::vector<double> c_star(N, 0.0);
    std::vector<double> d_star(N, 0.0);
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for (int i = 1; i < N; i++)
    {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }
    for (int i = N - 1; i-- > 0;)
    {
        f[i] = d_star[i] - c_star[i] * d[i + 1];
    }
};
double Double_Barrier_knock_out_call_option::option_price_call_black_scholes(double S)
{
    double time_sqrt = sqrt(T);
    double d1 = (log(S / K) + (r - q) * T) / (v * time_sqrt) + 0.5 * v * time_sqrt;
    double d2 = d1 - (v * time_sqrt);
    return S * exp(-r * T) * Norm(d1) - K * exp(-r * T) * Norm(d2);
};
double Double_Barrier_knock_out_call_option::Norm(const double &z)
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