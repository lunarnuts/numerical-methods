#include "finite_difference.h"
// Double Barrier knock-out option pricing by Finite Differene Method.
// Created by Aigerim Tursynbekova (C) 2021.
// Please do not copy to submit as your own work.

int main(int argc, char *argv[])
{
     double L = 80;    // Lower Boundary
     double U = 120;   // Upper Boundary
     int M = 800;      // Spatial discretization steps
     int N = 10000;    // Temporal discretization steps
     double q = 0.00;  // dividend yield
     double K = 100.0; // Strike price
     double r = 0.01;  // Risk-free rate
     double v = 0.20;  // volatility
     double T = 0.25;  // time to maturity

     Double_Barrier_knock_out_call_option *option = new Double_Barrier_knock_out_call_option(L, U, N, M, q, K, r, v, T);
     cout << "Temporal discretization steps: " << N << endl;
     cout << "Spatial discretization steps: " << M << endl;
     cout << "Lower Boundary: " << L << ", Xmin: " << log(L) << endl;
     cout << "Upper Boundary: " << U << ", Xmax: " << log(U) << endl;
     cout << "Strike price: " << K << endl;
     cout << "Risk-free interest rate: " << r * 100 << "%" << endl;
     cout << "Dividend yield: " << q * 100 << "%" << endl;
     cout << "Volatility: " << v * 100 << "%" << endl;
     cout << "Time to maturity: " << T << " years" << endl;
     cout << "******** Explicit Euler Schema FDE *********" << endl;
     auto start = chrono::steady_clock::now();
     option->explicit_euler_scheme();
     auto end = chrono::steady_clock::now();
     auto diff = end - start;
     double BSM = 3.913; //option->option_price_call_black_scholes(100); // for x441
     cout
         << "Estimated Price : " << option->call_price << endl;
     cout << "Estimated BSM price: " << BSM << endl;
     cout << "Estimated Standard Error: " << BSM - option->call_price << endl;
     cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
          << option->call_price + option->SE * 1.96 << "]" << endl;
     cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
     cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;

     cout << "******** Implicit Euler Schema FDE *********" << endl;
     start = chrono::steady_clock::now();
     option->implicit_euler_scheme();
     end = chrono::steady_clock::now();
     diff = end - start;
     cout << "Estimated Price : " << option->call_price << endl;
     cout << "Estimated BSM price: " << BSM << endl;
     cout << "Estimated Standard Error: " << BSM - option->call_price << endl;
     cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
          << option->call_price + option->SE * 1.96 << "]" << endl;
     cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
     cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;

     cout << "******** Crank-Nicolson Schema FDE *********" << endl;
     start = chrono::steady_clock::now();
     option->crank_nicolson();
     end = chrono::steady_clock::now();
     diff = end - start;
     cout << "Estimated Price : " << option->call_price << endl;
     cout << "Estimated BSM price: " << BSM << endl;
     cout << "Estimated Standard Error: " << BSM - option->call_price << endl;
     cout << "The 95% CI: [" << option->call_price - option->SE * 1.96 << ", "
          << option->call_price + option->SE * 1.96 << "]" << endl;
     cout << "Computational time: " << chrono::duration<double>(diff).count() << " seconds" << endl;
     cout << "Efficiency: " << option->SE * option->SE * chrono::duration<double>(diff).count() << endl;
     cout << "**************** end ********************" << endl;
}