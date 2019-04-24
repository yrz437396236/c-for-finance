#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"          // this defines the normal distribution from Odegaard's files
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price, downtick_prob, notick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

float max(float a, float b) {
	return (b < a) ? a : b;
}

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K * exp(-r * time)*N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S * N(d1) - K * exp(-r * time)*N(d2);
};

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double european_call_option(int k, int i, double **memo) {
	if (memo[k][i+k] != -1)
	{
		return memo[k][i+k];
	}
	else
	{
		if (k == no_of_divisions)
			return max(0.0, (initial_stock_price*pow(up_factor, ((float)i))) - strike_price);
		else
			memo[k][i + k] = (((uptick_prob*european_call_option(k + 1, i + 1, memo)) +
							 (notick_prob*european_call_option(k + 1, i, memo)) +
							 (downtick_prob*european_call_option(k + 1, i - 1, memo))) / R);
		return memo[k][i+k];
	}
}

double european_put_option(int k, int i, double **memo) {
	if (memo[k][i+k] != -1)
	{
		return memo[k][i+k];
	}
	else
	{
		if (k == no_of_divisions)
			return max(0.0, strike_price - (initial_stock_price*pow(up_factor, ((float)i))));
		else
			memo[k][i+k] = (((uptick_prob*european_put_option(k + 1, i + 1, memo)) +
						 (notick_prob*european_put_option(k + 1, i, memo)) +
				         (downtick_prob*european_put_option(k + 1, i - 1, memo))) / R);
		return memo[k][i+k];
	}
}

int main(int argc, char* argv[])
{
	clock_t start = clock();
	//sscanf_s(argv[1], "%f", &expiration_time);
	//sscanf_s(argv[2], "%d", &no_of_divisions);
	//sscanf_s(argv[3], "%f", &risk_free_rate);
	//sscanf_s(argv[4], "%f", &volatility);
	//sscanf_s(argv[5], "%f", &initial_stock_price);
	//sscanf_s(argv[6], "%f", &strike_price);
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);
	double **array1 = new double*[no_of_divisions+1];
	for (int i = 0; i <= no_of_divisions; i++)
		array1[i] = new double[2* no_of_divisions];
	for (int i = 0; i <= no_of_divisions; i++)
		for (int j = 0; j <= (2 * no_of_divisions - 1); j++)
		{
			array1[i][j] = -1;
		}
	double **array2 = new double*[no_of_divisions + 1];
	for (int i = 0; i <= no_of_divisions; i++)
		array2[i] = new double[2 * no_of_divisions];
	for (int i = 0; i <= no_of_divisions; i++)
		for (int j = 0; j <= (2 * no_of_divisions - 1); j++)
		{
			array2[i][j] = -1;
		}
	//pow(R,0.5) pow(up_factor, 0.5)
	up_factor = exp(volatility*sqrt((2 * expiration_time) / ((double)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((double)no_of_divisions));
	uptick_prob = pow((sqrt(R) - (1 / (sqrt(up_factor)))) / ((sqrt(up_factor)) - (1 / (sqrt(up_factor)))), 2);
	downtick_prob = pow((sqrt(up_factor) - (sqrt(R))) / ((sqrt(up_factor)) - (1 / (sqrt(up_factor)))), 2);
	notick_prob = 1 - uptick_prob - downtick_prob;
	cout << "(Memoized)Recursive Trinomial European Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = european_call_option(0, 0, array1);
	cout << "Trinomial Price of an European Call Option = " << call_price << endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	double put_price = european_put_option(0, 0, array2);
	cout << "Trinomial Price of an European Put Option = " << put_price << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " = " << strike_price * exp(-risk_free_rate * expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	clock_t end = clock();
	cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	system("pause");
}