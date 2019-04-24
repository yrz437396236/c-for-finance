#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include "normdist.h"
#include <math.h>
#include <chrono>

#define pi 3.141592653589793
using namespace std;

double risk_free_rate, strike_price, barrier_price, initial_stock_price;
double expiration_time, volatility;
double adjusted_call_price, adjusted_put_price;
int no_of_divisions, no_of_trials;
double call, put;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double max(double a, double b) {
	return (b < a) ? a : b;
}

double get_uniform()
{
	std::uniform_real_distribution<double>distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

void european_down_and_out_option_continous() {
	//generate 4 paths of price and calculate the average
	int k;
	double delta_T = expiration_time / ((double)no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility * sqrt(delta_T);

	double S1 = initial_stock_price;
	double S2 = initial_stock_price;
	double S3 = initial_stock_price;
	double S4 = initial_stock_price;
	int barrier1 = 1; int barrier2 = 1; int barrier3 = 1; int barrier4 = 1;

	for (k = 1; k <= no_of_divisions; k++) {
		double x = get_uniform();
		double y = get_uniform();
		double a = sqrt(-2.0*log(x)) * cos(2 * pi*y);
		double b = sqrt(-2.0*log(x)) * sin(2 * pi*y);
		if (S1 <= barrier_price) {
			barrier1 = 0;
		}
		else {
			S1 = S1 * exp(delta_R + delta_SD * a);
		}
		if (S2 <= barrier_price) {
			barrier2 = 0;
		}
		else {
			S2 = S2 * exp(delta_R - delta_SD * a);
		}
		if (S3 <= barrier_price) {
			barrier3 = 0;
		}
		else {
			S3 = S3 * exp(delta_R + delta_SD * b);
		}
		if (S4 <= barrier_price) {
			barrier4 = 0;
		}
		else {
			S4 = S4 * exp(delta_R - delta_SD * b);
		}
	}
	call = (max(0.0, S1 - strike_price)*barrier1 + max(0.0, S2 - strike_price)*barrier2
		+ max(0.0, S3 - strike_price)*barrier3 + max(0.0, S4 - strike_price)*barrier4) / 4.0;
	put = (max(0.0, strike_price - S1)*barrier1 + max(0.0, strike_price - S2)*barrier2
		+ max(0.0, strike_price - S3)*barrier3 + max(0.0, strike_price - S4)*barrier4) / 4.0;
//if the price touch the barrier, then put zero in it
}

void simulation_adjusted() {
	double PC1, PC2, PC3, PC4;
	double R = (risk_free_rate - 0.5*volatility*volatility)*expiration_time;
	double SD = volatility * sqrt(expiration_time);
	double S1 = initial_stock_price;
	double S2 = initial_stock_price;
	double S3 = initial_stock_price;
	double S4 = initial_stock_price;
	double x = get_uniform();
	double y = get_uniform();
	double a = sqrt(-2.0*log(x)) * cos(2 * pi*y);
	double b = sqrt(-2.0*log(x)) * sin(2 * pi*y);
	S1 = S1 * exp(R + SD * a);
	if (initial_stock_price <= barrier_price || S1 <= barrier_price) {
		PC1 = 1;
	}
	else {
		PC1 = exp(-(2 * log(initial_stock_price / ((double)barrier_price))*log(S1 / ((double)barrier_price))) / (expiration_time*pow(volatility, 2)));
	}
	S2 = S2 * exp(R - SD * a);
	if (initial_stock_price <= barrier_price || S2 <= barrier_price) {
		PC2 = 1;
	}
	else {
		PC2 = exp(-(2 * log(initial_stock_price / ((double)barrier_price))*log(S2 / ((double)barrier_price))) / (expiration_time*pow(volatility, 2)));
	}
	S3 = S3 * exp(R + SD * b);
	if (initial_stock_price <= barrier_price || S3 <= barrier_price) {
		PC3 = 1;
	}
	else {
		PC3 = exp(-(2 * log(initial_stock_price / ((double)barrier_price))*log(S3 / ((double)barrier_price))) / (expiration_time*pow(volatility, 2)));
	}
	S4 = S4 * exp(R - SD * b);
	if (initial_stock_price <= barrier_price || S4 <= barrier_price) {
		PC4 = 1;
	}
	else {
		PC4 = exp(-(2 * log(initial_stock_price / ((double)barrier_price))*log(S4 / ((double)barrier_price))) / (expiration_time*pow(volatility, 2)));
	}

	adjusted_call_price = (max(0.0, S1 - strike_price)*(1 - PC1) + max(0.0, S2 - strike_price)*(1 - PC2) +
		max(0.0, S3 - strike_price)*(1 - PC3) + max(0.0, S4 - strike_price)*(1 - PC4)) / 4.0;
	adjusted_put_price = (max(0.0, strike_price - S1)*(1 - PC1) + max(0.0, strike_price - S2)*(1 - PC2) +
		max(0.0, strike_price - S3)*(1 - PC3) + max(0.0, strike_price - S4)*(1 - PC4)) / 4.0;
}


double option_price_put_black_scholes(const double& S,      
	const double& K,      
	const double& r,      
	const double& sigma,  
	const double& time)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K * exp(-r * time)*N(-d2) - S * N(-d1);
}

double option_price_call_black_scholes(const double& S,       
	const double& K,       
	const double& r,       
	const double& sigma,   
	const double& time)      
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S * N(d1) - K * exp(-r * time)*N(d2);
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; 
	if (z < -6.0) { return 0.0; };
	//prevent extreme result
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
}

double closed_form_down_and_out_european_call_option()
{
	double K = (2 * risk_free_rate) / (volatility*volatility);
	double A = option_price_call_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	double B = (barrier_price*barrier_price) / initial_stock_price;
	double C = pow(initial_stock_price / barrier_price, -(K - 1));
	double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D * C);
}

double closed_form_down_and_in_european_put_option()
{
	double S = initial_stock_price;
	double r = risk_free_rate;
	double T = expiration_time;
	double sigma = volatility;
	double H = barrier_price;
	double X = strike_price;

	double lambda = (r + ((sigma*sigma) / 2)) / (sigma*sigma);
	double temp = 2 * lambda - 2.0;
	double x1 = (log(S / H) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y = (log(H*H / (S*X)) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	double y1 = (log(H / S) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
	return (-S * N(-x1) + X * exp(-r * T)*N(-x1 + sigma * sqrt(T)) +
		S * pow(H / S, 2 * lambda)*(N(y) - N(y1)) -
		X * exp(-r * T)*pow(H / S, temp)*(N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

double closed_form_down_and_out_european_put_option()
{
	double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	double put_down_in = closed_form_down_and_in_european_put_option();
	return (vanilla_put - put_down_in);
}


int main(int argc, char* argv[])
{
	clock_t start = clock();
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_divisions);
	sscanf_s(argv[8], "%lf", &barrier_price);

	double continous_call_option = 0.0;
	double continous_put_option = 0.0;
	double call_option_price = 0.0;
	double put_option_price = 0.0;
	cout << "--------------------------------------" << endl;
	cout << "European Down—and-Out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;

	for (int i = 0; i < no_of_trials; i++) {
		european_down_and_out_option_continous();
		simulation_adjusted();

		continous_call_option = continous_call_option + call;
		continous_put_option = continous_put_option + put;
		call_option_price = call_option_price + adjusted_call_price;
		put_option_price = put_option_price + adjusted_put_price;
	}

	call_option_price = exp(-risk_free_rate * expiration_time)*(call_option_price / ((double)no_of_trials));
	put_option_price = exp(-risk_free_rate * expiration_time)*(put_option_price / ((double)no_of_trials));
	continous_call_option = exp(-risk_free_rate * expiration_time)*(continous_call_option / ((double)no_of_trials));
	continous_put_option = exp(-risk_free_rate * expiration_time)*(continous_put_option / ((double)no_of_trials));

	cout << "The average Call Price by explicit simulation  = " << continous_call_option << endl;
	cout << "The call price using the (1-p)-adjustment term = " << call_option_price << endl;
	cout << "Theoretical Call Price                         = " << closed_form_down_and_out_european_call_option() << endl;
	cout << "--------------------------------------" << endl;
	cout << endl;
	cout << "The average Put Price by explicit simulation  = " << continous_put_option << endl;
	cout << "The put price using the (1-p)-adjustment term = " << put_option_price << endl;
	cout << "Theoretical Put Price                         = " << closed_form_down_and_out_european_put_option() << endl;
	cout << "--------------------------------------" << endl;
	cout << endl;
	clock_t end = clock();
	cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	system("pause");
}