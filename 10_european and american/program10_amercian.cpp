#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price, downtick_prob, notick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

float max(float a, float b) {
	return (b < a) ? a : b;

}

float american_call_option(int k, int i, float current_stock_price, double **memo) {
	if (memo[k][i + k] != -1)
	{
		return memo[k][i + k];
	}
	else
	{
		if (k == no_of_divisions)
			return max(0.0, (current_stock_price - strike_price));
		else
			memo[k][i + k] = max((current_stock_price - strike_price),
							 (uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor,memo) +
				             notick_prob*american_call_option(k + 1, i, current_stock_price,memo) +
				             downtick_prob*american_call_option(k + 1, i - 1, current_stock_price / up_factor,memo)) / R);
		    return memo[k][i + k];
	}
}

float american_put_option(int k, int i, float current_stock_price, double **memo) {
	if (memo[k][i + k] != -1)
	{
		return memo[k][i + k];
	}
	else
	{
		if (k == no_of_divisions)
			return max(0.0, (strike_price - current_stock_price));
		else
			memo[k][i + k] = max((strike_price - current_stock_price),
							 (uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor,memo) +
							 notick_prob*american_put_option(k + 1, i, current_stock_price,memo) +
							 downtick_prob*american_put_option(k + 1, i - 1, current_stock_price / up_factor,memo)) / R);
			return memo[k][i + k];
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
	double **array1 = new double*[no_of_divisions + 1];
	for (int i = 0; i <= no_of_divisions; i++)
		array1[i] = new double[2 * no_of_divisions];
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

	up_factor = exp(volatility*sqrt((2 * expiration_time) / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - (1 / (sqrt(up_factor)))) / ((sqrt(up_factor)) - (1 / (sqrt(up_factor)))), 2);
	downtick_prob = pow((sqrt(up_factor) - (sqrt(R))) / ((sqrt(up_factor)) - (1 / (sqrt(up_factor)))), 2);
	notick_prob = 1 - uptick_prob - downtick_prob;
	cout << "Recursive Binomial American-Asian Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = american_call_option(0, 0, initial_stock_price, array1);
	cout << "Trinomial Price of an American Call Option = " << call_price << endl;
	double put_price = american_put_option(0, 0, initial_stock_price, array2);
	cout << "Trinomial Price of an American Put Option = " << put_price << endl;
	clock_t end = clock();
	cout << "It took " << (long double)(end - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	system("pause");

}