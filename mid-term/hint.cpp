// Written by Prof. Sreenivas for IE523: Financial Computing

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>
#include "lp_lib.h"


using namespace std;

const double error = 1e-10;
int number_of_cash_flows;
vector <double> price_list;
vector <double> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;
vector <vector <double> > database;
double PV_of_obligation;
double average_value_of_the_YTMs;


double function(vector <double> cash_flow, double price, int maturity, double rate)
{
	double temp1 = 0; double temp2 = 0;
	for (int i = 0; i < cash_flow.size(); i++)
	{
		temp1 = cash_flow[i] * pow((1 + rate), (maturity - i - 1)) + temp1;
	}
	temp2 = (price*pow((1 + rate), maturity)) - temp1;
	return temp2;
	// write a function that computes f(r) in page 2 of lesson 3 of my notes
	// temp1 used to add all the Pt(1 + r)^(M-t) together and temp2 is the result of f(r)

}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
	double temp1 = 0; double temp2 = 0;
	for (int i = 0; i < (cash_flow.size() - 1); i++)
	{
		temp1 = (maturity - i - 1)* (cash_flow[i]) * pow((1 + rate), (maturity - i - 2)) + temp1;
	}
	temp2 = (maturity*price*pow((1 + rate), (maturity - 1))) - temp1;
	return temp2;
	// write a function that computes f'(r) in the bottom of page 2 of lesson 3 
	// of my notes
	// the same as f(r) part, temp1 add together and temp2 get the result
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
	while (function(cash_flow, price, maturity, rate) >= error)
	{
		rate = rate - (function(cash_flow, price, maturity, rate) / derivative_function(cash_flow, price, maturity, rate));
	}
	return rate;
	// write a function that finds the (only) +ve root of f(r) of page 2 of 
	// lesson 3 using Newton-Raphson method
	// call the two funcion above and use while loop to calculate the rate of each cash flow when f(rate)<error 
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
	double temp = 0;
	for (int i = 0; i < cash_flow.size(); i++)
	{
		temp = (cash_flow[i] * (i + 1)) / pow((1 + rate), (i + 1)) + temp;
	}
	temp = temp / price;
	return temp;
	// write a function that computes the duration of a cash flow
	// add together then divide by price from price list
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
	double temp = 0;
	for (int i = 0; i < cash_flow.size(); i++)
	{
		temp = (cash_flow[i] * (i + 2)*(i + 1)) / pow((1 + rate), (i + 3)) + temp;
	}
	temp = temp / price;
	return temp;
	// write a function that computes the convexity of a cash flow
	// add together then divide by price from price list
}

double present_value_of_debt(double average_YTMs, double amount, double year)
{
	double PV_of_debt = 0;
	PV_of_debt = amount / (pow((1 + average_YTMs), year));
	return PV_of_debt;
	// compute PV of future debt obligation
	// using the average-value-of-the-YTMs 
}

void print_data(char *filename)
{
	cout << "Input File: " << filename << endl;
	cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
	cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		cout << "---------------------------" << endl;
		cout << "Cash Flow #" << i + 1 << endl;
		cout << "Price = " << price_list[i] << endl;
		cout << "Maturity = " << maturity_list[i] << endl;
		cout << "Percentage of Face Value that would meet the obligation = " <<
			percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
		cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
		cout << "Duration = " << duration[i] << endl;
		cout << "Duration {to be used in LP_formulation below} = " <<
			(duration[i]/percentage_of_cash_flow_to_meet_debt_obligation[i]) << endl;
		cout << "{Note} " << duration[i] <<" = "<< 
			(duration[i] / percentage_of_cash_flow_to_meet_debt_obligation[i])<<
			" X "<< percentage_of_cash_flow_to_meet_debt_obligation[i]<< endl;
		cout << "Convexity = " << convexity[i] << endl;
		cout << "Convexity {to be used in LP_formulation below} = " <<
			(convexity[i] / percentage_of_cash_flow_to_meet_debt_obligation[i]) << endl;
		cout << "{Note} " << convexity[i] << " = " <<
			(convexity[i] / percentage_of_cash_flow_to_meet_debt_obligation[i]) <<
			" X " << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
	}
	cout << "***************************" << endl;
	cout << "Average YTM{which I use to computr PV of Debt} = "<< average_value_of_the_YTMs <<endl;
	cout << "Present value of debt = " << PV_of_obligation << endl;
	cout << "***************************" << endl;
}

void get_data(char* argv[])
{
	// the part use the logic like this:
	// first, I read the number of cash flows and use this number to form a vector of vector called database
	// database is with (2*number of cash flows+1)rows and later I will put what I read into this huge vector of vector
	// take the input1 as an example
	// row numbers     context
	// 0               the price and maturity of Cash_Flow1: 1131.27 10
	// 1               the cash flow: 67 67 67 67 67 67 67 67 67 1067
	// 2               the price and maturity of Cash_Flow2: 1069.88 15
	// 3               the cash flow: 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 69.88 1069.88 
	// 4               the price and maturity of Cash_Flow3: 863.5 30
	// 5               the cash flow: 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 59 1059 
	// 6               the price and maturity of Cash_Flow4: 1148.75 12
	// 7               the cash flow: 75 75 75 75 75 75 75 75 75 75 75 1075
	// 8               the price and maturity of Cash_Flow5: 1121.39 11
	// 9               the cash flow: 70 70 70 70 70 70 70 70 70 70 1070
	// 10              the price of obligation and the time: 1790.85 10
	// then I put each part into the vector or double where they actually need to go
	// by using many for loop

	double temp = 0;
	//ifstream input_filename;
	//input_filename.open("input1");
	ifstream input_filename(argv[1]);
	if (input_filename.is_open()) 
	{
	    input_filename >> temp;
	    number_of_cash_flows = temp;
		for (int j = 0; j < number_of_cash_flows; j++)
		{
			vector<double> cf_part1;
			for (int i = 0; i < 2; i++)
			{
				input_filename >> temp;
				cf_part1.push_back(temp);
			}
			database.push_back(cf_part1);
			vector<double> cf_part2;
			for (int i = 0; i < cf_part1[1]; i++)
			{
				input_filename >> temp;
				cf_part2.push_back(temp);
			}
			database.push_back(cf_part2);
		}
		vector<double> debt_obligation;
		for (int i = 0; i < 2; i++)
		{
			input_filename >> temp;
			debt_obligation.push_back(temp);
		}
		database.push_back(debt_obligation);
		debt_obligation_amount = database[2 * number_of_cash_flows][0];
		time_when_debt_is_due = database[2 * number_of_cash_flows][1];
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			price_list.push_back(database[2 * i][0]);
			maturity_list.push_back(database[2 * i][1]);
		}
		vector <vector <double> > all_cash_flow;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			vector <double> cash_flows_temp;
			for (int j = 0; j < (database[2 * i + 1].size()); j++)
			{
				cash_flows_temp.push_back(database[2 * i + 1][j]);
			}
			all_cash_flow.push_back(cash_flows_temp);
		}

		for (int i = 0; i < number_of_cash_flows; i++)
		{
			yield_to_maturity.push_back(Newton_Raphson(all_cash_flow[i], price_list[i], maturity_list[i], 0.1));
		}
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			duration.push_back(get_duration(all_cash_flow[i], price_list[i], maturity_list[i], yield_to_maturity[i]));
		}
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			convexity.push_back(get_convexity(all_cash_flow[i], price_list[i], maturity_list[i], yield_to_maturity[i]));
		}
		double sum = std::accumulate(std::begin(yield_to_maturity), std::end(yield_to_maturity), 0.0);
		average_value_of_the_YTMs = sum / yield_to_maturity.size();
		PV_of_obligation = debt_obligation_amount / (pow(1 + average_value_of_the_YTMs, time_when_debt_is_due));
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			percentage_of_cash_flow_to_meet_debt_obligation.push_back(PV_of_obligation / price_list[i]);
		}
		
	}
	else {
		cout << "Input file missing" << endl;
		exit(0);
	}
	// write the code that reads the data from the file identified 
	// on the command-line. 	
}

void get_optimal_portfolio()
{
	// I discuss the lpsolve part with my classmate Rain
	// in this part I add three kinds of constraints
	// 1. x1,x2.......xn > 0
	// 2. meet the need of duration
	// 3. x1+x2+.....+xn = 1
	// then I us lpsolve to the the max of the combination of convexity with this constraints
	lprec *lp;
	vector<double> lambda;
	double* row = NULL;
	double* solution = NULL;
	row = new double[number_of_cash_flows + 1];
	solution = new double[number_of_cash_flows];
	lp = make_lp(0, number_of_cash_flows);
	set_verbose(lp, 3);
	
	for (int i = 1; i < number_of_cash_flows + 1; i++)
	{
		for (int j = 0; j < number_of_cash_flows + 1; j++)
		{
			row[j] = 0;
		}
		row[i] = 1;
		add_constraint(lp, row, GE, 0);
	}

	for (int i = 0; i < number_of_cash_flows; i++)
	{
		row[i + 1] = duration[i];
	}
	add_constraint(lp, row, EQ, time_when_debt_is_due);

	for (int i = 0; i < number_of_cash_flows; i++)
	{
		row[i + 1] = 1;
	}
	add_constraint(lp, row, EQ, 1);

	for (int i = 0; i < number_of_cash_flows; i++)
	{
		row[i + 1] = convexity[i];
	}
	set_obj_fn(lp, row);
	set_maxim(lp);
	if (solve(lp) == 0)
	{
		get_variables(lp, solution);
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			lambda.push_back(solution[i]);
		}
		double cov = 0;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			cov = cov + lambda[i] * convexity[i];
		}
		cout << "Largest Convexity we can get is: " << cov << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			cout << "%Cash Flow:" << (i + 1) << "  " << ((lambda[i] * PV_of_obligation) / price_list[i]) << endl;
		}
		cout << "That is, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			if (lambda[i] > 0)
			{
				cout << "$" << (lambda[i] * PV_of_obligation) << " of Cash Flow#" << i+1 << endl;
			}
		}
	}
	else
	{
		cout << "There is no portfolio that meets the duration constraint of " << time_when_debt_is_due << "years" << endl;
	}
	delete_lp(lp);
}
// write the lp_solve specific material that 
// computes the optimal_portfolio



int main(int argc, char* argv[])
{
	if (argc == 1) {
		cout << "Input filename missing" << endl;
	}
	else 
	{
	get_data(argv);

	print_data(argv[1]);

	get_optimal_portfolio();
	}
	system("pause");
	return (0);
}

