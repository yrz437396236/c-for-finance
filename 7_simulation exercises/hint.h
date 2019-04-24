/*
 *  alice_and_bob.h
 *  Loosing as little as possible
 *
 *  Created by Ramavarapu Sreenivas on 9/2/12.
 *  Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved.
 *
 */
#ifndef ALICE_AND_BOB
#define ALICE_AND_BOB

#include <cmath>
#include <chrono>
#include <random>

using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
	double alice_probability, bob_probability;
	
	// private member function: uniform RV generator
	double get_uniform()
	{
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		return(distribution(generator));
		// write the appropriate code here
	}
	
	// private member function: nCi (i.e. n-take-i) 
	int take(int n, int i)
	{
			if (i <= 1)
			{
				if(i==1)
				{
					return(n);
				}
				else
				{
					return(1);
				}
			}
			else 
			{
				return(((double)(n - i + 1) / i)*take(n, i - 1));
			}		
		// write a **RECURSIVE** implementation of n-take-i. 
		// If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it 
		// will take too long for large sizes 
	}
	
	// this routine implements the probability that Alice has more 
	// heads than Bob after n-many coin tosses
	double theoretical_value(double q, double p, int n)
	{
		double temp2 = 0;
		for (int r = 0 ; r < n; r++)
		{
			double temp = 0;
			for (int s = r+1; s <= n; s++)
			{
				temp = take(n, s)*pow(q, s)*pow((1 - q), (n - s))+temp;
			}
			temp2= take(n, r)*pow(p, r)*pow((1 - p), (n - r))*temp + temp2;
			//cout << take(n, r) <<" "<<pow(p, r) <<" "<< pow((1 - p), (n - r)) <<  endl;
		}
		return(temp2);
		// implement equation 1.1 of Addona-Wagon-Wilf paper
	}

public: 
	// public function: 
	void set_probability(double alice_p, double bob_p)
	{
		alice_probability = alice_p;
		bob_probability = bob_p;
	}
	
	// probability of Alice winning the game.
	double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
	{
		int no_of_wins_for_alice = 0;
		for (int i = 0; i < no_of_trials; i++) 
		{
			int number_of_heads_for_alice = 0;
			int number_of_heads_for_bob = 0;
			for (int j = 0; j < number_of_coin_tosses_in_each_game; j++) 
			{
				if (get_uniform() < alice_probability) 
					number_of_heads_for_alice++;
				if (get_uniform() < bob_probability)
					number_of_heads_for_bob++;
			}
			if (number_of_heads_for_alice > number_of_heads_for_bob)
				no_of_wins_for_alice++;
		}
		return (((double) no_of_wins_for_alice)/((double) no_of_trials));
	}
		
	int search_result()
	{
		for (int i = 1; i < 100; i++)
		{
			double temp=0, temp1=0, temp2 = 0;
			temp = theoretical_value(alice_probability, bob_probability, i-1);
			temp1 = theoretical_value(alice_probability, bob_probability, i);
			temp2 = theoretical_value(alice_probability, bob_probability, i+1);
			//cout << temp1 << endl;
			if ((temp1 >= temp )& (temp1 >= temp2))
			{
				return(i);
				break;
			}
		}
		return(0);
		// implememt a discrete-search procedure for the optimal n-value. 
		// start with n = 1 and find the discrete-value of n that has 
		// the largest probability for Alice winning.  Why would this work?
		// See Theorem 2.2 of the paper for the reason!
	}
	void print_data(int no_of_trials)
	{
		for (int i = 1; i < 31; i++)
		{
			cout << theoretical_value(alice_probability, bob_probability, i) <<" "<< simulated_value(i, no_of_trials) << endl;
		}
	}

};
#endif









