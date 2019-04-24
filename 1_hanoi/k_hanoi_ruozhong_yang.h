// k_hanoi.h
// Written by Ruozhong Yang for IE523

#ifndef K_HANOI_H
#define K_HANOI_H

#include <cmath>
#include <deque>
using std::deque;

class Towers 
{
	// private data member: #steps 
	int number_of_steps;
	
	// private data member: #pegs
	int number_of_pegs;
	
	// private data member: #disks on 1st Peg
	int number_of_DISKS;
	
	// private data member: k-Deques-of-Deques
	// towers[0] is deque for the 1st Peg, towers[1] is the deque 
	// for the 2nd Peg,..., towers[k-1] is the deque for the kth Peg
	deque < deque <int> > towers; 
	
	// private member function:  returns 'false' if a larger disk
	// is moved on top of a smaller disk during a move
	bool is_everything_legal() 
	{
		bool result=true;
		for(int i=0;i<number_of_pegs;i++)
		    for(int j=0;j<towers[i].size();j++)
		         for(int k=j;k<towers[i].size();k++)
					 if (towers[i][k] < towers[i][j]) {
						 result = false;
					 }
		return (result);
		// go through each peg and each disks to find if there is a big on a small one 
	}
	
	// private member function: initializes the pegs and #disks on Peg 1
	void initialize(int m, int n)
	{
		number_of_steps = 0;
		number_of_DISKS = n;
		number_of_pegs = m;
		for (int i = 0; i < m; i++)
		{
			deque <int> peg;
			towers.push_back(peg);
		}
		for (int i = 0; i < n; i++)
			towers[0].push_back(i + 1);
		// form the matrix of deque of deque and put disks on the first peg 
	}
	// private member function:  picks the #disks to be moved based on
	// the #free pegs that are available at a given instant
	int pick_the_right_number_to_move(int number_of_disks, int number_of_free_pegs)
	{
		if (1 == number_of_free_pegs)
			return (number_of_disks-1);
		else 
			return (number_of_disks/2);
	}
	
	// private member function: primitive/basic move of the top disk
	// from the source peg to the destination peg -- and legality-check
	void move_top_disk(int source, int dest)
	{		
		std::cout << "Move disk " << towers[source].front() << " from Peg ";
		std::cout << source + 1 << " to Peg " << dest + 1;
		towers[dest].push_front(towers[source].front());
		towers[source].pop_front();
		if (true == is_everything_legal())
			std::cout << " (Legal)" << std::endl;
		else
			std::cout << " (Illegal)" << std::endl;
	}
	
	// private member function: recursive solution to the m Peg Tower of Hanoi
	void Move(int number_of_disks, int source_peg_number, int destination_peg_number, deque <int> free_peg_numbers)
	{
		if (free_peg_numbers.size() != 0)
		{
			if (number_of_disks > 1)
			{
				
					int m = pick_the_right_number_to_move(number_of_disks, free_peg_numbers.size());
					int a = free_peg_numbers.at(free_peg_numbers.size() - 1);
					free_peg_numbers.pop_back();
					free_peg_numbers.push_back(destination_peg_number);
					Move(m, source_peg_number, a, free_peg_numbers);
					free_peg_numbers.pop_back();
					Move(number_of_disks - m, source_peg_number, destination_peg_number, free_peg_numbers);
					free_peg_numbers.push_back(source_peg_number);
					Move(m, a, destination_peg_number, free_peg_numbers);
				
			}
			
				// basicially according to the process from hints, the first thing is to choose a destination for the first move, then move the m disks to that peg 
			    //later put the rest disks on the first peg to the destination, then put the disks on the a peg to the distination   
		
			else
			{
				// number_of_disks == 1
				number_of_steps++;
				move_top_disk(source_peg_number, destination_peg_number);
			}
		}
		else 
		{
			
			// there are no free pegs
			number_of_steps++;
			if(number_of_disks=1)
			move_top_disk(source_peg_number, destination_peg_number);
			else {

			}
		}
	}
	
	// private member function: prints the state of Peg m (m = 0, 1, 2)
	void print_peg_state(int m)
	{
		std::cout << "-----------------------------" << std::endl;
		std::cout << "State of Peg " << m + 1 << " (Top to Bottom): ";
		for (int i = 0; i < towers[m].size(); i++)
			std::cout << towers[m][i] << " ";
		std::cout << std::endl;
		std::cout << "Number of Steps = " << number_of_steps << std::endl;
		std::cout << "-----------------------------" << std::endl;
		// the same as three pegs
	}
public:
	// Solves the K Peg Tower of Hanoi
	void solve(int number_of_pegs, int number_of_disks)
	{
		initialize(number_of_pegs, number_of_disks);
		print_peg_state(0);
		
		deque <int> free_peg_numbers;
		// source peg # is 0; destination peg # is no_of_pegs-1; 
		// free pegs {1,.., no_of_pegs-2}
		for (int i = 1; i <= number_of_pegs-2; i++)
			free_peg_numbers.push_back(i);
		Move(number_of_disks, 0, number_of_pegs-1, free_peg_numbers);
		print_peg_state(number_of_pegs-1);
	}
};
#endif