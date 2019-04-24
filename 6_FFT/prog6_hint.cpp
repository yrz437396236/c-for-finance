﻿// IE523: Financial Computation
// Processing tick-data using the Fast Fourier Transform (FFT)
// Written by Prof. R.S. Sreenivas

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmat.h"
#include "newmatap.h" // Need this for the FFT algorithm
#include"newmatio.h"

using namespace std;

// function that is used for sorting a 2-col array based on the value of the entries in the 2nd column.
// taken from http://stackoverflow.com/questions/3041897/sorting-a-2-dimensional-array-on-multiple-columns
bool compareTwoRows2(double* rowA, double* rowB){
	return ( (rowA[1]>rowB[1]) || ((rowA[1]==rowB[1])&&(rowA[0]>rowB[0])) );
}

class Filtering_Instance 
{ 
	// Private
	int no_of_terms, no_of_data_points;
	
	// Private member function that computes the mean of a data array
	double compute_mean(ColumnVector data)
	{
		double temp = 0;
		for (int i = 1; i <= data.size(); i++)
		{
			temp = data(i) + temp;
		}
		return(temp / data.size());			
	// write the code to compute the mean of "data"
	}
	
	// Private member function that computes the magnitude of (an array of) complex #s
	void compute_magnitude(ColumnVector &magnitude, ColumnVector real_part, ColumnVector imag_part)
	{
		for (int i = 1; i <= real_part.size(); i++)
		{
			magnitude(i)=sqrt(pow(real_part(i),2) + pow(imag_part(i), 2));
		}
	// write the code to compute sqrt(real_part(i)^2 + imag_part(i)^i)
	// is this imag_part(i)^i a mistake? 
	}	
	
	// Private member function that reads the data from the input file 
	// and stores it in "data" 
	void get_data(char* file_name, int no_of_data_points,ColumnVector &data)
	{
		//char* file_name
		int i = 1;
		double element = 0;
		ifstream input_file(file_name);
		if (!input_file)
		{
			cout << "cannot open input" << endl;
		}
		else
		{
			// read the data vector
			for (int i = 1; i <= no_of_data_points; i++)
			{
				input_file >> element;
				data(i) = element;
			}
			input_file.close();
		}
	// write code that reads the ticker-data from the input file
    // and stores it in "data"
	}
	
	// private member function that writes the data file into a file
	void write_data(char* file_name, ColumnVector &data)
	{
		//char* file_name
		ofstream output_file(file_name);
		if (!output_file)
		{
			cout << "cannot open output" << endl;
		}
		else
		{
			for (int i = 1; i <= data.size(); i++)
			{
			output_file << data(i) << endl;
			}
		}
		output_file.close();
	}
	// write code that writes "data" to file_name.

	
	// private member function that filters data using the FFT 
	// The filtered data is computed using the top "no_of_terms"-many 
	// magnitude-components of the orginal data
	void filter_the_data(ColumnVector &data, ColumnVector &filtered_data, int no_of_terms)
	{

		ColumnVector fft_real_part(data.nrows()), fft_imag_part(data.nrows());
		ColumnVector mean_adjusted_data(data.nrows()), magnitude(data.nrows());
		
		double mean = compute_mean(data);
		for (int i = 1; i <= data.nrows(); i++)
			mean_adjusted_data(i) = data(i) - mean;
		
		RealFFT(mean_adjusted_data, fft_real_part, fft_imag_part);
		compute_magnitude(magnitude, fft_real_part, fft_imag_part);
		
		// creating a two dimensional array: first col is the index; second col is the 
		// magnitude.  The plan is to have this 2-D array sorted using the 2nd col (ie.
		// sorted based on magnitude). Then we pick the top "no_of_terms"-many of these
		// components to reconstitute/reconstruct the signal back. 
		
		double** two_dimensional_array = new double*[magnitude.nrows()];
		for (int i = 0; i < magnitude.nrows(); i++)
			two_dimensional_array[i] = new double[2];
		
		for (int i = 0; i < fft_imag_part.nrows(); i++)
		{
			two_dimensional_array[i][0] = i;
			two_dimensional_array[i][1] = magnitude(i+1);
		}
		std::sort(two_dimensional_array,two_dimensional_array+fft_imag_part.nrows(),&compareTwoRows2);
		
		// if do_we_pick_this(i) == 1, then we keep that component for reconstruction
		// of the filtered signal.  The rest of the array-names should be self-explanatory
		ColumnVector do_we_pick_this(fft_imag_part.nrows());
		ColumnVector filtered_fft_real_part(fft_imag_part.nrows());
		ColumnVector filtered_fft_imag_part(fft_imag_part.nrows());
		for (int i = 0; i < no_of_terms; i++)
		{
			do_we_pick_this((two_dimensional_array[i][0]+1)) = 1;
		}
		for (int i = 1; i < fft_imag_part.nrows() + 1; i++)
		{
			if (do_we_pick_this(i) == 1)
			{
				filtered_fft_real_part(i) = fft_real_part(i);
				filtered_fft_imag_part(i) = fft_imag_part(i);
			}
			else
			{
				filtered_fft_real_part(i) = 0;
				filtered_fft_imag_part(i) = 0;
			}
		}
		// write the code for picking the top "no_of_terms" many magnitudes
		// put the real-part in "filtered_fft_real_part" and imaginary-part in 
		// "filtered_fft_imag_part" -- and reconstruct the filtered signal as 
		// shown below. 
		// reconstructed signal using just the "no_of_terms"-many top-magnitude 
		// components.
		//cout << setw(10) << fft_real_part << endl;
		//cout << setw(10) << filtered_fft_real_part << endl;
		RealFFTI(filtered_fft_real_part, filtered_fft_imag_part, filtered_data);	
		//cout << setw(10) << filtered_data << endl;
		// write code to add the mean-back to the reconstructed, filtered-signal
		for (int i = 1; i <= filtered_data.nrows(); i++)
		{
			filtered_data(i) = filtered_data(i) + mean;
		}
		//cout << setw(10) << filtered_data << endl;
	}
	
public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void run_the_filter(int argc, char * const argv[])
	{
		sscanf_s (argv[1], "%d", &no_of_terms);
		sscanf_s (argv[2], "%d", &no_of_data_points);
		//no_of_terms = 5;
		//no_of_data_points = 378;
		std::cout << "Input File Name: " << argv[3] << std::endl;
		std::cout << "Number of data points in the input file = " << no_of_data_points << endl;
		std::cout << "Number of dominant terms in the FFT = " << no_of_terms << endl;
		std::cout << "Output File Name: " << argv[4] << std::endl;
		
		ColumnVector data(no_of_data_points), filtered_data(no_of_data_points);
		
		// get ticker data
		get_data(argv[3], no_of_data_points,data);
		//get_data(no_of_data_points,data);
		//filter the ticker data
		filter_the_data(data, filtered_data, no_of_terms);
		
		// write the filtered data
		write_data(argv[4], filtered_data);
		//write_data(filtered_data);
	}
};

																	   
int main (int argc, char* argv[])
{
	Filtering_Instance x;
	x.run_the_filter(argc,argv);
}

