#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "newmatap.h"           
#include "newmat.h"
#include "newmatio.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>
#include <windows.h>

Matrix A;
double Random_Number()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(-5.0, 5.0);
	return(distribution(generator));
}

Matrix form_the_Matrix(int no_rows)
{	
	Matrix B(no_rows, no_rows);
	double temp = 0;
	for (int i = 1; i <= no_rows; i++)
	{
		for (int j = 1; j <= no_rows; j++)
		{
			temp = Random_Number();
			B(i,j) = temp;
		}
	}
	return(B);
}
void output(char* file_name, Matrix B, int no_rows)
{
	ofstream output_file(file_name);
	if (!output_file)
	{
		cout << "cannot open output" << endl;
	}
	else
	{
		for (int i = 1; i <= no_rows; i++)
		{
			for (int j = 1; j <= no_rows; j++)
			{
				output_file << B(i,j)<<" " ;
			}
			output_file << endl;
		}
	}
	output_file.close();
}

Matrix Matrix_repeated_squaring(Matrix A, int exponent, char* file_name, int no_rows)
{
	IdentityMatrix I(no_rows);

	if (exponent == 0) 
		return I;
	if (exponent == 1)
		return A;
	{
		if (exponent % 2 == 1) 
			// if exponent is odd
			return(A * Matrix_repeated_squaring(A*A, (exponent - 1) / 2, file_name, no_rows));
		else 
			//if exponent is even
			return (Matrix_repeated_squaring(A*A, exponent / 2, file_name, no_rows));
	}
}

void Brute_Force(Matrix A, int exponent, char* file_name, int no_rows)
{
	Matrix B = A;
	for (int i = 0; i < exponent - 1; i++)
	{
		B = A * B;
	}
	output(file_name, B, no_rows);
}
// with output for one exponent and can print the matrix

void Brute_Force2(Matrix A, int exponent, char* file_name, int no_rows)
{
	Matrix B = A;
	for (int i = 0; i < exponent - 1; i++)
	{
		B = A * B;
	}
}
// without return for count the time
int main(int argc, char * argv[])
{
	int exponent = 0;
	int no_rows = 0;
	sscanf_s(argv[1], "%d", &exponent);
	sscanf_s(argv[2], "%d", &no_rows);
	Matrix A(no_rows, no_rows);
	cout << "The number of rows/columns in the square matrix is: " << no_rows << endl;
	cout << "The exponent is: " << exponent << endl;
	cout << "Output File Name(Repeaded Squaring Algorithm): " << argv[3] << std::endl;
	cout << "Output File Name(Brute Force): " << argv[4] << std::endl;
	
	// form the Matrix
	A=form_the_Matrix(no_rows);

	clock_t start = clock();
	// Smart Way
	output(argv[3], Matrix_repeated_squaring(A, exponent, argv[3], no_rows), no_rows);
	clock_t end1 = clock();
	cout << "Repeated Squaring Result:" << endl;
	cout << "It took " << (long double)(end1 - start) / CLOCKS_PER_SEC <<" seconds to complete"<< endl;
	
	// Brute Force
	Brute_Force(A, exponent, argv[4], no_rows);
	clock_t end2 = clock();
	cout << "Direct Mutiplication Result:" << endl;
	cout << "It took " << (long double)(end2 - start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	ofstream output_file("Timer_k.txt");
	for (int i = 1; i <= 300; i++)
	{
		LARGE_INTEGER t1,t2,t3,t4,tc;
		QueryPerformanceFrequency(&tc);
		//clock_t start1 = clock();
		QueryPerformanceCounter(&t1);
		Matrix_repeated_squaring(A, i, argv[3], no_rows);
		QueryPerformanceCounter(&t2);
		//clock_t end1 = clock();
		//clock_t start2 = clock();
		QueryPerformanceCounter(&t3);
		Brute_Force2(A, i, argv[4], no_rows);
		//clock_t end2 = clock();
		QueryPerformanceCounter(&t4);
		output_file << scientific << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << " " << (t4.QuadPart - t3.QuadPart)*1.0 / tc.QuadPart << endl;
		//output_file<<scientific<<(long double)(end1 - start1) / CLOCKS_PER_SEC << " "<< (long double)(end2 - start2) / CLOCKS_PER_SEC<<endl;
		// If I use the clock_t function for the plot, I usually get 0.000 because the number is really small, but MAC may not face this problem.
	    // For my method:https://stackoverflow.com/questions/1739259/how-to-use-queryperformancecounter
	}
	ofstream output_file2("Timer_n.txt");
	for (int i = 2; i <= 300; i++)
	{
		Matrix B(i, i);
		B = form_the_Matrix(i);
		LARGE_INTEGER t1, t2, t3, t4, tc;
		QueryPerformanceFrequency(&tc);
		//clock_t start1 = clock();
		QueryPerformanceCounter(&t1);
		Matrix_repeated_squaring(B, exponent, argv[3], i);
		QueryPerformanceCounter(&t2);
		//clock_t end1 = clock();
		//clock_t start2 = clock();
		QueryPerformanceCounter(&t3);
		Brute_Force2(B, exponent, argv[4], i);
		//clock_t end2 = clock();
		QueryPerformanceCounter(&t4);
		output_file2 << scientific << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << " " << (t4.QuadPart - t3.QuadPart)*1.0 / tc.QuadPart << endl;
		//output_file<<scientific<<(long double)(end1 - start1) / CLOCKS_PER_SEC << " "<< (long double)(end2 - start2) / CLOCKS_PER_SEC<<endl;
		// If I use the clock_t function for the plot, I usually get 0.000 because the number is really small, but MAC may not face this problem.
		// For my method:https://stackoverflow.com/questions/1739259/how-to-use-queryperformancecounter
	}
	output_file.close();
	cout << "It's done" << endl;
	system("pause");
	return(0);
}
