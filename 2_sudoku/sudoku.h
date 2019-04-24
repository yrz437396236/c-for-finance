/*
 *  sudoku.h
 *  Sudoku
 *  Created by Ruozhong Yang 
 *  Inspired by: http://web.eecs.utk.edu/courses/spring2012/cs140/Notes/Sudoku/index.html
*/

#ifndef sudoku
#define sudoku
#include <vector>
#include <fstream>
#include <iostream>
using std::vector;
using namespace std;
class Sudoku 
{ 
	// Private
	int puzzle[9][9];
    int m=0;
	
	// Private member function that checks if the named row is valid
	bool row_valid(int row)
	{
		int row1[10] = { 0,0,0,0,0,0,0,0,0,0 };
		for (int i = 0; i < 9; i++)
		{
            if (puzzle[row][i]!=0)
            {
                row1[(puzzle[row][i])] = row1[(puzzle[row][i])] + 1;
                if (row1[(puzzle[row][i])] > 1)
                {                    
                        return false;     
                }
                
            }
		}
		return true;
	}
		// write code that checks if "row" is valid
	// there is one more place in the row1, this place is for the "0" in the puzzle and will not be uesd

	
	// Private member function that checks if the named column is valid
	bool col_valid(int col)
	{
		int col1[10] = { 0,0,0,0,0,0,0,0,0,0 };
		for (int i = 0; i < 9; i++)
		{
            if (puzzle[i][col]!=0) {
                col1[(puzzle[i][col])] = col1[(puzzle[i][col])] + 1;
                if (col1[(puzzle[i][col])] > 1)
                {
                        return false; 
                }
            }
                
		}
		return true;
	}
		// check validity of "col" 
	//the same as row
	
	
	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col)
	{
		int a, b;
		a = row / 3;
		b = col / 3;
		int temp[9] = { puzzle[3 * a][3 * b], puzzle[3 * a + 1][3 * b], puzzle[3 * a + 2][3 * b], puzzle[3 * a][3 * b + 1], puzzle[3 * a + 1][3 * b + 1], puzzle[3 * a + 2][3 * b + 1], puzzle[3 * a][3 * b + 2], puzzle[3 * a + 1][3 * b + 2], puzzle[3 * a + 2][3 * b + 2] };
		int block1[10] = { 0,0,0,0,0,0,0,0,0,0 };
		for (int i = 0; i < 9; i++)
		{
			block1[(temp[i])] = block1[(temp[i])] + 1;
            if (block1[temp[i]] > 1)
            {
                if (temp[i] == 0)
                {     
                }else
                {
                    return false;   
                }
            }
		}
		return true;
	}	// check 3 x 3 block validity 
	//the same as row
	
	
public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void read_puzzle(int argc, char * const argv[])
	{
		ifstream infile(argv[1]);
		//ifstream infile;
		//infile.open("input4.txt");
		//this two lines is to read from inputs
        if(infile.fail())
        {
            cout<<"error"<<endl;
            system("pause");
        }
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				infile >> puzzle[i][j];
			}
		}
		infile.close();
        //the lines below is a function to check if it succeed
        //for (int i = 0; i < 9; i++)
		//{
			//for (int j = 0; j < 9; j++)
			//{
				//std::cout << puzzle1[i][j]<<endl;
			//}
		//}

	}
	// Public member function that prints the puzzle when called
	void print_puzzle()
	{
		std::cout << std::endl << "Board Position" << std::endl;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					std::cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					std::cout << "X ";
				}
			}
			std::cout << std::endl;
		}
	}
    void print_puzzle1()
    {
        std::cout << std::endl << "Final Solution" << std::endl;
        for (int i = 0; i < 9; i++)
        {
            for (int j = 0; j < 9; j++)
            {
                // check if we have a legitimate integer between 1 and 9
                if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
                {
                    // printing initial value of the puzzle with some formatting
                    std::cout << puzzle[i][j] << " ";
                }
                else {
                    // printing initial value of the puzzle with some formatting
                    std::cout << "X ";
                }
            }
            std::cout << std::endl;
        }
    }
	//change the cout of print_puzzle
	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	bool Solve(int row, int col)
    {
        for(int i=row;i < 9;i++)
        {
            for(int j=0;j < 9;j++)
            {
                if(puzzle[i][j]==0)
                {
                    for(int k=1;k<10;k++)
                    {
                        puzzle[i][j]=k;
                        if( block_valid(i, j)&&col_valid(j) &&row_valid(i) && Solve(i, j) )
                        { return true;}
                    }
                puzzle[i][j]=0;
                return false;
                }
            }
            
        }
        // this part of the code identifies the row and col number of the
		// first incomplete (i.e. 0) entry in the puzzle.  If the puzzle has
		// no zeros, the variable row will be 9 => the puzzle is done, as 
		// each entry is row-, col- and block-valid...
		
		// use the pseudo code of figure 3 of the description
		return true;
	}
    bool alternate_Solve(int row, int col)
    {
        int i=row,j=col;
        while(puzzle[i][j]!=0)
        {
            i++;
            if(i>8)
            {
                i=0;
                j++;
            }
            if(j>8){
                m++;
                cout<<endl<<"Solution#"<<m<<endl;
                print_puzzle1();
                return false;}
        }
        for(int k=1;k<10;k++)
    {
        puzzle[i][j]=k;
        if( block_valid(i, j)&&col_valid(j) &&row_valid(i) && alternate_Solve(i, j) )
        return true;
    }
    puzzle[i][j]=0;
    return false;
        }
};
// if we return false, the code will search again
#endif
