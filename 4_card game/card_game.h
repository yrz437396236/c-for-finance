#ifndef	CARDGAME_H
#define CARDGAME_H
#include <algorithm>
#include <vector>
using namespace std;
double value(int r, int b, double **memo)
{
	if (0 == r)
		return ((double) b);
	if (0 == b)
		return (0); 
	if (memo[r-1][b-1] !=-1) 
	{
		return memo[r - 1][b - 1];
    }
	else
	{
		double term1 = ((double)r / (r + b)) * value(r - 1, b, memo);

		double term2 = ((double)b / (r + b)) * value(r, b - 1, memo);

		memo[r - 1][b - 1] = std::max((term1 + term2), (double)(b - r));

		return memo[r - 1][b - 1];
		//use a matrix to remember the solution,if there is an solution, read it, dont caculate it again, this will save time
	}

}
#endif