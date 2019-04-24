#include <iostream>
#include "card_game.h"
#include <ctime>
#include <vector>

using namespace std;
	
int main (int argc, char * const argv[]) 
{
	int total_number_of_cards;
	sscanf (argv[1], "%d", &total_number_of_cards);
	cout << "Total Number of Cards = " << total_number_of_cards << endl;
	clock_t start = clock();
	double **array = new double*[total_number_of_cards/2];
	for (int i = 0; i < total_number_of_cards / 2; i++)
		array[i] = new double[total_number_of_cards / 2];
	for (int i = 0; i < total_number_of_cards / 2; i++)
		for (int j = 0; j < total_number_of_cards / 2; j++)
		{
			array[i][j] = -1;
		}
	cout << "Value of the game = " << value(total_number_of_cards/2,total_number_of_cards/2, array) << endl;
	clock_t end = clock();
	cout << "time(s):" << (double)(end - start) / CLOCKS_PER_SEC << endl;
	system("pause");
	return 0;
}
