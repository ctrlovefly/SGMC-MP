#include <algorithm>
#include <iostream>
#include <time.h> 
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <sstream>
#include <cstddef>
#include "GA.h"
#include <math.h>

using namespace std;

string getTime()
{
	time_t timep;
	time(&timep);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
	return tmp;
}

int main() {
	Test ta;

	for(int i=0;i<10;i=i+1)//10 trials
	{
		srand((unsigned)time(NULL));
		int fitcount = 0;     
		GA ga;
		cout << getTime()<<" ";
		ga.CreateGenomes(ta);//generate a initial chromosome
		ga.RankPopulation(ta); 
		ga.m_preGeneration_fitness=ga.m_thisGeneration[m_populationSize - 1].m_fitness;
		for (int i = 0; i < m_generationLimit; i++)
		{		
			ga.CreateNextGeneration(ta,i);   
			cout << i << " ";
			cout << ga.m_thisGeneration[m_populationSize - 1].m_fitness<< " ";
		}
		cout << getTime() << endl;
		int stationbool[PIPS];
		int satbool[satcandidate];
		cout << ga.m_thisGeneration[m_populationSize - 1].m_fitness<<endl;		
	}
	system("pause");	
}




