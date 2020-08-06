#include <string>
#include"all.h"
#include"Test.h"
using namespace std;

struct chromosat
{
	int numsat;
	int numroll[period];	
	double swap;

};//one column of the satellite swath segment

class Genome//one chromosome
{
public:
	Genome::Genome()
	{
	
	}
	Genome::~Genome()
	{
	}
	int chromostation[station];
	chromosat chromosat[sat];
	double m_fitness;
	double SumFitness;
	double FitnessFunction(Test& TA);
	void Genome::random(int a[], int n);
	void Genome::CreateGenesIntVector(Test& TA);
};

double Genome::FitnessFunction(Test& TA)
{
	double best_fitness_1;
	double best_fitness_2;
	double best_fitness_3;
	double ind_fitness = 0;
	int m;
	bool jump=true;

	for (int i = 0; i <= demandsize; i++)
	{
		best_fitness_1 = 0;
		best_fitness_2 = 0;
		best_fitness_3 = 0;
		for (int j = 0; j <station; j++)
		{
			if (TA.cover1[chromostation[j]][i] >0)
			{
				for (int q = 0; q < sat; q++)
				{
					if (TA.cover2[0][chromosat[q].numsat][chromosat[q].numroll[0]][i] >0)
					{
						best_fitness_1 =TA.cover2[0][chromosat[q].numsat][chromosat[q].numroll[0]][i];//Cover2
					}	
					if (TA.cover2[1][chromosat[q].numsat][chromosat[q].numroll[1]][i] >0)
					{
						best_fitness_2 =TA.cover2[1][chromosat[q].numsat][chromosat[q].numroll[1]][i];//Cover2
					}	
					if (TA.cover2[2][chromosat[q].numsat][chromosat[q].numroll[2]][i] >0)
					{
						best_fitness_3 =TA.cover2[2][chromosat[q].numsat][chromosat[q].numroll[2]][i];//Cover2
					}
				}
			}		
		}
		ind_fitness += best_fitness_1+best_fitness_2+best_fitness_3;
	}
	return ind_fitness;
}

void Genome::random(int a[], int n)
{
   int index, tmp, i;
   srand(time(NULL));
   for (i = 0; i <n; i++)
    {
       index = rand() % (n - i) + i;
       if (index != i)
         {
            tmp = a[i];
            a[i] = a[index];
            a[index] = tmp;
         }
    }
 }

void Genome::CreateGenesIntVector(Test& TA)
{
	/*******************************
	*****ground station segment*****
	********************************/	
	int num[PIPS];
	int r = PIPS - 1;
	int n;
	int tmp;
	for (int i = 0; i < PIPS; i++)
	{
		num[i] = i;
	}

	for (int i = 0; i < station; i++)
	{
		n = rand() % r + 1;   
		chromostation[i] = num[n];        
		tmp = num[n];               //avoid repetition
		num[n] = num[r];
		num[r] = tmp;
		r--;                       
	}
	/*******************************
	*****satellite swath segment****
	********************************/
	int num_1[satcandidate];
	int num_2[swath];
	int p = satcandidate-1;
	int q = swath-1;
	for (int i = 0; i < satcandidate; i++)
	{
		num_1[i] = i;
	}
	for(int i=0; i<swath;i++)
	{
		num_2[i]=i;
	}

	for (int i=0;i<sat;i++)
	{
		if (p==0)
		{
			n=0;
		}
		else
		{
			n=rand()%p+1;
		}
		chromosat[i].numsat=num_1[n];//satellite type index
		tmp = num_1[n];               
		num_1[n] = num_1[p];
		num_1[p] = tmp;
		p--; 
		for (int j = 0; j < period; j++)
		{
			if( TA.f[j][chromosat[i].numsat].size()!=0)
			{
				n=rand()% TA.f[j][chromosat[i].numsat].size();
				chromosat[i].numroll[j]=TA.f[j][chromosat[i].numsat][n];//satellite swath index
			}
			else
			{
				n=rand()% q+1;
				chromosat[i].numroll[j]=num_2[n];
				tmp = num_2[n];               
				num_2[n] = num_2[q];
				num_2[q] = tmp;
				q--; 
			}			
		}
	}	
}