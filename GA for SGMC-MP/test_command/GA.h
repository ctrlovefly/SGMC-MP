#include<vector>
#include<stdlib.h>
#include<algorithm>
#include"chromo.h"
using namespace std;


struct MyPoint
{
	int chromostation;
	double swap_s;
};
class GA
{
public:
	GA()
	{
		pc=m_crossoverRate;
		pm=m_mutationRate;
	}
	~GA()
	{
	}
	int *already_add;
	double pc;
	double pm;
	void CreateGenomes(Test& TA);//generate a initial population
	void RankPopulation(Test& TA);//rank a population in ascending order
	void CreateNextGeneration(Test& GA,int generation);//generate the next population
	void select_match();//tournament selection
	Genome* GA::CrossoverOnePoint_Intvector_se(Genome genome1, Genome genome2, Test& GA,  Genome ta[]);//single-point crossover for the satellite swath segment
	Genome* GA::CrossoverUniform_Intvector(Genome genome1,Genome genome2, Test& GA, Genome ta[]);//uniform crossover for the ground station segment
	void GA::Mutate_normal(Genome &chromo,Test& TA);//mutation for the satellite swath segment	
	void Mutate_station(Genome &chromo,Test& TA);//mutation for the ground station segment

	template<typename T>
	T randT(T Lower, T Upper); 
	
	Genome m_wholeGeneration[m_populationSize*2];
	Genome m_thisGeneration[m_populationSize];
	Genome m_nextGeneration[m_populationSize];
	double m_preGeneration_fitness;
};
void GA::CreateGenomes(Test& TA)
{	
	already_add=new int[PIPS];
	int begin_num=0;
	for (int i = 0; i < m_populationSize; i++)
	{
		Genome gene; 
		gene.CreateGenesIntVector(TA);
		gene.m_fitness = gene.FitnessFunction(TA);
		m_thisGeneration[i]=gene;		                                  
	}
}
 bool comp(Genome x, Genome y)
 {
	 return ((Genome)x).m_fitness < ((Genome)y).m_fitness;
 }
 
 void GA::RankPopulation(Test& TA)
 {
	 sort(m_thisGeneration, m_thisGeneration + m_populationSize,comp);//ascending order
 }
 void GA::CreateNextGeneration(Test& TA,int generation)
 {
	 for (int i = 0; i < m_populationSize; i++)
		 m_wholeGeneration[i]=m_thisGeneration[i];
	 select_match();//tounament selection
	 for (int i = 0; i < m_populationSize; i = i + 2)
	 {			 
		 Genome parent1, parent2, child1, child2;
		 Genome child1_1,child2_1;
		 parent1 = m_thisGeneration[i];
		 parent2 = m_thisGeneration[i + 1];
		 //**************
		 //**crossover***
		 //**************
	    int flag=0;
		if ((double)(rand() % 101) / 101 < pc)
		{
			Genome te[2]={Genome(),Genome()};
			CrossoverUniform_Intvector(parent1, parent2, TA,  te);//uniform crossover
			parent1 = te[0];
			parent2 = te[1];
			CrossoverOnePoint_Intvector_se(parent1, parent2, TA,  te);//single-point crossover
			child1 = te[0];
			child2 = te[1];
		}
		else
		{
			child1 = parent1;
			child2 = parent2;
		}

		child1.m_fitness = child1.FitnessFunction(TA);
		child2.m_fitness = child2.FitnessFunction(TA);
		//**************
		//***mutation***
		//**************
		int mutate_time_1=1;
		for(int i=0;i<mutate_time_1;i++)
		{
			if ((double)(rand() % 101) / 101< pm)
			{
					Mutate_station(child1,TA);		
					Mutate_station(child2,TA);	 
			}
		}		 
		int mutate_time=sat;
		for(int i=0;i<mutate_time;i++)
		{
			if ((double)(rand() % 101) / 101< pm)
			{
				Mutate_normal(child1,TA);				
				Mutate_normal(child2,TA);			 
			}	
		}
		child1.m_fitness = child1.FitnessFunction(TA);
		child2.m_fitness = child2.FitnessFunction(TA);	 

		m_nextGeneration[i]=child1;
		m_nextGeneration[i+1]=child2;
	 }

	 //Replacement population method: steady state 
	 for (int i = 0; i < m_populationSize; i++)
		 m_wholeGeneration[i + m_populationSize] = m_nextGeneration[i];
	 sort(m_wholeGeneration, m_wholeGeneration + m_populationSize*2,comp);
	 
	 for (int i = 0; i < m_populationSize; i++)
		 m_thisGeneration[i]=m_wholeGeneration[i + m_populationSize];
 }

 void GA::select_match()		//use the match method
{
	Genome NewPopulation[m_populationSize ];
	int i;
	int p1,p2;
	for(i=0;i<m_populationSize;i++)
	{
		p1 = rand()%m_populationSize;
		p2 = rand()%m_populationSize;
		while(p1 == p2)
			p2 = rand()%m_populationSize;
		if(m_thisGeneration[p1].m_fitness < m_thisGeneration[p2].m_fitness)	//select the best one
			NewPopulation[i] = m_thisGeneration[p2];
		else
			NewPopulation[i] = m_thisGeneration[p1];
	}
	for (i = 0; i < m_populationSize; i++)
		m_thisGeneration[i] = NewPopulation[i];
}

Genome* GA::CrossoverOnePoint_Intvector_se(Genome genome1, Genome genome2, Test& GA,  Genome ta[])
{
	Genome child1;
	Genome child2;
	for(int i = 0;i < station; i++ )
	{
		child1.chromostation[i]=genome1.chromostation[i];
		child2.chromostation[i]=genome2.chromostation[i];
	}

	int *same_genes_1 = new int[sat];
	int *ps1_genes_1 = new int[sat];
	int *ps2_genes_1 = new int[sat];

	int **ps1_genes_roll= new int* [satcandidate];
	for (int i = 0; i < satcandidate; i++)
	{
		ps1_genes_roll[i] = new int[period];
	}
	int **ps2_genes_roll= new int* [satcandidate];
	for (int i = 0; i < satcandidate; i++)
	{
		ps2_genes_roll[i] = new int[period];
	}
	//extract the same genes
	int same_genes_roll[satcandidate][period];

	int same_1 = 0;

	for (int i = 0; i < sat; i++)
	{
		for (int j = 0; j <sat; j++)
		{
			if (genome1.chromosat[i].numsat == genome2.chromosat[j].numsat)
			{
				same_genes_1[same_1] = genome1.chromosat[i].numsat;
				for(int z=0;z<period;z++)
				{
					same_genes_roll[genome1.chromosat[i].numsat][z]=genome1.chromosat[i].numroll[z];
				}
				same_1++;
				break;
			}
		}
	}
	//extract the distinct genes for parent 1 
	int different = sat - same_1;
	int m = 0;
	int sameCount_1;
	for (int i = 0; i < sat; i++)
	{
		if (same_1 != 0)
		{
			sameCount_1 = 0;
			for (int j = 0; j < same_1; j++)
			{
				if (genome1.chromosat[i].numsat  == same_genes_1[j])
				{
					sameCount_1++;
				}
			}
			if (sameCount_1 == 0)
			{
				ps1_genes_1[m++] = (genome1.chromosat[i].numsat);
				for(int z=0;z<period;z++)
				{
					ps1_genes_roll[genome1.chromosat[i].numsat][z]=genome1.chromosat[i].numroll[z];
				}
			}
		}
		else
		{
	 		ps1_genes_1[m++] = genome1.chromosat[i].numsat;
			for(int z=0;z<period;z++)
			{
				ps1_genes_roll[genome1.chromosat[i].numsat][z]=genome1.chromosat[i].numroll[z];
			}
		}
	}
	//extract the distinct genes for parent 2 
	m = 0;
	for (int i = 0; i < sat; i++)
	{
		if (same_1 != 0)
		{
			sameCount_1 = 0;
			for (int j = 0; j < same_1; j++)
			{
				if (genome2.chromosat[i].numsat  == same_genes_1[j])
				{
					sameCount_1++;
				}
			}
			if (sameCount_1 == 0)
			{
				ps2_genes_1[m++] = (genome2.chromosat[i].numsat);
				for(int z=0;z<period;z++)
				{
					ps2_genes_roll[genome2.chromosat[i].numsat][z]=genome2.chromosat[i].numroll[z];
				}
			}
				
		}
		else
		{
			ps2_genes_1[m++] = genome2.chromosat[i].numsat;
			for(int z=0;z<period;z++)
			{
				ps2_genes_roll[genome2.chromosat[i].numsat][z]=genome2.chromosat[i].numroll[z];
			}
		}
			
	}
	int **cs1_genes_1= new int* [sat];
	for (int i = 0; i < sat; i++)
	{
		cs1_genes_1[i] = new int[4];
	}
	int **cs2_genes_1= new int* [sat];
	for (int i = 0; i < sat; i++)
	{
		cs2_genes_1[i] = new int[4];
	}

	//swap

	if (same_1 < sat)
	{

		int pos_1 = rand() % different-1;
		for (int i = 0; i < different; i++)
		{
			if (i < pos_1)
			{
				cs1_genes_1[i][0]=ps1_genes_1[i];
				cs2_genes_1[i][0]=ps2_genes_1[i];
				for(int j=0;j<period;j++)
				{
					cs1_genes_1[i][j+1]=ps1_genes_roll[ps1_genes_1[i]][j];
					cs2_genes_1[i][j+1]=ps2_genes_roll[ps2_genes_1[i]][j];
				}
			}
			else
			{
				cs1_genes_1[i][0] = ps2_genes_1[i];
				cs2_genes_1[i][0]= ps1_genes_1[i];
				for(int j=0;j<period;j++)
				{
					try
					{
					cs1_genes_1[i][j+1]=ps2_genes_roll[ps2_genes_1[i]][j];
					cs2_genes_1[i][j+1]=ps1_genes_roll[ps1_genes_1[i]][j];
					}
					catch(int e)
					{
						continue;
					}				
				}
			}
		}
		//get back the common genes
		for (int i = 0; i < same_1; i++)
		{
			cs1_genes_1[different + i][0] = same_genes_1[i];
			cs2_genes_1[different + i][0]= same_genes_1[i];
			for(int j=0;j<period;j++)
			{
				cs1_genes_1[different + i][j+1]=same_genes_roll[same_genes_1[i]][j];
				cs2_genes_1[different + i][j+1]=same_genes_roll[same_genes_1[i]][j];
			}
		}

		for (int i = 0; i < sat; i++)
		{
			child1.chromosat[i].numsat = cs1_genes_1[i][0];
			child2.chromosat[i].numsat = cs2_genes_1[i][0];
			for(int j=0;j<period;j++)
			{
				child1.chromosat[i].numroll[j]=cs1_genes_1[i][j+1];
				child2.chromosat[i].numroll[j]=cs2_genes_1[i][j+1];
			}
		}
	}
	else
	{
		for (int i = 0; i < sat; i++)
		{
			child1.chromosat[i].numsat = genome1.chromosat[i].numsat;
			child2.chromosat[i].numsat = genome2.chromosat[i].numsat;
			for(int j=0;j<period;j++)
			{
				child1.chromosat[i].numroll[j]=genome1.chromosat[i].numroll[j];
				child2.chromosat[i].numroll[j]=genome2.chromosat[i].numroll[j];
			}
		}
	}
	ta[0] = child1;
	ta[1] = child2;
	delete []cs1_genes_1;
	delete []cs2_genes_1;
	delete []ps2_genes_roll;
	delete []ps1_genes_roll;
	delete []ps2_genes_1;
	delete []ps1_genes_1;
	delete []same_genes_1;
	return ta;
}

void GA::Mutate_station(Genome &chromo,Test& TA)
{
	int index = rand() % (station);
	int candidates[PIPS];
	for (int i = 0; i < PIPS; i++)
	{
		candidates[i] = i;
	}

	for (int i = 0; i <station; i++)
	{
		candidates[chromo.chromostation[i]] = -1;	//delete the existing station index
	}
	int* genes=new int[PIPS - station];
	int k = 0;
	for (int i = 0; i < PIPS; i++)
	{
		if (candidates[i] != -1)
			genes[k++]=candidates[i];//avoid repetition
	}
	int site = rand() % (PIPS-station);
	chromo.chromostation[index] = genes[site];
	delete []genes;
}

Genome* GA::CrossoverUniform_Intvector(Genome genome1,Genome genome2, Test& GA, Genome ta[])
{
	Genome child1;
	Genome child2;

	//extract the same genes
	int *same_genes = new int[station];
	int *ps1_genes = new int[station];
	int *ps2_genes = new int[station];
	int same = 0;
	for (int i = 0; i < station; i++)
	{
		for (int j = 0; j < station; j++)
		{
			if (genome1.chromostation[i] == genome2.chromostation[j])
			{
				same_genes[same] = genome1.chromostation[i];
				same++;
				break;
			}
		}
	}
	int different = station - same;

	//extract the distinct genes for parent 1
	int m = 0;
	int sameCount;
	for (int i = 0; i < station; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome1.chromostation[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps1_genes[m++] = (genome1.chromostation[i]);
		}
		else
			ps1_genes[m++] = genome1.chromostation[i];
	}

	//extract the distinct genes for parent 2
	m = 0;
	for (int i = 0; i < station; i++)
	{
		if (same != 0)
		{
			sameCount = 0;
			for (int j = 0; j < same; j++)
			{
				if (genome2.chromostation[i] == same_genes[j])
				{
					sameCount++;
				}
			}
			if (sameCount == 0)
				ps2_genes[m++] = genome2.chromostation[i];
		}
		else
			ps2_genes[m++] = genome2.chromostation[i];
	}

	//swap

	int *cs1_genes=new int[station];
	int *cs2_genes=new int[station];
	if (same <station)
	{
		for (int i = 0; i < different; i++)
		{
			if (rand()%3== 1)
			{
				cs1_genes[i] = ps1_genes[i];
			}
			else
			{
				cs1_genes[i] = ps2_genes[i];
			}
		}
		for (int i = 0; i < different; i++)
		{
			if (rand() % 3== 1)
			{
				cs2_genes[i]=ps1_genes[i];
			}
			else
			{
				cs2_genes[i]=ps2_genes[i];
			}
		}

		//return the same genes
		for (int i = 0; i < same; i++)
		{
			cs1_genes[different + i] = same_genes[i];
			cs2_genes[different + i] = same_genes[i];
		}

		for (int i = 0; i < station; i++)
		{
			child1.chromostation[i] = cs1_genes[i];
			child2.chromostation[i] = cs2_genes[i];
		}
	}
	else
	{
		//if the two parents are the same, the swap is not conducted.
		for (int i = 0; i < station; i++)
		{
			child1.chromostation[i] = genome1.chromostation[i];
			child2.chromostation[i] = genome2.chromostation[i];
		}
	}

	//add the satellite swath segment
	for(int i = 0;i < sat; i++)
	{
		child1.chromosat[i]=genome1.chromosat[i];
		child2.chromosat[i]=genome2.chromosat[i];
	}
	ta[0] = child1;
	ta[1] = child2;
	delete []cs1_genes;
	delete []cs2_genes;
	delete []same_genes;
	delete []ps1_genes;
	delete []ps2_genes ;

	return ta;
}

void GA::Mutate_normal(Genome &chromo,Test& TA)
{
	//mutation for the satellite type index in a satellite swath segment 
	int index = rand() % (sat);
	int candidates_1[satcandidate];
	for (int i = 0; i < satcandidate; i++)
	{
		candidates_1[i] = i;
	}
	for (int i = 0; i <sat; i++)
	{
		candidates_1[chromo.chromosat[i].numsat] = -1;
	}
	int* genes_1=new int[satcandidate - sat];
	int k = 0;
	for (int i = 0; i < satcandidate; i++)
	{
		if (candidates_1[i] != -1)
			genes_1[k++]=candidates_1[i];
	}
	if((double)(rand() % 101) / 101<0.2)
	{
		if (satcandidate - sat>0)
		{
	int site = rand() % (satcandidate - sat);
	chromo.chromosat[index].numsat = genes_1[site];
		}
	}
	delete []genes_1;

	//mutation for the satellite swath index in a satellite swath segment 
	index = rand() % (sat);
	int period_index=rand() % (period);
	int candidates_2[swath];
	for (int i = 0; i < swath; i++)
	{
		candidates_2[i] = i;
	}
	candidates_2[chromo.chromosat[index].numroll[period_index]] = -1;

	if( TA.f[period_index][chromo.chromosat[index].numsat].size()!=0)
    {
		int n=rand()% TA.f[period_index][chromo.chromosat[index].numsat].size();
		chromo.chromosat[index].numroll[period_index] =TA.f[period_index][chromo.chromosat[index].numsat][n];			
	}
	else
	{
		while (1)
		{
			int site = rand() % (swath);
			if(candidates_2[site]!=-1)
			{
				chromo.chromosat[index].numroll[period_index] = candidates_2[site];
				break;
			}	
		}	
	}	
}