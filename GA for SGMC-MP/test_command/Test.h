#include <SDKDDKVer.h>
#include <stdio.h>
#include <tchar.h>
#include <fstream>
#include <numeric>
#include "all.h"
using namespace std;
struct StationPoint
{
	double x;
	double y;	
	int index;
	double STrate;
};
struct SatelliteSwath
{
	int satNum;
	int rollNum;
	int periodNum;

};

class Test
{
public:
	Test();
	~Test();
	double sum[period][satcandidate][swath];
    StationPoint *PIPSpoint;
	double **cover1;
	double *DemandArea;
	double cover2[period][satcandidate][swath][demandsize];
	vector<int> split2(string str,char del);
	vector<int> f[3][15];
};
vector<int> Test:: split2(string str,char del) //£¨::192:168:ABC::416£©->£¨192 168 ABC 416£©
{
    stringstream ss(str);
    string tok;
    vector<int> ret;
    while(getline(ss,tok,del))
    {
        if(tok > "")
            ret.push_back(stoi(tok));
    }
    return ret;
};

Test::~Test()
{
}
Test::Test()
{	
	PIPSpoint = new StationPoint[PIPS];
	cover1 = new double*[PIPS];
	for (int i = 0; i < PIPS; i++)
	{
		cover1[i] = new double[demandsize];
	}
	DemandArea = new double[demandsize];

	//read the coverage relationship between candidate stations and demand objects
	ifstream Co1("PIPS.txt"); 
	for (int i = 0; i < PIPS; i++)
	{
		for (int j = 0; j < demandsize; j++)
		{
			Co1 >> cover1[i][j];
		}
	}

	//read the coverage relationship between strips and demand objects 
	ifstream Co1swath("period1.txt"); 
	for (int i = 0; i < satcandidate; i++)
	{
		for (int j = 0; j < swath; j++)
		{
			for(int z = 0;z < demandsize; z++)
			{
				Co1swath >> cover2[0][i][j][z];
			}
		}
	}

	ifstream Co1swath1("period2.txt"); 
	for (int i = 0; i < satcandidate; i++)
	{
		for (int j = 0; j < swath; j++)
		{
			for(int z = 0;z < demandsize; z++)
			{
				Co1swath1 >> cover2[1][i][j][z];
			}
		}
	}

	ifstream Co1swath2("period3.txt"); 
	for (int i = 0; i < satcandidate; i++)
	{
		for (int j = 0; j < swath; j++)
		{
			for(int z = 0;z < demandsize; z++)
			{
				Co1swath2 >> cover2[2][i][j][z];
			}
		}
	}

	for(int z=0;z<period;z++)
	{
		stringstream ss;
		ss<<z+1; 
		ifstream file(ss.str()+".txt");
		string str[15];
		int arraynum1=0;
		if(file.is_open())
		{
			for(int i=0;i<15;i++)
			{
				getline(file,str[arraynum1]);
				arraynum1++;
			}
		}
		//vector<int> f[15];
		int n=0;
		for(int j=0;j<15;j++)
		{
			if(!str[j].empty())
			{
				f[z][j]=(split2(str[j],','));
					n++;
			}
			else
			{
				
			}
		}			
	}

	for (int i = 0; i < period; i++)
	{
		for (int j = 0; j < satcandidate; j++)
		{
			for (int z = 0; z < swath; z++)
			{
				sum[i][j][z]=accumulate(cover2[i][j][z], cover2[i][j][z] + demandsize, 0);
			}
		}
	}
}