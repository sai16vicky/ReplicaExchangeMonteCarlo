#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Conformation.h"
#include "const.h"
#include "random.h"
#include "searchparameters.h"
#include "stopwatch.h"
#define sz(x) (int)x.size()

using namespace std;

typedef Conformation Replica;
typedef vector< Replica *> Replicas;
typedef vector< Replica *>::iterator ReplicaIterator;

class ReplicaExchange
{
public:
    Replicas v;
    searchparameters s;
    int no_of_replicas;
    double probchangevalues[1001][1001];
    double tempval[1001];
    
    ReplicaExchange(searchparameters &_s)
    {
        s = _s;
        double mintemp = _s.mintemp;
        double maxtemp = _s.maxtemp;
        no_of_replicas = _s.no_of_replicas;
        
        double cur = mintemp, step = 1.0;
        
        if(no_of_replicas > 1)
            step = (maxtemp - mintemp)/(no_of_replicas - 1);
        
        for(int i=0;i<no_of_replicas;i++)
            tempval[i] = cur, cur += step;
        
        for(int i=0;i<no_of_replicas;i++)
        {
            Replica *tmp;
            tmp = new Replica(_s.seq, tempval[i], _s.experimentID, i);
            
            if(sz(_s.initial) > 0)
            {
                tmp->setSequence(_s.initial);
                tmp->getFreeEnergy();
            }
        }
        
        
        for(int i=0;i<(no_of_replicas - 1);i++)
        {
            double x1 = 1.0/(K_b*tempval[i]), x2 = 1.0/(K_b2*tempval[i+1]);
            double diff = x2 - x1;
            
            for(int j=0;j<=abs(s.energy);j++)
            {
                double change = diff * j * (-1.0);
                probchangevalues[i][j] = exp(-change);
            }
        }
        
        Conformation::initProbAcceptValues(tempval, no_of_replicas, _s.energy*-1);
    }
    int conditionswap(int x,int y)
    {
        
        int ex = v[x]->getFreeEnergy(), ey = v[y]->getFreeEnergy(), diff = ex - ey;
        
        if(diff > 0)
            return 1;
        
        double rd = ran01();
        
        if(probchangevalues[x][-diff] > rd)
            return 1;
        
        return 0;
    }
    void swapreplicas(int x, int y)
    {
        
        Conformation *tmp = v[x];
        v[x] = v[y];
        v[y] = tmp;
        
        v[x]->setReplicaID(x);
        v[y]->setReplicaID(y);
    }
    void runreplicaexchange(vector<solution> &best)
    {
        int prevans, bestans = 0, fg = 0;
        string bdir;
        
        clock_t start = clock();
        
        while(bestans > s.energy and (double)(clock() - start)/(CLOCKS_PER_SEC) < s.runtime)
        {
            int tid = 0, tbest = 0, tenergy = 0, badval = 0, badenergy = -1000;
            
            for(ReplicaIterator i= v.begin();i != v.end();i++)
            {
                Replica *tmp = *i;
                searchresult &s1 = tmp->runSimulation(s.steps);
                solution &sol = s1.getminsol();
                
                if(tmp->getFreeEnergy() < tenergy)
                {
                    tenergy = tmp->getFreeEnergy();
                    tbest = tid;
                }
                
                if(tmp->getFreeEnergy() > badenergy)
                {
                    badenergy = tmp->getFreeEnergy();
                    badval = tid;
                }
                
                if(sol.energy < bestans)
                {
                    bestans = sol.energy;
                    sol.rt = (double)(clock() - start)/(CLOCKS_PER_SEC);
                }
                
                
                tid++;
            }
            
            int offset = fg;
            
            fg = (fg + 1)%2;
            
            int width = no_of_replicas - offset - 1;
            
            int st = offset;
            
            while(st < width + offset)
            {
                int next = st + 1;
                
                if(conditionswap(st, next))
                    swapreplicas(st, next);
                
                st += 2;
            }
            
            
            if(bestans != prevans)
            {
                cout<<(double)(clock() - start)/(CLOCKS_PER_SEC)<<" : Best Solution Now: "<<bestans<<"\n";
                
                prevans = bestans;
            }
            
            Conformation::swapOffset();
            
        }
        
        
        cout<<bdir<<"\n";
        
        Conformation *final = new Conformation(s.seq, 200, "", 1);
		
        final->loadConformation(bdir);
        
        final->printConformation();
    }
};

