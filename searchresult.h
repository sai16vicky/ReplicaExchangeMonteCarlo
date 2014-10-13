

#ifndef SEARCH_RESULT_H
#define SEARCH_RESULT_H

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include "const.h"

using namespace std;


struct solution
{
    long double rlen;
    double rt, t;
    int energy;
    string dir;
};

class searchresult
{
public:
    
    solution sol;
    string id;
    int sampleid;
    long double ma, mr, mi, steps, nr;
    //ostream sol;
    //vector< sol > trace;
    
    searchresult(const string &_id, const int &_sampleid)
    {
        id = _id, sampleid = _sampleid;
        ma = 0, mr = 0, mi = 0, sol.rt = 0, sol.rlen = 0, sol.energy = 0, sol.dir = "";
    }
    void write()
    {
        return ;
    }
    void addsol(const long double rlen, const double rt, const int energy, const double t, const int rdir[], const int dlen)
    {
        if(energy < sol.energy)
        {
            sol.rlen = rlen;
            sol.rt = rt;
            sol.energy = energy;
            sol.t = t;
            
            ostringstream ss;
            
            for(int i=0;i<dlen;i++)
            {
                if(rdir[i] == 0)ss << "s";
                else if(rdir[i] == 1)ss<<"l";
                else if(rdir[i] == 2)ss<<"r";
                else if(rdir[i] == 3)ss<<"u";
                else if(rdir[i] == 4)ss<<"d";
                else if(rdir[i] == 5)ss<<"g";
            }
            
            sol.dir = ss.str();
            
        }
    }
    solution& getminsol(){return sol;}
    void setid(const string &_id){id = _id;}
    void setsampleid(const int _sampleid){sampleid = _sampleid;}
    string &getid(){return id;}
    int getsampleid(){return sampleid;}
    void setma(long double tmp){ma = tmp;}
    void setmr(long double tmp){mr = tmp;}
    void setmi(long double tmp){mi = tmp;}
    void setsteps(long double tmp){steps = tmp;}
    void setresidue(long double tmp){nr = tmp;}
    void resetminsol(){sol.rlen = 0, sol.rt = 0, sol.energy = 0;}
    void setminsol(const int rdir[], const int _energy, int _dlen)
    {
        sol.rlen = 0, sol.rt = 0, sol.energy = _energy, sol.t = 0;
        ostringstream ss;
        
        for(int i=0;i<_dlen;i++)
        {
            if(rdir[i] == 0)ss << "s";
            else if(rdir[i] == 1)ss<<"l";
            else if(rdir[i] == 2)ss<<"r";
            else if(rdir[i] == 3)ss<<"u";
            else if(rdir[i] == 4)ss<<"d";
            else if(rdir[i] == 5)ss<<"g";
        }
        
        sol.dir = ss.str();
    }
    long double getma(){return ma;}
    long double getmr(){return mr;}
    long double getmi(){return mi;}
    long double getsteps(){return steps;}
    long double getnr(){return nr;}
    // vector<sol> getsol();
    
};

#endif
