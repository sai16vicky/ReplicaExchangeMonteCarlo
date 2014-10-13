#ifndef SEARCHPARAMETERS_H
#define SEARCHPARAMETERS_H
#include <string>
#include <vector>
#include "Const.h"

using namespace std;

class searchparameters
{
public:
    static int runs,no_of_replicas,steps,energy;
    static long long seedvalue;
    static double mintemp, maxtemp, runtime, pmweight;
    static string seq, experimentID, initial;
    static bool showwwbonds;
};

int searchparameters::no_of_replicas = 10, searchparameters::steps = 25, searchparameters::energy = 0, searchparameters::runs= 1;
long long searchparameters::seedvalue = 0;
double searchparameters::mintemp = 1.0,searchparameters::maxtemp = 225.0, searchparameters::pmweight = 0.4;
double searchparameters::runtime = 36000;
string searchparameters::seq = "wbwbwbwbwbwbwbw";
string searchparameters::experimentID = "default_experiment";
string searchparameters::initial = "";
bool searchparameters::showwwbonds = true;

#endif

