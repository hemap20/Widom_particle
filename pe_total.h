#ifndef PE_TOTAL_H
#define PE_TOTAL_H

#include "pe_i.h"
#include "pe_total.h"
#include <vector>
#include <iostream>

using namespace std;

double total_e (vector<vector<double> >& box_dim,const vector<vector<double> >& positions, const int N, double * vir);

#endif