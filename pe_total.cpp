#include "pe_i.h"
#include "pe_total.h"
#include <vector>
#include <iostream>

using namespace std;

double total_e (const double e, vector<vector<double> >& box_dim, const double s, const vector<vector<double> >& positions, const int N) {
    int i;
    double pe_t = 0.0;
    for (i=0; i<N-1; i++) {
        pe_t += e_i(i+1, e, box_dim, s, positions, i, N);
    }
    return pe_t;
}
