
#include"SimplexMethod_.h"
#include<iostream>
#include<vector>
using namespace std;

int M = 3;
int N = 5;
const vector<double> A{
         -2,6,1,0,0,
         3,2,0,1,0,
         2,-1,0,0,1
};
const vector<double> b{
    40,
    28,
    14
};
const vector<double> c{
    -2,
    -3,
    0,
    0,
    0
};

int main()
{
    SimplexMethod simp = SimplexMethod(A,b,c,M,N);
    simp.run();
    return 0;
}

