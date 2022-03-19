
#include"SimplexMethod_.h"
#include<iostream>
#include<vector>
using namespace std;

int M = 4;
int N = 7;
const vector<double> A{
         20,40,25,1,0,0,0,
         25,50,0,0,1,0,0,
         30,20,40,0,0,1,0,
         10,10,15,0,0,0,1
};
const vector<double> b{
    2000,
    5000,
    4000,
    3000
};
const vector<double> c{
    30,
    15,
    40,
    0,
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

