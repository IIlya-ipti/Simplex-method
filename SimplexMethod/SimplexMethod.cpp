
#include"SimplexMethod_.h"
#include "Matrix.h"

int M = 3;
int N = 7;
const std::vector<double> A{
         -1,1,0,0,1,0,0,
         0,0,1,-1,0,1,0,
         -2,2,-3,3,0,0,1
        
};
const std::vector<double> b{
    0.4,
    5,
    0
};
const std::vector<double> c{
    10,
   -10,
    2,
    -2,
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

