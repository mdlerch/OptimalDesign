//#include <RcppArmadillo.h>
//#include <algorithm>

#include<armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

//using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

//double get_delta_g(double, arma::mat, arma::mat);

// [[Rcpp::export]]
int main()
{
    int size = 10;

    // cells that contain zeros and ones
    int zeros[size*size*(size-1)][3];
    int ones[size*size][3];

    // generate the most basic latin square
    // there's gotta be a more susinct way to do this
    int o_fill = 0;
    int z_fill = 0;
    int i, j, k;
    
    for(i = 0; i<size; i++)
    {
        for(j = 0; j<size; j++)
        {
            for(k = 0; k<size; k++)
            {
                if((i + j) % size == k)
                {
                    ones[o_fill][0] = i;
                    ones[o_fill][1] = j;
                    ones[o_fill][2] = k;

                    o_fill++;
                }
                else
                {
                    zeros[z_fill][0] = i;
                    zeros[z_fill][1] = j;
                    zeros[z_fill][2] = k;

                    z_fill++;
                }
            }
        }
    }

    // this random number generator is kinda lame on my machine
    // i can watch it increase
    srand(time(NULL));
    rand();

    int z_index;
    int negone[3];
    int iterations = 100;
    
    for(i = 0; i < iterations; i++)
    {
        if(negone[2] == -1){
            // choose a random zero-cell
            z_index = rand() % (size * size * (size - 1));

            // find the one-cells that have the same row, column, or letter as z
            for(j = 0; j < size; j++)
            {
            }
        }
    }
}
