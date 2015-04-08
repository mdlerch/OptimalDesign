//#include <RcppArmadillo.h>
//#include <algorithm>

#include<armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>

//using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

//double get_delta_g(double, arma::mat, arma::mat);


// use this to print cube
// cube[(int)(iter / pow(size,2))][(int)(iter/size) % size][iter % size];

// [[Rcpp::export]]
int main()
{
    int size = 3;

    // cells that contain zeros and ones
    int cube[size][size][size];

    // generate the most basic latin square
    // i is row
    // j is column
    // k is letter
    int i, j, k;
    
    for(i = 0; i<size; i++)
    {
        for(j = 0; j<size; j++)
        {
            for(k = 0; k<size; k++)
            {
                if((i + j) % size == k)
                {
                    cube[i][j][k] = 1;
		}
                else
                {
		    cube[i][j][k] = 0;
		}
            }
        }
    }

    // this random number generator is kinda lame on my machine
    // i can watch it increase
    srand(time(NULL));
    rand();

    // indexing values
    int z_row;
    int z_col;
    int z_let;

    int w1_let;
    int w2_col;
    int w3_row;

    bool negone = 0;
    int iterations = pow(size, 3);
    int iter;
    
    for(iter = 0; iter < iterations; iter++)
    {
        if(negone == 0) // if we have a proper incidence cube with no -1
        {
	    // choose a random zero-cell
            do
	    {
		z_row = rand() % size;
		z_col = rand() % size;
		z_let = rand() % size;
	    }
	    while(cube[z_row][z_col][z_let]);

	    // find the 1-cells that will make the permutation sub-cube
	    for(k = 0; k < size; k++){
		if(cube[z_row][z_col][k])
		{
		    w1_let = k;
		}
		if(cube[z_row][k][z_let])
		{
		    w2_col = k;
		}
		if(cube[k][z_col][z_let])
		{
		    w3_row = k;
		}
	    }

	    if(cube[w3_row][w2_col][w1_let])
	    {
		cube[z_row][z_col][z_let] = 1;
		cube[z_row][z_col][w1_let] = 0;
		cube[z_row][w2_col][z_let] = 0;
		cube[w3_row][z_col][z_let] = 0;
		cube[z_row][w2_col][w1_let] = 1;
		cube[w3_row][z_col][w1_let] = 1;
		cube[w3_row][w2_col][z_let] = 1;
		cube[w3_row][w2_col][w1_let] = 0;
	    }
	    else
	    {
		cube[z_row][z_col][z_let] = 1;
		cube[z_row][z_col][w1_let] = 0;
		cube[z_row][w2_col][z_let] = 0;
		cube[w3_row][z_col][z_let] = 0;
		cube[z_row][w2_col][w1_let] = 1;
		cube[w3_row][z_col][w1_let] = 1;
		cube[w3_row][w2_col][z_let] = 1;
		cube[w3_row][w2_col][w1_let] = -1;
	    }
        }
	else
	{
	    
	}
    }
}
