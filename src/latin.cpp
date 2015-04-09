#include <RcppArmadillo.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::umat opt_genLatin(int size, int iterations)
{
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

    // i think that this RNG is suspect
    srand(time(NULL));
    rand();

    // indexing values
    // z is your 'focus point' on the sub-cube that gets switched
    int z_row;
    int z_col;
    int z_let;
    // w1, w2, w3 are the 1-cells that are the other corners of the sub-cube
    int w1_let;
    int w2_col;
    int w3_row;

    // indicator of whether the current cube is proper (has only 0's and 1's)
    // or is improper (has a -1). A proper cube represents a latin square
    bool negone = 0;
    
    int iter;
    for(iter = 0; iter < iterations || negone; iter++)
    {
        // if we have a proper incidence cube
        if(negone == 0)
        {
            // choose a random 0-cell
            do
            {
                z_row = rand() % size;
                z_col = rand() % size;
                z_let = rand() % size;
            }
            while(cube[z_row][z_col][z_let]);

            // find the 1-cells that will make the permutation sub-cube
            for(k = 0; k < size; k++){
                if(cube[z_row][z_col][k] == 1)
                {
                    w1_let = k;
                }
                if(cube[z_row][k][z_let] == 1)
                {
                    w2_col = k;
                }
                if(cube[k][z_col][z_let] == 1)
                {
                    w3_row = k;
                }
            }

            // if opposing ?-cell is a 1-cell then make this permutation
            if(cube[w3_row][w2_col][w1_let] == 1)
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
            // if opposing ?-cell is a 0-cell then make this permutation
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
                negone = 1;
            }
        }
        // if we have an improper incidence cube
        else
        {
            // the last permutation we made led us to this state. so
            // let z now be the guy that got switched to a -1
            z_row = w3_row;
            z_col = w2_col;
            z_let = w1_let;
            
            // find the 1-cells that will make the permutation sub-cube.
            // there are two for each row, column, and letter now so you gotta
            // randomly choose for each
            int a = rand()%2;
            int b = rand()%2;
            int c = rand()%2;
            
            for(k = 0; k < size; k++)
            {
                if(cube[z_row][z_col][k] == 1)
                {
                    if(a)
                    {
                        w1_let = k;
                        a = 0;
                    } else
                    {
                        a = 1;
                    }
                }
                if(cube[z_row][k][z_let] == 1)
                {
                    if(b)
                    {
                        w2_col = k;
                        b = 0;
                    } else
                    {
                        b = 1;
                    }
                }
                if(cube[k][z_col][z_let] == 1)
                {
                    if(c)
                    {
                        w3_row = k;
                        c = 0;
                    } else
                    {
                        c = 1;
                    }
                }
            }

            // if opposing ?-cell is a 1-cell then make this permutation
            if(cube[w3_row][w2_col][w1_let] == 1)
            {
                cube[z_row][z_col][z_let] = 0;
                cube[z_row][z_col][w1_let] = 0;
                cube[z_row][w2_col][z_let] = 0;
                cube[w3_row][z_col][z_let] = 0;
                cube[z_row][w2_col][w1_let] = 1;
                cube[w3_row][z_col][w1_let] = 1;
                cube[w3_row][w2_col][z_let] = 1;
                cube[w3_row][w2_col][w1_let] = 0;
                negone = 0;
            }
            // if opposing ?-cell is a 0-cell then make this permutation
            else
            {
                cube[z_row][z_col][z_let] = 0;
                cube[z_row][z_col][w1_let] = 0;
                cube[z_row][w2_col][z_let] = 0;
                cube[w3_row][z_col][z_let] = 0;
                cube[z_row][w2_col][w1_let] = 1;
                cube[w3_row][z_col][w1_let] = 1;
                cube[w3_row][w2_col][z_let] = 1;
                cube[w3_row][w2_col][w1_let] = -1;
            }
        }
    }
    
    arma::umat latin_square(size, size);

    for(i = 0; i<size; i++)
    {
        for(j = 0; j<size; j++)
        {
            for(k = 0; k<size; k++)
            {
                if(cube[i][j][k]) latin_square(i,j) = k;
            }
        }
    }

    return latin_square;
}
