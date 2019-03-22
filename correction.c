/******************************************************************
 * Program to correct the data for the sensor
 * Note: Each student gets unique data, so this exact equation will
 * differ for each students solution
 * Be sure to use Honer's factorization  
 * ***************************************************************/

#include <stdio.h>
#include <stdlib.h>

/* Runs the data through the fitting line */

int main(int argc, char *argv[])
{
    int res, real, ideal;
    
	double i = 0.0;
    while(scanf("%d %d", &ideal, &real) != EOF)
    {

	i = (((-8.1489e-12*real + 4.85272e-08)*real + -0.00010079)*real + 0.231156)*real+ 13.0225;

	res = i >= 0 ? (int)(i+0.5): (int)(i-0.5);
	
	
        printf("%d %d\n", ideal, real - res);
    }
    return 0;
}
