/********************************************************************
 * Applied Programming:                                             
 *    Solution of Overdetermined System of Equations Ax=b arising   
 *    in least square problems via QR factorizations using the GSL
 *                                                                  
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "ClassErrors.h"
#include "DynamicArrays.h"

/* Defines */
#define BUFFER (100000)
#define Initial_Size (100)

/* Function prototypes */
void readPoints (FILE *file, DArray *DArrayPtr);
void Norm_FindPoint(int nr, int nc, const DArray *points, gsl_vector *x_ls, 
                                                                 int verbose);
double RSquareError(int nr, int nc, const DArray *points, 
                                                      const gsl_vector *x_ls);
double normOfResiduals(int nr, int nc, const DArray *points, 
                                                      const gsl_vector *x_ls);
double pearson_correl(int nr, int nc, const DArray *points,
                                                      const gsl_vector *x_ls);
double evalPoly(int nc, double x, const gsl_vector *x_ls);


/*****************************************************************************
    Uses least squares to generate approximate functions
    usage: hw8  -norm -order  num   -points  file  [-verbose] \n");
 
 Returns: 0 for success, non-zero for errors
 *****************************************************************************/
int main(int argc, char *argv[])
{
  DArray points;     
  FILE *pointsFile;
  gsl_vector *x_ls; /* least squares solution   */
  
   /* your code here  */
  int verbose = 0;
  int N = 0;
  int order = 0;
  int rc;
  int nc;
  int nr;
  int option_index = 0; // getopt_long stores the option index here.
  char *getoptOptions = "vgno:p:";   //This contains the short command line parameters list 
                             

                  
   struct option long_options[] = {
      /* These options don’t set a flag.
         We distinguish them by their indices. */
      {"verbose",        no_argument,       0, 'v'},
      {"verb",           no_argument,       0, 'v'},
      {"N",         	 no_argument,       0, 'n'},
      {"order",          required_argument,       0, 'o'},
      {"pointsFile",          required_argument,       0, 'p'},
      {0, 0, 0, 0}
      };
   opterr = 1;           /* Enable automatic error reporting */
   while ((rc = getopt_long_only(argc, argv, getoptOptions, long_options, 
                                                &option_index)) != -1) {
        
      /* Detect the end of the options. */
      switch (rc)
        {
         case 'v':                    /* Verbose */
         verbose = 1;
         break;

	 case 'n':                    /* Input */
	 N = 1;          
	 break;

	 case 'o':                    /* Input */
	 order = atof(optarg);     
	
	 if (order < 0){
	 fprintf(stderr, "Order can not be less then 0 /n");
	 return(-1);	
	 }

	 nc = order + 1;     
	 break;

	 case 'p':                    /* Input */
	 pointsFile = fopen(optarg,"r");           
	 break;

	 case '?':  /* Handled by the default error handler */
         break;

       default:
          printf ("Internal error: undefined option %0xX\n", rc);
          exit(PGM_INTERNAL_ERROR);
       }

}
 if ((optind < argc)  ){
      fprintf(stderr, "Tests matrix methods\n");
      fprintf(stderr, "usage: hw8    ‐n[orm]      ‐o[rder] num   ‐p[oints]  file   [‐v[er[bose]]]\n");
      fprintf(stderr, " e.g:  hw8 -n -o 3 -p data.txt \n"); 
      fflush(stderr);
      return(PGM_INTERNAL_ERROR);
	}

	x_ls =gsl_vector_alloc(nc);
	CreateDArray(&points,  Initial_Size );

	readPoints(pointsFile , &points);
	nr = points.EntriesUsed;
	
	if (N == 1){
		Norm_FindPoint(nr, nc, &points, x_ls, verbose);

		printf(" Least Squares Solution via Norm factorization: \n");

		
		for (int k = nc-1 ; k >= 0 ; k--){

			if (k != 0){
				printf("%g x^%d + ",gsl_vector_get(x_ls, k), k);
			}
			else {
				printf("%g", gsl_vector_get(x_ls, k));
			}


		}
		printf("\n");

		printf(" Norm of Residuals:  %f \n ",normOfResiduals(nr, nc, &points, x_ls));
		printf(" RSquare Error: :  %f \n ",RSquareError(nr, nc, &points, x_ls));
		printf(" Pearson Correl:  %f \n ",pearson_correl( nr, nc, &points, x_ls));	
	}

	else {
	}
  
 
 /* Clean up */
   gsl_vector_free(x_ls);  
   DestroyDArray(&points);
   
   fclose(pointsFile);
   
  return PGM_SUCCESS; /* main */
}

/*---------------------------------------------------------------------------
  Find the least squares approximation to data "points" of order "nc" using
  the "Normal equations" approach.
        
                        A'Az = A'b
  
  Where: int nr           - The number of points (rows) of the input file
         int nc           - The number of columns (order) of the solution
         DArray *points   - Pointer to the x,y data
         gsl_vector *x_ls - The solution is returned here
         int verbose      - Verbose flag

  Returns: nothing
  Errors: Assumes the standard GSL error handler
---------------------------------------------------------------------------*/
void Norm_FindPoint(int nr, int nc, const DArray *points, gsl_vector *x_ls, int verbose) {
   //double x;
   //int i, j;         /* counters                 */
   gsl_matrix *A;    /* coefficient matrix A     */
   gsl_matrix *AT;   /* coefficient matrix A'    */
   gsl_matrix *ATA;  /* coefficient matrix A'A   */
   gsl_vector *b;    /* coefficient vector b     */
   gsl_vector *ATB;  /* coefficient vector A'b   */
   gsl_vector *tau;  /* Householder coefficients */
   gsl_vector *res;  /* vector of residuals      */

   /* Allocate space for Matrices and vectors */
   ATA  = gsl_matrix_alloc(nc, nc); /* Data matrix */
   AT   = gsl_matrix_alloc(nc, nr); /* Data matrix */
   A    = gsl_matrix_alloc(nr, nc); /* Data matrix */
   b    = gsl_vector_alloc(nr);     /* Data vector */
   ATB  = gsl_vector_alloc(nc);     /* Data vector */
   tau  = gsl_vector_alloc(nc);     /* required place holder for GSL */
   res  = gsl_vector_alloc(nc);     /* Contains the residual */
   
  /* your code here  */
 	double LOL = 0.0; 
       for (int i=0;i<nr;i++) {
		gsl_matrix_set(A,i,0,1);
         for (int j=1;j<nc;j++) {
		LOL = points -> Payload[i].x;
		for(int z = 1; z < j ; z++){


           	LOL = LOL * points -> Payload[i].x;
		}
			
		gsl_matrix_set(A,i,j,LOL);
	
         } /* for j */
         
         gsl_vector_set(b,i,points -> Payload[i].y); 

       }

 if (verbose == 1){    
   printf("a matrix : ");
for (int k = 0; k < A -> size1; k++){
	for (int i = 0; i < A -> size2; i++){
	printf(" %f ", gsl_matrix_get(A, k, i));
	
	}
	putchar('\n');
	}
} 
 
  gsl_matrix_transpose_memcpy(AT,A);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AT, A, 0.0, ATA); 
  gsl_blas_dgemv(CblasNoTrans, 1.0, AT, b, 0.0, ATB); 
  
  gsl_linalg_QR_decomp(ATA, tau);
  gsl_linalg_QR_lssolve(ATA, tau, ATB, x_ls, res);
 
  /* Free memory  */
  gsl_matrix_free(A);
  gsl_matrix_free(AT);
  gsl_matrix_free(ATA);
  gsl_vector_free(b);
  gsl_vector_free(ATB);
  gsl_vector_free(tau); 
  gsl_vector_free(res);
} /* end Norm_FindPoint() */



/****************************************************************************
  This calculate the norm of residuals given the points and the solution
  Note: You must apply the correction polynomial f(x) to the input data to 
        simulate the correction effect 
  
                   normR = squareRoot [sum {( yi - f(xi))**2} ]

  Where: int nr           - The number of points (rows) of the input file
         int nc           - The number of columns (order) of the solution
         DArray *points   - Pointer to the x,y data
         gsl_vector *x_ls - The solution vector, small power first

  Errors: Assumes the standard GSL error handler
  
  Returns: double norm of residuals
****************************************************************************/
double normOfResiduals(int nr, int nc, const DArray *points, const gsl_vector *x_ls) {

	double eval= 0.0;
	double ans= 0.0;
	double norm= 0.0;

	for(int z = 0; z < nr ; z++){
		eval = evalPoly(nc, points -> Payload[z].x, x_ls);
           	ans += (points -> Payload[z].y - eval)*(points -> Payload[z].y - eval);
	}
	norm = sqrt(ans);

	
	return(norm);
  
} /* normOfResiduals */


/****************************************************************************
  This calculate the R2 coefficient of Determination error between the points  
  and the solution.
  Note: You must apply the correction polynomial f(x) to the input data to 
        simulate the correction effect 
                   
  Where: int nr           - The number of points (rows) of the input file
         int nc           - The number of columns (order) of the solution
         DArray *points   - Pointer to the x,y data
         gsl_vector *x_ls - The solution vector, small power first
 
  Errors: Assumes the standard GSL error handler
  
  Returns: R squared error
****************************************************************************/
double RSquareError(int nr, int nc, const DArray *points, const gsl_vector *x_ls) {

	double temp = 0.0;
	double u= 0.0;
	double eval= 0.0;
	double top= 0.0;
	double bot= 0.0;
	double r2= 0.0;

	for(int z = 0; z < nr ; z++){
           	temp += (points -> Payload[z].y);
		}

		u = (1/nr) * temp;

	for(int i = 0; i < nr ; i++){
		eval = evalPoly(nc, points -> Payload[i].x, x_ls);
		top += (points -> Payload[i].y - eval)*(points -> Payload[i].y - eval);
		bot += ((points -> Payload[i].y - u)*(points -> Payload[i].y - u));
		}
		
		r2 = (1 - (top/bot));
		return(r2);
  
} /* End RSquareError */


/*****************************************************************************
 This calculates the Pearson's Correlation, or the excel function correl()
  Note: You must apply the correction polynomial f(x) to the input data to 
        simulate the correction effect 
                    
  Where: int nr           - The number of points (rows) of the input file
         int nc           - The number of columns (order) of the solution
          DArray *points   - Pointer to the x,y data
         gsl_vector *x_ls - The solution vector, small power first

  Errors: Assumes the standard GSL error handler       
       
 Returns: double pearson_srq 
*****************************************************************************/
double pearson_correl(int nr, int nc, const DArray *points,
                                              const gsl_vector *x_ls) {
	double eval= 0.0;
	double top1= 0.0;
	double top2= 0.0;
	double top3= 0.0;
	double bot1= 0.0;
	double bot2= 0.0;
	double bot3= 0.0;
	double bot4= 0.0;
	double R= 0.0;

	for (int i = 0; i < nr ; i++){
	eval = evalPoly(nc, points -> Payload[i].x, x_ls);
	top1 += nr*(points -> Payload[i].y * eval);
	top2 += (points -> Payload[i].y);
	top3 += eval;
	bot1 += nr*((points -> Payload[i].y)*(points -> Payload[i].y));
	bot3 += nr*(eval*eval);
	}
	
	bot2 = top2*top2;
	bot4 = top3*top3; 
	
	R = (top1 - (top2*top3))/sqrt((bot1 - bot2)*(bot3-bot4));
	return(R);
 
} /* End pearson_correl */




/***************************************************************************************
 Evaluates a polynomial at a point, assumes low order term first.  
 Must use Horner's factorization 
 
 Where: int nc           - The number of columns in the solution
        double x         - Point at which splines should be evaluated
        gsl_vector *x_ls - The solution vector, small power first
          
 Returns: double - The value at the desired point
 Errors:  none
*****************************************************************************************/
double evalPoly(int nc, double x, const gsl_vector *x_ls) {
  
   	double final =  gsl_vector_get(x_ls, nc-1);
	for (int i = nc-2; i >= 0 ; i--){
		
		final = final * x + (gsl_vector_get(x_ls, i));
		
	}
	return(final);
  
} /* End evalPoly */


/***************************************************************************************
 Reads the points data from file and returns it in a Darray
 
 Where: FILE *file     - open handle to read from
                         of the form:     22.0      6.7
                                          23.4      18.8
        DArray *DArrayPtr - Pointer to a dynamic array to store the data
  Returns: nothing
  Errors:  none
*****************************************************************************************/
void readPoints(FILE *file, DArray *DArrayPtr)
{

Data dat;
char String[BUFFER];
char str[] = " \r\n\t";
char *ch;
int count;

   while (fgets(String, BUFFER, file) != NULL)
   {		
	count = 1;
	ch = strtok( String , str);
	while (ch != NULL){
		if( count == 1){
			dat.x = atof(ch);
		}
		if( count == 2){
			dat.y = atof(ch);
		}

		count++;		
		ch = strtok( NULL , str);

	}
	PushToDArray( DArrayPtr, &dat);
     }
     
  return;
} /* readPoints */
