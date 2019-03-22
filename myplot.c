/****************************************************************
 myPlot devicedata output

 Generates a png figure from columnar data.

 Author: R. Repka
 Author: Abhimanyu Gupta
 ***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include "ClassErrors.h"

/* Remove comment to enable debug output */
#define DEBUG 0

/* For posix compatibility */
FILE *popen(const char* cmd, const char *mode);
int pclose(FILE *stream);


int main (int argc, char *argv[]) 
{
  FILE *pipe;                    /* pipe to gnuplot  */

   int rc;
   char *inData = NULL;
   char *outPNG = NULL;
   
   /* getopt_long stores the option index here. */
   int option_index = 0;
  
   /* This contains the short command line parameters list */
   char *getoptOptions = "i:o:";
  
   /* This contains the long command line parameter list, it should mostly 
     match the short list                                                  */
   struct option long_options[] = {
      {"in",   required_argument, 0, 'i'},
      {"out",  required_argument, 0, 'o'},
      {0, 0, 0, 0}
   };
  
   opterr = 1;           /* Enable automatic error reporting */
   while ((rc = getopt_long_only(argc, argv, getoptOptions, long_options, 
                                                    &option_index)) != -1) {
      /* Detect the end of the options. */
      switch (rc)
        {
        case 'i':                   /* input */
          inData = optarg;
          break;

        case 'o':                   /* output */
          outPNG = optarg;
          break;

        case '?':  /* Handled by the default error handler */
          break;

       default:
          printf ("Internal error: undefined option %0xX  '%c'\n", rc, rc);
          exit(PGM_INTERNAL_ERROR);
       } // End switch 
   } /* end while */

   /* Check for command line syntax errors */
   if ((optind < argc) || (inData == NULL) || (outPNG == NULL)){
      /* wrong number of args, display program usage */
      printf("Usage: %s -i[n] datafile -o[ut] plotfile\n", argv[0]);
      printf("    Converts 'datafile' to a PNG file called 'plotfile'\n");
      printf("    datafile format:  realVolt  origData  realVolt newData\n");
      return PGM_SYNTAX_ERROR;
   } /* End if error */   
   printf("Generating plot for '%s'\n", argv[1]); 
      
   /* Open pipe to gnuplot */
   #ifndef DEBUG
   printf(" **** Debug mode enabled ***\n\n");
   pipe = stdout;

   #else
   if( ( pipe = (FILE *) popen("gnuplot -persist","w") ) == NULL ) {
         printf("Unable to open pipe to gnuplot\n");
         return PGM_SYNTAX_ERROR;
   }
   #endif
 
   
   /* Gnuplot commands */
   /* Set output to png file, can change to svg or x11 */
 	fprintf(pipe, "set terminal png enhanced font 'DejaVuSans.ttf' 12 \n");

   /* Set the output file name */
	fprintf(pipe, "set output '%s' \n", outPNG);

   /* Set plot attributes */ 
	fprintf(pipe, "set key left \n");
	fprintf(pipe, "set key on \n ");
	fprintf(pipe, "set border 3 \n");
	fprintf(pipe, "set style data lines \n");
    
   /* Title, labels and legend */
  	fprintf(pipe, "set title 'Sensor Response Curve'\n");
	fprintf(pipe, "set ylabel 'Digital Output'\n");
	fprintf(pipe, "set xlabel 'Analog Input'\n");
	
   /* Read data from file and generate plot */ 
   	fprintf(pipe,"plot '%s' using 1:2 lt 1 lw 2 t'Real', '%s' using 1:4 lt 2 lw 2 t'Correction', '%s' using 1:3 lt 3 lw 2 t'Ideal'\n", inData, inData, inData);
 
   /* clean up: close pipe and return success code */
   pclose (pipe);
    
  return 0;
} /* main */
