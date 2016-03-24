//#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define CHUNKSIZE 50
#define num_threads 32
#define N     1000
#define PI  3.14159621234161928
/* 
COMPILE
============================
CC=icc
CFLAGS=-Wall -O3 -openmp

HW3ParallelCorrelationFunc1_0 : HW3ParallelCorrelationFunc1_0.c
   $(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
   rm -f HW3ParallelCorrelationFunc1_0
============================

RUN
============================
#!/bin/bash

if [[ $# -ne 1 ]] ; then
    echo "usage: ./HW3ParallelCorrelationFunc1_0 num_threads"
    exit 1
fi
export OMP_NUM_THREADS=$1
./HW3ParallelCorrelationFunc1_0
=============================
*/

/*
 * ------------------------------------------------------------------------
 *       ___           _     __  __      _   _           
 *      |   \ __ _ _ _| |__ |  \/  |__ _| |_| |_ ___ _ _ 
 *      | |) / _` | '_| / / | |\/| / _` |  _|  _/ -_) '_|
 *      |___/\__,_|_| |_\_\ |_|  |_\__,_|\__|\__\___|_|  
 *                     
#  .-----..-.  .-. .---.              .-.-.  .---. .-..-. .-..-----.    
#  `-' '-'| {  } |/ {-. \     ___     | } }}/ {-. \{ ||  \{ |`-' '-'    
#    } {  {  /\  }\ '-} /    {___}    | |-' \ '-} /| }| }\  {  } {      
#    `-'  `-'  `-' `---'              `-'    `---' `-'`-' `-'  `-'      
#  .----. .---. .---. .---. .----..-.     .--. .-----..-. .---. .-. .-. 
#  | }`-'/ {-. \} }}_}} }}_}} |__}} |    / {} \`-' '-'{ |/ {-. \|  \{ | 
#  | },-.\ '-} /| } \ | } \ } '__}} '--./  /\  \ } {  | }\ '-} /| }\  { 
#  `----' `---' `-'-' `-'-' `----'`----'`-'  `-' `-'  `-' `---' `-' `-' 
#  .----..-. .-..-. .-..----..-----..-. .---. .-. .-.                   
#  } |__}| } { ||  \{ || }`-'`-' '-'{ |/ {-. \|  \{ |                   
#  } '_} \ `-' /| }\  {| },-.  } {  | }\ '-} /| }\  {                   
#  `--'   `---' `-' `-'`----'  `-'  `-' `---' `-' `-'                   
#                                                                       
# This program will make the correlation function given files with [RA, Dec, Z]. 
# It will use the estimator:
#      {   Xi(r) = (N_r/N_d)**2 (DD(r)/RR(r) - 1)   }
#
# The Sample data is separated into 15 logarithmic bins of separation in range 
# 0.1 - 20 h**-1 Mpc ( the first bin at log r = -1 and the last bin is at 
# log r = 1.301 )
#
# The Random Points are redshifts in range 0.02 <= z <= 0.06. Using the SDSS sky
# coverage. 
#
#
# FILE NAMES: 
#    Part 1:
#      Galaxy Samples: 
#      { SDSS_Mr21_rspace.dat, SDSS_Mr20_rspace.dat, SDSS_Mr20_zspace.dat }
#      Random Points: 
#      { SDSS_random.dat }
#
#     A.  log Xi vs. log r  for two real-space galaxies. {-20, -21}
#     B.  log Xi vs. log r  for real-space vs. redshift-space sample. {-20}
#     
#  Part 2: Use the DM files to find Xi(r).
#     DM Samples: Points in Cartesian Coordinates {Xi, Yi, Zi,} Vol=141.3h-1Mpc
#     {DM.dat}
#     Random Points:
#     {DM_random.dat}
#     
#     B. Compute the bias function: b(r) = SQRT[ Xi_gal / Xi_DM  ]

------------------------------------------------------------------------
*/

int main ()  
{
   printf("FIRST of MAIN\n");
   long int i, j;


   int chunk;
   //float a[N], b[N], c[N];



   double dec_to_rad = (2. * PI / 360.);
   double deg_to_rad = (2. * PI / 360.);

   //Create the Distance Array. 
   const int num_bins = 15;
   double logr_min = -1.0;
   double logr_max = 1.3011;
   long int errorcounts = 0;
   //    r                 0.1,                 ...                  , 20 
   // logr                 -1,                  ...                  ,1.3011
   // index                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
   long int distance_counts[15] = { 0 };
   long int randdistance_counts[15] = { 0 };
   
   double Xi_func[15] = { 0.0 };
   int dist_index;

   //This should be set below 5495 to take a limited sample. 
   //int FILELENGTH = 4000;

   //long int FILELENGTH20r  = 28162; 
   //long int FILELENGTH20z  = 28383; 
   //long int FILELENGTH21r  = 5495;  
   //long int FILELENGTHrand = 42654; 
   long int FILELENGTHDM = 50000;
   //long int FILELENGTHDMrand = 100000;
   long int FILELENGTHDMrand = 100000;


   //string datafilename = 'SDSS_random.dat'
   //static const char random_datafile[] = "SDSS_random.dat"; 

   // =========================================================
   /***
    *      ___           _     __  __      _   _           
    *     |   \ __ _ _ _| |__ |  \/  |__ _| |_| |_ ___ _ _ 
    *     | |) / _` | '_| / / | |\/| / _` |  _|  _/ -_) '_|
    *     |___/\__,_|_| |_\_\ |_|  |_\__,_|\__|\__\___|_|  
    *                                                      
    */

   
   static const char DM_datafile[] = "/home/chasonn/LargeScaleHW/HW3Codes/DMCorrelationFunctionParallel/DM.dat"; 
   FILE *myfile = fopen ( DM_datafile, "r" );
   if (myfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   printf("Opened file - Begining assignment.\n");

   //long int N_data = 1000;
   long int N_data = FILELENGTHDM;
   
   double X_LIST[N_data];
   double Y_LIST[N_data];
   double Z_LIST[N_data];
   

   //==================
   // READ IN Mr21 DATA. 
   //==================
   //i = 0;
   //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
   //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF && i < N_data){
   
   for(i = 0; i < (N_data); i++){
      fscanf(myfile, "%lf", &X_LIST[i]);
      fscanf(myfile, "%lf", &Y_LIST[i]);
      fscanf(myfile, "%lf", &Z_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N_data-2)){
         printf("Close or exceeded N_data limit. X: %lf \n", X_LIST[i]);
         //break;

      }
      //++i;
   }

   fclose(myfile);
   printf("Closing File.\n");

   //Begin Nested For Loop to Count Pairs. 

   
   //==========================
   // COMPUTE Mr 21 DATA COUNTS
   //==========================
   printf("Beginning Nested Loops...\n");
   
   double D, logD; 
   double r = 0; 
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;  

   double rj = 0; 
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;  

   chunk = CHUNKSIZE;

   #pragma omp parallel shared( Z_LIST, Y_LIST, X_LIST, N_data, chunk) private (D,\
   logD, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      
      //OMP_NUM_THREADS = 16; 
      omp_set_num_threads(num_threads);
      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for schedule(guided, chunk)
      for(i=0; i < (N_data-1); ++i){
         
         
         x1 = X_LIST[i]; 
         y1 = Y_LIST[i]; 
         z1 = Z_LIST[i];

         
         for(j=0; j < (N_data-1); ++j){
            if ( j!=i ){
               
               x2 = X_LIST[j]; 
               y2 = Y_LIST[j]; 
               z2 = Z_LIST[j];

               //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
               D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
               logD = log10(D);
               
               //dist_index = (logD +1)*(num_bins/(logr_max +1))
               dist_index = (int) (floor((logD - logr_min)*(num_bins/(logr_max - logr_min))));            
               if (dist_index >= 0 && dist_index < num_bins){
                  //Increment the appropiate bin.
                  if (dist_index > 14)
                     printf("YELLING!");
                  /*
                  //OLD NON-MULTITHREADED:
                  distance_counts[dist_index] += 1;
                  */
                  sum_local_counts[dist_index] += 1;
               } 
            }
         
         //end inner for loop
         }
      //end outer for loop
      }

      #pragma omp critical
      {
         //Sum up over the local counts on each thread to get total distance count.
         
         for(i=0 ; i < num_bins; ++i)
         {
            distance_counts[i] += sum_local_counts[i];
         }

      }

   //END PRAGMA PARALLEL BLOCK
   }

   printf("\n*");
   printf("\n   *");
   printf("\n     *");
   printf("\n       *");
   printf("\n      *");
   printf("\n     *");
   printf("\n    *");
   printf("\n   *        *");
   printf("\n   *");
   printf("\n     *");
   printf("\n      *");
   printf("\n       **");
   printf("\n        * *");
   printf("\n        * * *");
   printf("\n       * * * *\n");

   printf("****************************\n");
   printf("FINISHED PRAGMA OMP CRITICAL\n");
   printf("****************************\n");


   printf("FINISHED DM NESTED LOOP. \n");

   printf("Dividing Counts by two to correct double counting...");
   printf("Counts: ");
   for(i=0 ; i < num_bins; ++i)
   {
      distance_counts[i] = (long long) (floor(distance_counts[i]/2.)) ;
      printf("%ld ", distance_counts[i]);
   }
   /*
   for(i =0; i< num_bins; ++i){
      printf("%ld ", distance_counts[i]);
   }
   */
   printf("\n") ;  


 

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   /***
    *      ___    _   _  _ ___   ___  __  __ 
    *     | _ \  /_\ | \| |   \ / _ \|  \/  |
    *     |   / / _ \| .` | |) | (_) | |\/| |
    *     |_|_\/_/ \_\_|\_|___/ \___/|_|  |_|
    *                                        
    */


   // INIITIALIZE FOR RANDOM COUNTS
   static const char random_datafile[] = "/home/chasonn/LargeScaleHW/HW3Codes/DMCorrelationFunctionParallel/DM_random.dat"; 
   //static const char r21_datafile[] = "SDSS_Mr21_rspace.dat"; 
   // Stores the number of points in random set. 
   long int N_rand = FILELENGTHDMrand;

   FILE *myrandfile = fopen ( random_datafile, "r" );
   if (myrandfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }
   
   //double * randX_LIST[N_rand] = malloc(N_rand * sizeof(double));
   //double * randY_LIST[N_rand] = malloc(N_rand * sizeof(double));
   //double * randZ_LIST[N_rand] = malloc(N_rand * sizeof(double));
   double randX_LIST[N_rand]; 
   double randY_LIST[N_rand]; 
   double randZ_LIST[N_rand]; 
   
   //===============
   // READ IN DATA. 
   //===============
   for(i = 0; i < (N_rand); ++i){
      fscanf(myrandfile, "%lf", &randX_LIST[i]);
      fscanf(myrandfile, "%lf", &randY_LIST[i]);
      fscanf(myrandfile, "%lf", &randZ_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N_rand-2)){
         printf("Close or exceeded N_data limit. X: %lf \n", randX_LIST[i]);
         //break;
      }
      //++i;
   }

   fclose(myrandfile);
   printf("Closing File.\n");

   //Begin Nested For Loop to Count Pairs. 
   //=========================
   // COMPUTE RAND DATA COUNTS
   //=========================
   printf("Beginning Random Nested Loops...\n");


   #pragma omp parallel shared( randZ_LIST, randY_LIST, randX_LIST, N_rand,chunk) private (D,\
   logD, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      
      omp_set_num_threads(num_threads);
      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for schedule(guided, chunk) 
      for(i=0; i< (N_rand-1); ++i){
         

         x1 = randX_LIST[i];
         y1 = randY_LIST[i];
         z1 = randZ_LIST[i];

         
         for(j=0; j< (N_rand-1); ++j){
            if (j!=i){
               x2 = randX_LIST[j];
               y2 = randY_LIST[j];
               z2 = randZ_LIST[j];

               //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
               D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
               logD = log10(D);
               
               //dist_index = (logD +1)*(num_bins/(logr_max +1))
               dist_index = (int) (floor((logD - logr_min)*(num_bins/(logr_max - logr_min))));            
               if (dist_index >= 0 && dist_index < num_bins){
                  //Increment the appropiate bin.
                  if (dist_index > 14)
                     printf("YELLING!");
                  //randdistance_counts[dist_index] += 1;
                  sum_local_counts[dist_index] += 1;
               } 
            }
         
         //end inner for loop
         }
      //end outer for loop
      }

      #pragma omp critical
      {
         //Sum up over the local counts on each thread to get total distance count.
         
         for(i=0 ; i < num_bins; ++i)
         {
            randdistance_counts[i] += sum_local_counts[i];
         }

      }

   //END PRAGMA PARALLEL BLOCK
   }

   printf("FINISHED RANDOM NESTED LOOPS! \n");
   printf("Counts: ");

   for(i =0; i< num_bins; ++i){
      randdistance_counts[i] = (long long) (floor(randdistance_counts[i]/2.)) ;
      printf("%ld ", randdistance_counts[i]);

   }
   printf("\n") ;  

   //COUNTS COMPLETED.  {random, 21r, 20r, 20z}
   //==============================================================================
   /***
    *       ___ ___  _   _ _  _ _____ ___    ___ ___  __  __ ___ _    ___ _____ ___ 
    *      / __/ _ \| | | | \| |_   _/ __|  / __/ _ \|  \/  | _ \ |  | __|_   _| __|
    *     | (_| (_) | |_| | .` | | | \__ \ | (_| (_) | |\/| |  _/ |__| _|  | | | _| 
    *      \___\___/ \___/|_|\_| |_| |___/  \___\___/|_|  |_|_| |____|___| |_| |___|
    *      -------------------------------   ---------------------------------------                                                                         
    */
   //==============================================================================

   //==============================================================================
   /***
    *       ___ ___  __  __ ___ _   _ _____ ___                                 
    *      / __/ _ \|  \/  | _ \ | | |_   _| __|                                
    *     | (_| (_) | |\/| |  _/ |_| | | | | _|                                 
    *      \___\___/|_|_ |_|_|__\___/  |_|_|___|___ ___  _  _                   
    *      / __/ _ \| _ \ _ \ __| |    /_\_   _|_ _/ _ \| \| |                  
    *     | (_| (_) |   /   / _|| |__ / _ \| |  | | (_) | .` |                  
    *      \___\___/|_|_\_|_\___|____/_/ \_\_|_|___\___/|_|\_|_   ___  __  __   
    *     | __| | | | \| |/ __|_   _|_ _/ _ \| \| |   | _ )/ _ \ / _ \|  \/  |  
    *     | _|| |_| | .` | (__  | |  | | (_) | .` |_  | _ \ (_) | (_) | |\/| |_ 
    *     |_|  \___/|_|\_|\___| |_| |___\___/|_|\_(_) |___/\___/ \___/|_|  |_(_)
    *                                                                           
    */
   //==============================================================================

   // =================================
   // COMPUTE THE CORRELATION FUNCTION
   // =================================

   //DM.
   printf("Calculating DM Correlation Function...\n");
   double ratio = (double) N_rand / (double) N_data;
   for(i=0; i<num_bins; i++){
      //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
      Xi_func[i] = ratio * ratio *  ( (double) distance_counts[i] / (double) randdistance_counts[i]) - 1.0 ;
      printf("%f ", Xi_func[i]); 
   }
   printf("\n");

   //----------------------------------------------------
   /***
    *      ___                 _         ___ ___ _    ___ 
    *     / __| __ ___ _____  | |_ ___  | __|_ _| |  | __|
    *     \__ \/ _` \ V / -_) |  _/ _ \ | _| | || |__| _| 
    *     |___/\__,_|\_/\___|  \__\___/ |_| |___|____|___|
    *                                                     
    */
   //----------------------------------------------------

   // ==========================
   // WRITING DM RESULTS TO FILE
   // ==========================
   printf("Saving DM counts to file.\n");
   FILE *fp_out;
   fp_out = fopen("output_DMcounts.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", distance_counts[i]);
   } 
   fclose(fp_out);


   // ==========================
   // WRITING random RESULTS TO FILE
   // ==========================
   printf("Saving Random counts to file.\n");
   //FILE *fp_out;
   fp_out = fopen("output_counts_DMrandom.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", randdistance_counts[i]);
   } 
   fclose(fp_out);



   // ===============================
   // WRITING DM Xi RESULTS TO FILE
   // ===============================
   printf("Saving Xi-DM to file.\n");
   //FILE *fp_out;
   fp_out = fopen("output_Xi_DM.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%f \n", Xi_func[i]);
   } 
   fclose(fp_out);


   return 0;

}



//END
