//#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define CHUNKSIZE 50
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
------------------------------------------------------------------------

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
   printf("Line 95 - FIRST of MAIN\n");
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
   printf("line104\n");
   //    r                 0.1,                 ...                  , 20 
   // logr                 -1,                  ...                  ,1.3011
   // index                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
   long int distance_counts[15] = { 0 };
   long int distance20_counts[15] = { 0 };
   long int distance20z_counts[15] = { 0 };
   long int randdistance_counts[15] = { 0 };
   
   double Xi_func[15] = { 0.0 };
   double Xi20_func[15] = { 0.0 };
   double Xi20z_func[15] = { 0.0 };
   printf("Line 118\n");

   int dist_index;

   //This should be set below 5495 to take a limited sample. 
   //int FILELENGTH = 4000;

   long int FILELENGTH20r  = 28162; 
   long int FILELENGTH20z  = 28383; 
   long int FILELENGTH21r  = 5495;  
   long int FILELENGTHrand = 42654; 

   //string datafilename = 'SDSS_random.dat'
   //static const char random_datafile[] = "SDSS_random.dat"; 

   // =========================================================
   /***
    *      __  __       ___ _   ___          _ 
    *     |  \/  |_ _  |_  ) | | _ \___ __ _| |
    *     | |\/| | '_|  / /| | |   / -_) _` | |
    *     |_|  |_|_|   /___|_| |_|_\___\__,_|_|
    *                                          
    */

   
   static const char r21_datafile[] = "SDSS_Mr21_rspace.dat"; 
   FILE *myfile = fopen ( r21_datafile, "r" );
   if (myfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   printf("Opened file - Begining assignment.\n");

   //long int N_data = 1000;
   long int N_data = FILELENGTH21r;
   
   double RA_LIST[N_data];
   double DEC_LIST[N_data];
   double Z_LIST[N_data];
   
   /*
   float * Z_LIST   = malloc(N_data * sizeof(float));
   float * DEC_LIST = malloc(N_data * sizeof(float));
   float * RA_LIST  = malloc(N_data * sizeof(float));
   float RA;
   float DEC;
   float Z;
   */

   //==================
   // READ IN Mr21 DATA. 
   //==================
   //i = 0;
   //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
   //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF && i < N_data){
   
   for(i = 0; i < (N_data); i++){
      fscanf(myfile, "%lf", &RA_LIST[i]);
      fscanf(myfile, "%lf", &DEC_LIST[i]);
      fscanf(myfile, "%lf", &Z_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N_data-2)){
         printf("Close or exceeded N_data limit. RA: %lf \n", RA_LIST[i]);
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
   double DEC_rad = 0;
   double RA_rad = 0;
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;  

   double rj = 0; 
   double DEC_radj = 0;
   double RA_radj = 0;
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;  

   chunk = CHUNKSIZE;

   #pragma omp parallel shared( Z_LIST, DEC_LIST, RA_LIST, N_data, deg_to_rad ,chunk) private (D,\
   logD, r, rj, DEC_rad, DEC_radj, RA_rad, RA_radj, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      

      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for
      for(i=0; i < (N_data-1); ++i){
         
         r = 2998. * Z_LIST[i];
         DEC_rad = DEC_LIST[i] * (deg_to_rad);
         RA_rad  = RA_LIST[i]  * (deg_to_rad);

         x1 = r * cos( DEC_rad ) * cos( RA_rad ); 
         y1 = r * cos( DEC_rad ) * sin( RA_rad ); 
         z1 = r * sin( DEC_rad );

         
         for(j=0; j < (N_data-1); ++j){
            if ( j!=i ){
               rj = 2998. * Z_LIST[j];
               DEC_radj = DEC_LIST[j] * (deg_to_rad);
               RA_radj  = RA_LIST[j]  * (deg_to_rad);

               x2 = rj * cos( DEC_radj ) * cos( RA_radj ); 
               y2 = rj * cos( DEC_radj ) * sin( RA_radj ); 
               z2 = rj * sin( DEC_radj );

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


   printf("FINISHED Mr21 NESTED LOOP. \n");

   printf("Dividing Counts by two to correct double counting...");
   for(i=0 ; i < num_bins; ++i)
   {
      distance_counts[i] = (long long) (floor(distance_counts[i]/2.)) ;
      printf("%ld ", distance_counts[i]);
   }
   printf("Counts: ");

   /*
   for(i =0; i< num_bins; ++i){
      printf("%ld ", distance_counts[i]);
   }
   */
   printf("\n") ;  


   //===================================================
   /***
    *      __  __       ___ __    ___          _ 
    *     |  \/  |_ _  |_  )  \  | _ \___ __ _| |
    *     | |\/| | '_|  / / () | |   / -_) _` | |
    *     |_|  |_|_|   /___\__/  |_|_\___\__,_|_|
    *                                            
    */

   static const char r20_datafile[] = "SDSS_Mr20_rspace.dat"; 
   FILE *my20file = fopen ( r20_datafile, "r" );
   if (my20file == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   printf("Opened file - Begining assignment.\n");

   //long int N_data = 1000;
   long int N20_data = FILELENGTH20r;
   
   double RA20_LIST[N20_data];
   double DEC20_LIST[N20_data];
   double Z20_LIST[N20_data];
   
   //==================
   // READ IN Mr20 DATA. 
   //==================
   //i = 0;
   //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
   //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF && i < N_data){
   
   for(i = 0; i < (N20_data); i++){
      fscanf(my20file, "%lf", &RA20_LIST[i]);
      fscanf(my20file, "%lf", &DEC20_LIST[i]);
      fscanf(my20file, "%lf", &Z20_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N20_data-2)){
         printf("Close or exceeded N20_data limit. RA: %lf \n", RA20_LIST[i]);
         //break;

      }
      //++i;
   }

   fclose(my20file);
   printf("Closing File.\n");

   //Begin Nested For Loop to Count Pairs. 

   
   //==========================
   // COMPUTE Mr 21 DATA COUNTS
   //==========================
   printf("Beginning Nested Loops...\n");
   /*
   double D, logD; 
   double r = 0; 
   double DEC_rad = 0;
   double RA_rad = 0;
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;  

   double rj = 0; 
   double DEC_radj = 0;
   double RA_radj = 0;
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;  
   */                      
   #pragma omp parallel shared( Z20_LIST, DEC20_LIST, RA20_LIST, N20_data, deg_to_rad ,chunk) private (D,\
   logD, r, rj, DEC_rad, DEC_radj, RA_rad, RA_radj, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      
      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for
      for(i=0; i< (N20_data-1); ++i){
         
         r = 2998. * Z20_LIST[i];
         DEC_rad = DEC20_LIST[i] * (deg_to_rad);
         RA_rad  = RA20_LIST[i]  * (deg_to_rad);

         x1 = r * cos( DEC_rad ) * cos( RA_rad ); 
         y1 = r * cos( DEC_rad ) * sin( RA_rad ); 
         z1 = r * sin( DEC_rad );

         for(j=0; j< (N20_data-1); ++j){
            if (j!=i){
               rj = 2998. * Z20_LIST[j];
               DEC_radj = DEC20_LIST[j] * (deg_to_rad);
               RA_radj  = RA20_LIST[j]  * (deg_to_rad);

               x2 = rj * cos( DEC_radj ) * cos( RA_radj ); 
               y2 = rj * cos( DEC_radj ) * sin( RA_radj ); 
               z2 = rj * sin( DEC_radj );

               //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
               D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
               logD = log10(D);
               
               //dist_index = (logD +1)*(num_bins/(logr_max +1))
               dist_index = (int) (floor((logD - logr_min)*(num_bins/(logr_max - logr_min))));            
               if (dist_index >= 0 && dist_index < num_bins){
                  //Increment the appropiate bin.
                  if (dist_index > 14)
                     printf("YELLING!");
                  //distance20_counts[dist_index] += 1;
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
            distance20_counts[i] += sum_local_counts[i];
         }

      }

   //END PRAGMA PARALLEL BLOCK
   }

   printf("FINISHED Mr20 NESTED LOOP. \n");
   printf("Counts: ");

   for(i =0; i< num_bins; ++i){
      distance20_counts[i] = (long long) (floor(distance20_counts[i]/2.)) ;
      printf("%ld ", distance20_counts[i]);
   }
   printf("\n") ;  



   //=========================================================
   /***
    *      __  __       ___ __    ____                        
    *     |  \/  |_ _  |_  )  \  |_  /__ ____ __  __ _ __ ___ 
    *     | |\/| | '_|  / / () |  / /___(_-< '_ \/ _` / _/ -_)
    *     |_|  |_|_|   /___\__/  /___|  /__/ .__/\__,_\__\___|
    *                                      |_|                
    */

   static const char r20z_datafile[] = "SDSS_Mr20_zspace.dat"; 
   FILE *my20zfile = fopen ( r20z_datafile, "r" );
   if (my20zfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   printf("Opened file - Begining assignment.\n");

   //long int N_data = 1000;
   long int N20z_data = FILELENGTH20z;
   
   double RA20z_LIST[N20z_data];
   double DEC20z_LIST[N20z_data];
   double Z20z_LIST[N20z_data];
   
   //==================
   // READ IN Mr20z DATA. 
   //==================
   //i = 0;
   //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
   //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF && i < N_data){
   
   for(i = 0; i < (N20z_data); i++){
      fscanf(my20zfile, "%lf", &RA20z_LIST[i]);
      fscanf(my20zfile, "%lf", &DEC20z_LIST[i]);
      fscanf(my20zfile, "%lf", &Z20z_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N20z_data-2)){
         printf("Close or exceeded N20z_data limit. RA: %lf \n", RA20z_LIST[i]);
         //break;

      }
      //++i;
   }

   fclose(my20zfile);
   printf("Closing File.\n");

   //Begin Nested For Loop to Count Pairs. 

   
   //==========================
   // COMPUTE Mr 20z DATA COUNTS
   //==========================
   printf("Beginning Nested Loops...\n");
   /*
   double D, logD; 
   double r = 0; 
   double DEC_rad = 0;
   double RA_rad = 0;
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;  

   double rj = 0; 
   double DEC_radj = 0;
   double RA_radj = 0;
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;  
   */                      

   #pragma omp parallel shared( Z20z_LIST, DEC20z_LIST, RA20z_LIST, N20z_data, deg_to_rad ,chunk) private (D,\
   logD, r, rj, DEC_rad, DEC_radj, RA_rad, RA_radj, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      
      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for
      for(i=0; i< (N20z_data-1); ++i){
         
         r = 2998. * Z20z_LIST[i];
         DEC_rad = DEC20z_LIST[i] * (deg_to_rad);
         RA_rad  = RA20z_LIST[i]  * (deg_to_rad);

         x1 = r * cos( DEC_rad ) * cos( RA_rad ); 
         y1 = r * cos( DEC_rad ) * sin( RA_rad ); 
         z1 = r * sin( DEC_rad );

         for(j=0; j< (N20z_data-1); ++j){
            if (j!=i){
               rj = 2998. * Z20z_LIST[j];
               DEC_radj = DEC20z_LIST[j] * (deg_to_rad);
               RA_radj  = RA20z_LIST[j]  * (deg_to_rad);

               x2 = rj * cos( DEC_radj ) * cos( RA_radj ); 
               y2 = rj * cos( DEC_radj ) * sin( RA_radj ); 
               z2 = rj * sin( DEC_radj );

               //D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
               D = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
               logD = log10(D);
               
               //dist_index = (logD +1)*(num_bins/(logr_max +1))
               dist_index = (int) (floor((logD - logr_min)*(num_bins/(logr_max - logr_min))));            
               if (dist_index >= 0 && dist_index < num_bins){
                  //Increment the appropiate bin.
                  if (dist_index > 14)
                     printf("YELLING!");
                  //distance20z_counts[dist_index] += 1;
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
            distance20z_counts[i] += sum_local_counts[i];
         }

      }

   //END PRAGMA PARALLEL BLOCK
   }

   printf("FINISHED Mr20z NESTED LOOP. \n");
   printf("Counts: ");

   for(i =0; i< num_bins; ++i){
      distance20z_counts[i] = (long long) (floor(distance20z_counts[i]/2.)) ;
      printf("%ld ", distance20z_counts[i]);
   }
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
   static const char random_datafile[] = "SDSS_random.dat"; 
   //static const char r21_datafile[] = "SDSS_Mr21_rspace.dat"; 
   // Stores the number of points in random set. 
   long int N_rand = FILELENGTHrand;

   FILE *myrandfile = fopen ( random_datafile, "r" );
   if (myrandfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }
   
   double randRA_LIST[N_rand]; 
   double randDEC_LIST[N_rand]; 
   double randZ_LIST[N_rand]; 
   
   //===============
   // READ IN DATA. 
   //===============
   for(i = 0; i < (N_rand); ++i){
      fscanf(myrandfile, "%lf", &randRA_LIST[i]);
      fscanf(myrandfile, "%lf", &randDEC_LIST[i]);
      fscanf(myrandfile, "%lf", &randZ_LIST[i]);
      //printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );
      
      if (i >= (N_rand-2)){
         printf("Close or exceeded N_data limit. RA: %lf \n", randRA_LIST[i]);
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

   /*
   //ALREADY IINITIALIZED FROM ABOVE.
   double D, logD; 
   double r = 0; 
   double DEC_rad = 0;
   double RA_rad = 0;
   double x1 = 0;
   double y1 = 0;
   double z1 = 0;  

   double rj = 0; 
   double DEC_radj = 0;
   double RA_radj = 0;
   double x2 = 0;
   double y2 = 0;
   double z2 = 0;  
   */


   #pragma omp parallel shared( randZ_LIST, randDEC_LIST, randRA_LIST, N_rand, deg_to_rad ,chunk) private (D,\
   logD, r, rj, DEC_rad, DEC_radj, RA_rad, RA_radj, x1, y1, z1, x2, y2, z2, dist_index, i, j )
   {
      

      long int sum_local_counts[15];
      memset(sum_local_counts, 0, 15 * sizeof(sum_local_counts[0]) );


      #pragma omp for
      for(i=0; i< (N_rand-1); ++i){
         
         //printf("Inside first: %lf %lf\n", RA_LIST[i], x1);
         //x1,y1,z1 = XYZ_given_RADECZ(RA_LIST[i], DEC_LIST[i], Z_LIST[i]);
         //XYZ_given_RADECZ(RA_LIST[i], DEC_LIST[i], Z_LIST[i], x1, y1, z1, dec_to_rad);
         //r = radius_given_z(Z_LIST[i]);
         r = 2998. * randZ_LIST[i];
         DEC_rad = randDEC_LIST[i] * (deg_to_rad);
         RA_rad  = randRA_LIST[i]  * (deg_to_rad);

         x1 = r * cos( DEC_rad ) * cos( RA_rad ); 
         y1 = r * cos( DEC_rad ) * sin( RA_rad ); 
         z1 = r * sin( DEC_rad );

         
         for(j=0; j< (N_rand-1); ++j){
            if (j!=i){
               //x2,y2,z2 = XYZ_given_RADECZ(RA_LIST[j], DEC_LIST[j], Z_LIST[j]);
               //XYZ_given_RADECZ(RA_LIST[j], DEC_LIST[j], Z_LIST[j], x2, y2, z2, dec_to_rad);
               rj = 2998. * randZ_LIST[j];
               DEC_radj = randDEC_LIST[j] * (deg_to_rad);
               RA_radj  = randRA_LIST[j]  * (deg_to_rad);

               x2 = rj * cos( DEC_radj ) * cos( RA_radj ); 
               y2 = rj * cos( DEC_radj ) * sin( RA_radj ); 
               z2 = rj * sin( DEC_radj );

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

   //Mr 21.

   printf("Calculating Mr21 Correlation Function...\n");
   printf("Nrand: %ld\n", N_rand);
   printf("Ndata: %ld\n", N_data);
   /*
   for(i=0; i<num_bins; i++){
      //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
      Xi_func[i] = (N_rand / N_data)*(N_rand / N_data) * \
              ((distance_counts[i]/randdistance_counts[i]) - 1.0);

      printf("%f ", Xi_func[i]); 
   }
   */
   
   double ratio = (double) N_rand / (double) N_data;
   for(i=0; i<num_bins; i++){
      //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
      Xi_func[i] = ratio * ratio *  ( (double) distance_counts[i] / (double) randdistance_counts[i]) - 1.0 ;
      printf("%f ", Xi_func[i]); 
   }
   printf("\n");

   //Mr 20.
   printf("Calculating Mr20 Correlation Function...\n");
   printf("Nrand: %ld", N_rand);
   printf("N20data: %ld", N20_data);

   for(i=0; i<num_bins; i++){
      //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
      Xi20_func[i] = (N_rand / N20_data)*(N_rand / N20_data) * \
              ((distance20_counts[i]/randdistance_counts[i]) - 1.0);

      printf("%f ", Xi20_func[i]); 
   }
   printf("\n");
   //Mr 20z
   printf("Calculating Mr20z Correlation Function...\n");
   for(i=0; i<num_bins; i++){
      //Compute the Correlation Function: Xi = (Ngal/Nrand)^2 * (DD/RR  - 1)
      Xi20z_func[i] = (N_rand / N20z_data)*(N_rand / N20z_data) * \
              ((distance20z_counts[i]/randdistance_counts[i]) - 1.0);

      printf("%f ", Xi20z_func[i]); 
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
   // WRITING 21r RESULTS TO FILE
   // ==========================
   printf("Saving Mr21r counts to file.\n");
   FILE *fp_out;
   fp_out = fopen("output_Mr21counts.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", distance_counts[i]);
   } 
   fclose(fp_out);


   // ==========================
   // WRITING 20r RESULTS TO FILE
   // ==========================
   printf("Saving Mr20r counts to file.\n");

   //FILE *fp_out;
   fp_out = fopen("output_Mr20counts.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", distance20_counts[i]);
   } 
   fclose(fp_out);


   // ==========================
   // WRITING 20z RESULTS TO FILE
   // ==========================
   printf("Saving Mr20z counts to file.\n");
   //FILE *fp_out;
   fp_out = fopen("output_Mr20zcounts.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", distance20z_counts[i]);
   } 
   fclose(fp_out);


   // ==========================
   // WRITING random RESULTS TO FILE
   // ==========================
   printf("Saving Random counts to file.\n");
   //FILE *fp_out;
   fp_out = fopen("output_counts_random.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%ld \n", randdistance_counts[i]);
   } 
   fclose(fp_out);



   // ===============================
   // WRITING Mr21 Xi RESULTS TO FILE
   // ===============================
   printf("Saving Xi21 to file.\n");
   //FILE *fp_out;
   fp_out = fopen("output_Xi_21r.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%f \n", Xi_func[i]);
   } 
   fclose(fp_out);

   // ===============================
   // WRITING Mr20 Xi RESULTS TO FILE
   // ===============================
   printf("Saving Xi20 to file.\n");

   //FILE *fp_out;
   fp_out = fopen("output_Xi_20r.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%f \n", Xi20_func[i]);
   } 
   fclose(fp_out);

   // ================================
   // WRITING Mr20z Xi RESULTS TO FILE
   // ================================
   printf("Saving Xi20z to file.\n");

   //FILE *fp_out;
   fp_out = fopen("output_Xi_20z.txt","w");
   if ( fp_out == NULL ) {
      printf("output_file.txt not opened, exiting...\n");
      exit(0);
   }
   for ( i=0 ; i < num_bins ; i++ ) {
      fprintf(fp_out,"%f \n", Xi20z_func[i]);
   } 
   fclose(fp_out);


   /*
   free(RA_LIST);
   free(DEC_LIST);
   free(Z_LIST);
   */


   return 0;

}



//END
