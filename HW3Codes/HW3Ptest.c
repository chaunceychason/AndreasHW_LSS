//#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define CHUNKSIZE 100
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

float radius_given_z(float z_val);

void XYZ_given_RADECZ(float ra, float dec, float z_val, float *x, float *y, float *z, float deg2rad);

float distance_given_2points(float x1,float  y1, float z1, float x2, float y2, float z2 );

int main ()  
{
   printf("TEST - FIRST of MAIN\n");
   long int i, j;
   int chunk;
   //float a[N], b[N], c[N];

   float dec_to_rad = (2. * PI / 360.);
   //Create the Distance Array. 
   float logr_min = -1.0;
   float logr_max = 1.3011;
   long int errorcounts = 0;
   printf("line104\n");
   //    r                 0.1,                 ...                  , 20 
   // logr                 -1,                  ...                  ,1.3011
   // index                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
   long int distance_counts[15] = { 0 };
   long int randdistance_counts[15] = { 0 }; 
   float Xi_func[15] = { 0.0 };
   printf("Line 112\n");

   static const int num_bins = 15;
   int dist_index;

   //This should be set below 5495 to take a limited sample. 
   long int FILELENGTH = 4000;
   long int FILELENGTH20r  = 28161; 
   long int FILELENGTH20z  = 28382; 
   long int FILELENGTH21r  = 5494;  
   long int FILELENGTHrand = 42653; 

   //long int FILELENGTH20z  = 29000; 
   //long int FILELENGTH21r  = 6000;  //increased by 2
   //long int FILELENGTHrand = 45000; //increased by 2

   //string datafilename = 'SDSS_random.dat'
   //static const char random_datafile[] = "SDSS_random.dat"; 

   // =========================================================
   
   
   printf("Line 131\n");
   static const char r20_datafile[] = "SDSS_Mr20_rspace.dat"; 
   FILE *myfile = fopen ( r20_datafile, "r" );
   if (myfile == NULL) {
      printf("input_file.txt not opened, exiting...\n");
      exit(0);
   }

   printf("Opened file - Begining assignment.\n");
   
   long int N_data = FILELENGTH20r;
   N_data = 5000;
   //FILELENGTH = FILELENGTH21r;
   
   float RA_LIST[N_data];
   float DEC_LIST[N_data];
   float Z_LIST[N_data];
   
   //float *Z_LIST;
   //float *DEC_LIST;
   //float *RA_LIST;
   /*
   float * Z_LIST   = malloc(N_data * sizeof(float));
   float * DEC_LIST = malloc(N_data * sizeof(float));
   float * RA_LIST  = malloc(N_data * sizeof(float));
   */

   /*
   float * RA_LIST[N_data] = malloc(N_data * sizeof(float)); 
   float * DEC_LIST[N_data] = malloc(N_data * sizeof(float)); 
   float * Z_LIST[N_data] = malloc(N_data * sizeof(float)); 
   */

   /*
   float *RA;
   float *DEC;
   float *Z;
   */

   float RA;
   float DEC;
   float Z;
   //===============
   // READ IN DATA. 
   //===============
   i = 0;
   //while(fscanf(myfile, "%f %f %f", &RA,&DEC,&Z) != EOF){
   //while(fscanf(myfile, "%f", &RA_LIST[i]) != EOF){
   for(i = 0; i < N_data; i++){
      fscanf(myfile, "%f", &RA_LIST[i]);
      fscanf(myfile, "%f", &DEC_LIST[i]);
      fscanf(myfile, "%f", &Z_LIST[i]);
      printf("i: %ld Data Z: %f \n", i, Z_LIST[i] );

      if (i >= (N_data-5)){
         //printf("Exceeded N_data limit.\n");
         //break;   
      }
      ++i;
   }

   fclose(myfile);
   printf("Closing File.");

   //printf("%f", RA_LIST[N_data-1]);



   /*
   //Begin Nested For Loop to Count Pairs. 
   float *x1 = 0;
   float *x2 = 0;
   float *y1 = 0;
   float *y2 = 0;
   float *z1 = 0;
   float *z2 = 0;

   float D, logD; 
   //======================
   // COMPUTE DATA COUNTS
   //======================
   for(i=0; i< FILELENGTH21r; i++){
      //x1,y1,z1 = XYZ_given_RADECZ(RA_LIST[i], DEC_LIST[i], Z_LIST[i]);
      XYZ_given_RADECZ(RA_LIST[i], DEC_LIST[i], Z_LIST[i], x1, y1, z1, dec_to_rad);

      for(j=0; j< FILELENGTH21r; j++){
         if (j!=i){
            //x2,y2,z2 = XYZ_given_RADECZ(RA_LIST[j], DEC_LIST[j], Z_LIST[j]);
            XYZ_given_RADECZ(RA_LIST[j], DEC_LIST[j], Z_LIST[j], x2, y2, z2, dec_to_rad);
            D = distance_given_2points(*x1, *y1, *z1, *x2, *y2, *z2);
            logD = log10(D);
            //dist_index = (logD +1)*(num_bins/(logr_max +1))
            dist_index = floor(logD - logr_min)*(num_bins/(logr_max - logr_min));            
            if (dist_index >= 0 && dist_index < num_bins){
               //Increment the appropiate bin.
               if (dist_index > 14)
                  printf("YELLING!");
               distance_counts[dist_index] += 1;
            } 
         }
      //end inner for loop
      }
   //end outer for loop
   }

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   // =========================================================
   // INIITIALIZE FOR RANDOM COUNTS


   float *randZ_LIST;
   float *randDEC_LIST;
   float *randRA_LIST;

   randZ_LIST   = (float *) calloc(N_rand, sizeof(float));
   randDEC_LIST = (float *) calloc(N_rand, sizeof(float));
   randRA_LIST  = (float *) calloc(N_rand, sizeof(float));

   
   float * randRA_LIST[N_rand] = malloc(N_data * sizeof(float)); 
   float * randDEC_LIST[N_rand] = malloc(N_data * sizeof(float)); 
   float * randZ_LIST[N_rand] = malloc(N_data * sizeof(float)); 
   */


   /*
   free(RA_LIST);
   free(DEC_LIST);
   free(Z_LIST);

   
   free(randRA_LIST);
   free(randDEC_LIST);
   free(randZ_LIST);
   */




   return 0;

}

float radius_given_z(float z_val)
{
   float radius = 2950. * z_val;
   return radius;
}

void XYZ_given_RADECZ(float ra, float dec, float z_val, float *xloc, float *yloc, float *zloc, float deg2rad)
{
   float r = radius_given_z(z_val);
   float DEC_rad = dec * (deg2rad);
   float RA_rad  = ra * (deg2rad);

   *xloc = r * cos( DEC_rad ) * cos( RA_rad ); 
   *yloc = r * cos( DEC_rad ) * sin( RA_rad ); 
   *zloc = r * sin( DEC_rad );
}

float distance_given_2points(float x1,float  y1, float z1, float x2, float y2, float z2 )
{
   return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
}




//END
