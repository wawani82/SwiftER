/* *
 * Copyright (c) 2014, James S. Plank and Kevin Greenan
 * All rights reserved.
 *
 * Jerasure - A C/C++ Library for a Variety of Reed-Solomon and RAID-6 Erasure
 * Coding Techniques
 *
 * Revision 2.0: Galois Field backend now links to GF-Complete
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  - Neither the name of the University of Tennessee nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Jerasure's authors:

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan.
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank.
 */

/* 

This program takes as input an inputfile, k, m, a coding 
technique, w, and packetsize.  It creates k+m files from 
the original file so that k of these files are parts of 
the original file and m of the files are encoded based on 
the given coding technique. The format of the created files 
is the file name with "_k#" or "_m#" and then the extension.  
(For example, inputfile test.txt would yield file "test_k1.txt".)
*/

#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <gf_rand.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "cauchy.h"
#include "liberation.h"

//whcho added
#include "galois.h" 
#include <math.h>
#include "timing.h" 


#define N 10



enum Coding_Technique {Reed_Sol_Van, Reed_Sol_R6_Op, Cauchy_Orig, Cauchy_Good, Liberation, Blaum_Roth, Liber8tion, RDP, EVENODD, No_Coding};

char *Methods[N] = {"reed_sol_van", "reed_sol_r6_op", "cauchy_orig", "cauchy_good", "liberation", "blaum_roth", "liber8tion", "no_coding"};

/* Global variables for signal handler */
int readins, n;
enum Coding_Technique method;

/* Function prototypes */
int is_prime(int w);
void ctrl_bs_handler(int dummy);

int jfread(void *ptr, int size, int nmembers, FILE *stream)
{
  int nd;
  int *li, i;
  if (stream != NULL) return fread(ptr, size, nmembers, stream);

  MOA_Fill_Random_Region(ptr, size);
  return size;
}


int main (int argc, char **argv) {
	FILE *fp, *fp2;				// file pointers
	char *memblock;				// reading in file
	char *block;				// padding file
	int size, newsize;			// size of file and temp size 
	struct stat status;			// finding file size

	
	enum Coding_Technique tech;		// coding technique (parameter)
	int k, m, w, packetsize;		// parameters
	int buffersize;					// paramter
	int i, j;						// loop control variables
	int blocksize;					// size of k+m files
	int total;
	int extra;
	
	/* Jerasure Arguments */
	char **data;				
	char **coding;
	int *matrix;
	int *bitmatrix;
	int **schedule;
	int *erasure;
	int *erased;
	
	/* Creation of file name variables */
	char temp[5];
	char *s1, *s2, *extension;
	char *fname;
	int md;
	char *curdir;
	
	/* Timing variables */
	struct timeval t1, t2, t3, t4;
	struct timezone tz;
	double tsec;
	double totalsec;
	struct timeval start, stop;

	/* Find buffersize */
	int up, down;


//whcho added
int integer;
struct timeval t_bus_start, t_bus_end;
struct timeval t_read_start, t_read_end;
struct timeval t_write_start, t_write_end;
struct timeval t_cal_start, t_cal_end;
struct timeval t_whcho_start, t_whcho_end;
double t_total_time;
double t_total_read ;
double t_total_write ;
double t_total_cal ;
struct timing t_en_read_start, t_en_read_end ;
struct timing t_en_write_start, t_en_write_end ;
double t_total_en_read ;
double t_total_en_write ;
double t_whcho ;



	signal(SIGQUIT, ctrl_bs_handler);

	/* Start timing */
	gettimeofday(&t1, &tz);
	totalsec = 0.0;
	matrix = NULL;
	bitmatrix = NULL;
	schedule = NULL;
	
	/* Error check Arguments*/
	if (argc != 8) {
		fprintf(stderr,  "usage: inputfile k m coding_technique w packetsize buffersize\n");
		fprintf(stderr,  "\nChoose one of the following coding techniques: \nreed_sol_van, \nreed_sol_r6_op, \ncauchy_orig, \ncauchy_good, \nliberation, \nblaum_roth, \nliber8tion");
		fprintf(stderr,  "\n\nPacketsize is ignored for the reed_sol's");
		fprintf(stderr,  "\nBuffersize of 0 means the buffersize is chosen automatically.\n");
		fprintf(stderr,  "\nIf you just want to test speed, use an inputfile of \"-number\" where number is the size of the fake file you want to test.\n\n");
		exit(0);
	}
	/* Conversion of parameters and error checking */	
	if (sscanf(argv[2], "%d", &k) == 0 || k <= 0) {
		fprintf(stderr,  "Invalid value for k\n");
		exit(0);
	}
	if (sscanf(argv[3], "%d", &m) == 0 || m < 0) {
		fprintf(stderr,  "Invalid value for m\n");
		exit(0);
	}
	if (sscanf(argv[5],"%d", &w) == 0 || w <= 0) {
		fprintf(stderr,  "Invalid value for w.\n");
		exit(0);
	}
	if (argc == 6) {
		packetsize = 0;
	}
	else {
		if (sscanf(argv[6], "%d", &packetsize) == 0 || packetsize < 0) {
			fprintf(stderr,  "Invalid value for packetsize.\n");
			exit(0);
		}
	}
	if (argc != 8) {
		buffersize = 0;
	}
	else {
		if (sscanf(argv[7], "%d", &buffersize) == 0 || buffersize < 0) {
			fprintf(stderr, "Invalid value for buffersize\n");
			exit(0);
		}
		
	}


//whcho added
/******************** for regenerating code, 3*k and 3*m **************/

	//k=2*k;
	//m=2*m;
	k=3*k;
	m=3*m;



	/* Determine proper buffersize by finding the closest valid buffersize to the input value  */
	if (buffersize != 0) {
		if (packetsize != 0 && buffersize%(sizeof(long)*w*k*packetsize) != 0) { 
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k*packetsize) != 0 && (down%(sizeof(long)*w*k*packetsize) != 0)) {
				up++;
				if (down == 0) {
					down--;
				}
			}
			if (up%(sizeof(long)*w*k*packetsize) == 0) {
				buffersize = up;
			}
			else {
				if (down != 0) {
					buffersize = down;
				}
			}
		}
		else if (packetsize == 0 && buffersize%(sizeof(long)*w*k) != 0) {
			up = buffersize;
			down = buffersize;
			while (up%(sizeof(long)*w*k) != 0 && down%(sizeof(long)*w*k) != 0) {
				up++;
				down--;
			}
			if (up%(sizeof(long)*w*k) == 0) {
				buffersize = up;
			}
			else {
				buffersize = down;
			}
		}
	}

	/* Setting of coding technique and error checking */
	
	if (strcmp(argv[4], "no_coding") == 0) {
		tech = No_Coding;
	}
	else if (strcmp(argv[4], "reed_sol_van") == 0) {
		tech = Reed_Sol_Van;
		if (w != 8 && w != 16 && w != 32) {
			fprintf(stderr,  "w must be one of {8, 16, 32}\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "reed_sol_r6_op") == 0) {
		if (m != 2) {
			fprintf(stderr,  "m must be equal to 2\n");
			exit(0);
		}
		if (w != 8 && w != 16 && w != 32) {
			fprintf(stderr,  "w must be one of {8, 16, 32}\n");
			exit(0);
		}
		tech = Reed_Sol_R6_Op;
	}
	else if (strcmp(argv[4], "cauchy_orig") == 0) {
		tech = Cauchy_Orig;
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "cauchy_good") == 0) {
		tech = Cauchy_Good;
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
	}
	else if (strcmp(argv[4], "liberation") == 0) {
		if (k > w) {
			fprintf(stderr,  "k must be less than or equal to w\n");
			exit(0);
		}
		if (w <= 2 || !(w%2) || !is_prime(w)) {
			fprintf(stderr,  "w must be greater than two and w must be prime\n");
			exit(0);
		}
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
		if ((packetsize%(sizeof(long))) != 0) {
			fprintf(stderr,  "packetsize must be a multiple of sizeof(long)\n");
			exit(0);
		}
		tech = Liberation;
	}
	else if (strcmp(argv[4], "blaum_roth") == 0) {
		if (k > w) {
			fprintf(stderr,  "k must be less than or equal to w\n");
			exit(0);
		}
		if (w <= 2 || !((w+1)%2) || !is_prime(w+1)) {
			fprintf(stderr,  "w must be greater than two and w+1 must be prime\n");
			exit(0);
		}
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize.\n");
			exit(0);
		}
		if ((packetsize%(sizeof(long))) != 0) {
			fprintf(stderr,  "packetsize must be a multiple of sizeof(long)\n");
			exit(0);
		}
		tech = Blaum_Roth;
	}
	else if (strcmp(argv[4], "liber8tion") == 0) {
		if (packetsize == 0) {
			fprintf(stderr, "Must include packetsize\n");
			exit(0);
		}
		if (w != 8) {
			fprintf(stderr, "w must equal 8\n");
			exit(0);
		}
		if (m != 2) {
			fprintf(stderr, "m must equal 2\n");
			exit(0);
		}
		if (k > w) {
			fprintf(stderr, "k must be less than or equal to w\n");
			exit(0);
		}
		tech = Liber8tion;
	}
	else {
		fprintf(stderr,  "Not a valid coding technique. Choose one of the following: reed_sol_van, reed_sol_r6_op, cauchy_orig, cauchy_good, liberation, blaum_roth, liber8tion, no_coding\n");
		exit(0);
	}

	/* Set global variable method for signal handler */
	method = tech;

	/* Get current working directory for construction of file names */
	curdir = (char*)malloc(sizeof(char)*1000);	
	getcwd(curdir, 1000);

        if (argv[1][0] != '-') {

		/* Open file and error check */
		fp = fopen(argv[1], "rb");
		if (fp == NULL) {
			fprintf(stderr,  "Unable to open file.\n");
			exit(0);
		}
	
		/* Create Coding directory */
		i = mkdir("Coding", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Coding directory.\n");
			exit(0);
		}
		
//whcho added 
//Create NewNode directory !
		i = mkdir("Coding/NewNode1", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create NewNode1 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/NewNode2", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create NewNode2 directory.\n");
			exit(0);
		}


// Create Node1 ~ Node10 directory !!
	
		i = mkdir("Coding/Node1", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node1 directory.\n");
			exit(0);
		}
	
		i = mkdir("Coding/Node2", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node2 directory.\n");
			exit(0);
		}

		i = mkdir("Coding/Node3", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node3 directory.\n");
			exit(0);
		}

		i = mkdir("Coding/Node4", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node4 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node5", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node5 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node6", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node6 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node7", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node7 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node8", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node8 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node9", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node9 directory.\n");
			exit(0);
		}
		i = mkdir("Coding/Node10", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node10 directory.\n");
			exit(0);
		}


// whcho added for test !!!(150611)


		i = mkdir("Coding/Node11", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node10 directory.\n");
			exit(0);
		}

		i = mkdir("Coding/Node12", S_IRWXU);
		if (i == -1 && errno != EEXIST) {
			fprintf(stderr, "Unable to create Node10 directory.\n");
			exit(0);
		}





		/* Determine original size of file */
		stat(argv[1], &status);	
		size = status.st_size;
        } else {
        	if (sscanf(argv[1]+1, "%d", &size) != 1 || size <= 0) {
                	fprintf(stderr, "Files starting with '-' should be sizes for randomly created input\n");
			exit(1);
		}
        	fp = NULL;
		MOA_Seed(time(0));
        }

	newsize = size;
	
	/* Find new size by determining next closest multiple */
	if (packetsize != 0) {
		if (size%(k*w*packetsize*sizeof(long)) != 0) {
			while (newsize%(k*w*packetsize*sizeof(long)) != 0) 
				newsize++;
		}
	}
	else {
		if (size%(k*w*sizeof(long)) != 0) {
			while (newsize%(k*w*sizeof(long)) != 0) 
				newsize++;
		}
	}
	
	if (buffersize != 0) {
		while (newsize%buffersize != 0) {
			newsize++;
		}
	}


	/* Determine size of k+m files */
	blocksize = newsize/k;
//whcho added
	printf("blocksize=%d\n",blocksize);
	printf("buffersize=%d\n",buffersize);
	printf("size=%d\n",size);

	/* Allow for buffersize and determine number of read-ins */
	if (size > buffersize && buffersize != 0) {
		if (newsize%buffersize != 0) {
			readins = newsize/buffersize;
		}
		else {
			readins = newsize/buffersize;
		}
		block = (char *)malloc(sizeof(char)*buffersize);
		blocksize = buffersize/k;
	}
	else {
		readins = 1;
		buffersize = size;
		block = (char *)malloc(sizeof(char)*newsize);
	}

//whcho added
	printf("blocksize=%d\n",blocksize);
	printf("buffersize=%d\n",buffersize);
	printf("size=%d\n",size);


	/* Break inputfile name into the filename and extension */	
	s1 = (char*)malloc(sizeof(char)*(strlen(argv[1])+20));
	s2 = strrchr(argv[1], '/');
	if (s2 != NULL) {
		s2++;
		strcpy(s1, s2);
	}
	else {
		strcpy(s1, argv[1]);
	}
	s2 = strchr(s1, '.');
	if (s2 != NULL) {
          extension = strdup(s2);
          *s2 = '\0';
	} else {
          extension = strdup("");
        }
	
	/* Allocate for full file name */
	fname = (char*)malloc(sizeof(char)*(strlen(argv[1])+strlen(curdir)+20));
	sprintf(temp, "%d", k);
	md = strlen(temp);
	
	/* Allocate data and coding */
	data = (char **)malloc(sizeof(char*)*k);
	coding = (char **)malloc(sizeof(char*)*m);
	for (i = 0; i < m; i++) {
		coding[i] = (char *)malloc(sizeof(char)*blocksize);
                if (coding[i] == NULL) { perror("malloc"); exit(1); }
	}

	

	/* Create coding matrix or bitmatrix and schedule */
	gettimeofday(&t3, &tz);
	switch(tech) {
		case No_Coding:
			break;
		case Reed_Sol_Van:
			matrix = reed_sol_vandermonde_coding_matrix(k, m, w);
			
			//whcho added

			jerasure_print_matrix(matrix, m, k, w);
			printf("\n\n\n");
			break;

		case Cauchy_Orig:
			matrix = cauchy_original_coding_matrix(k, m, w);
			bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, matrix);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Cauchy_Good:
			matrix = cauchy_good_general_coding_matrix(k, m, w);
			bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, matrix);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;	
		case Liberation:
			bitmatrix = liberation_coding_bitmatrix(k, w);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Blaum_Roth:
			bitmatrix = blaum_roth_coding_bitmatrix(k, w);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
		case Liber8tion:
			bitmatrix = liber8tion_coding_bitmatrix(k);
			schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
			break;
	}
	gettimeofday(&start, &tz);	
	gettimeofday(&t4, &tz);
	tsec = 0.0;
	tsec += t4.tv_usec;
	tsec -= t3.tv_usec;
	tsec /= 1000000.0;
	tsec += t4.tv_sec;
	tsec -= t3.tv_sec;
	totalsec += tsec;

//whcho added
t_total_en_read=0.0;
t_total_en_write=0.0;


	

	/* Read in data until finished */
	n = 1;
	total = 0;
//whcho added
	printf("readins=%d\n",readins);


	while (n <= readins) {

//whcho added
timing_set(&t_en_read_start);
		

		/* Check if padding is needed, if so, add appropriate 
		   number of zeros */
		if (total < size && total+buffersize <= size) {
			total += jfread(block, sizeof(char), buffersize, fp);
		}
		else if (total < size && total+buffersize > size) {
			extra = jfread(block, sizeof(char), buffersize, fp);
			for (i = extra; i < buffersize; i++) {
				block[i] = '0';
			}
		}
		else if (total == size) {
			for (i = 0; i < buffersize; i++) {
				block[i] = '0';
			}
		}
	
		/* Set pointers to point to file data */
		for (i = 0; i < k; i++) {
			data[i] = block+(i*blocksize);
		}

//whcho added
timing_set(&t_en_read_end);



 	 	gettimeofday(&t3, &tz);
		/* Encode according to coding method */
		switch(tech) {	
			case No_Coding:
				break;
			case Reed_Sol_Van:
				jerasure_matrix_encode(k, m, w, matrix, data, coding, blocksize);
				break;
			case Reed_Sol_R6_Op:
				reed_sol_r6_encode(k, w, data, coding, blocksize);
				break;
			case Cauchy_Orig:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Cauchy_Good:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Liberation:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Blaum_Roth:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
			case Liber8tion:
				jerasure_schedule_encode(k, m, w, schedule, data, coding, blocksize, packetsize);
				break;
		}
		gettimeofday(&t4, &tz);


//whcho added
timing_set(&t_en_write_start);

	
		/* Write data and encoded data to k+m files */
		for	(i = 1; i <= k; i++) {
			if (fp == NULL) {
				bzero(data[i-1], blocksize);
 			} else {
//whcho added
/* /mnt/node1 ~ /mnt/node8  */
/* Now,  3 data fragments in 1 node */


//integer=(i+1)/2;
integer=(i+2)/3;
printf("integer= %d\n",integer);


				//sprintf(fname, "%s/Coding/Node%d/%s_k%0*d%s", curdir, integer , s1, md, i, extension);
				sprintf(fname, "/mnt/node%d/%s_k%0*d%s",  integer    , s1, md, i, extension);
				if (n == 1) {
					fp2 = fopen(fname, "wb");
				}
				else {
					fp2 = fopen(fname, "ab");
				}
				fwrite(data[i-1], sizeof(char), blocksize, fp2);
				fclose(fp2);
			}
			
		}

		for	(i = 1; i <= m; i++) {
			if (fp == NULL) {
				bzero(data[i-1], blocksize);
 			} else {
//whcho added
/* /mnt/node9 ~ /mnt/node10 for parity */
				
/* Now,  3 data fragments in 1 node */

//integer=(k+1+i)/2;
integer=(k+2+i)/3;
printf("integer= %d\n",integer);

				//sprintf(fname, "%s/Coding/Node%d/%s_m%0*d%s", curdir, integer , s1, md, i, extension);
				sprintf(fname, "/mnt/node%d/%s_m%0*d%s", integer , s1, md, i, extension);
				if (n == 1) {
					fp2 = fopen(fname, "wb");
				}
				else {
					fp2 = fopen(fname, "ab");
				}
				fwrite(coding[i-1], sizeof(char), blocksize, fp2);
				fclose(fp2);
			}
		
		
		}
		n++;


//whcho added
timing_set(&t_en_write_end);



//whcho added
		//printf("n=%d\n",n);
		//printf("total=%d\n",total);
		printf("while() performed %d times \n",n-1);

		/* Calculate encoding time */
		tsec = 0.0;
		tsec += t4.tv_usec;
		tsec -= t3.tv_usec;
		tsec /= 1000000.0;
		tsec += t4.tv_sec;
		tsec -= t3.tv_sec;
		totalsec += tsec;
	}

	/* Create metadata file */
        if (fp != NULL) {
		sprintf(fname, "%s/Coding/%s_meta.txt", curdir, s1);
		fp2 = fopen(fname, "wb");
		fprintf(fp2, "%s\n", argv[1]);
		fprintf(fp2, "%d\n", size);
		fprintf(fp2, "%d %d %d %d %d\n", k, m, w, packetsize, buffersize);
		fprintf(fp2, "%s\n", argv[4]);
		fprintf(fp2, "%d\n", tech);
		fprintf(fp2, "%d\n", readins);
		fclose(fp2);
	}


	/* Free allocated memory */
	//whcho move to the end 
	//free(s1);
	//free(fname);
	//free(block);
	//free(curdir);
	
	/* Calculate rate in MB/sec and print */
	gettimeofday(&t2, &tz);
	tsec = 0.0;
	tsec += t2.tv_usec;
	tsec -= t1.tv_usec;
	tsec /= 1000000.0;
	tsec += t2.tv_sec;
	tsec -= t1.tv_sec;
	printf("Encoding (MB/sec): %0.10f\n", (((double) size)/1024.0/1024.0)/totalsec);
	printf("En_Total (MB/sec): %0.10f\n", (((double) size)/1024.0/1024.0)/tsec);
//whcho added
	printf("Encoding Time (sec): %0.5f\n", totalsec);
	printf("En_Total Time (sec): %0.5f\n", tsec);


//whcho added
	t_total_en_read += timing_delta(&t_en_read_start, &t_en_read_end);
	t_total_en_write += timing_delta(&t_en_write_start, &t_en_write_end);

	printf("Encoding Read Time (sec): %0.6f\n", t_total_en_read);
	printf("Encoding Write Time (sec): %0.6f\n", t_total_en_write);



//whcho added 
/*********************** Multipying factor and XOR     ***********************/
///******************* 1. k1, k2 file open and open parity file *************************************/
///******************* 2. allocate memory for k1, k2 */

	
/* Start timing */
struct timeval t5,t6;
gettimeofday(&t5, &tz);
tsec = 0.0;
	
	
char fname1[100];
char fname2[100];
char fname3[100];
char p_fname[100];
char *block1;
char *block2;
char **ori_data;
char **p_coding;
int counter=0;
FILE *f_p1, *f_p2 , *f_p3, *p_fp ;

int transfer_counter=0;
FILE *src_fp , *dst_fp ;
char t_fname[100];



/* for XOR parity 8 times at the end */
char *block3;
char *block4;
char *block5;
char *block6;
char *block7;
char *block8;
char **par_data;

int index;

//int buf[4][k];
int buf[6][k];

int buf_index;
int z;
//int integer;


char *trans_block1;
char *trans_block2;
char *trans_block3;
char *trans_block4;
char *trans_block5;
char *trans_block6;
char *trans_block7;
char *trans_block8;



block3 = (char *)malloc(sizeof(char)*blocksize);
block4 = (char *)malloc(sizeof(char)*blocksize);
block5 = (char *)malloc(sizeof(char)*blocksize);
block6 = (char *)malloc(sizeof(char)*blocksize);
block7 = (char *)malloc(sizeof(char)*blocksize);
block8 = (char *)malloc(sizeof(char)*blocksize);


//ori_data = (char **)malloc(sizeof(char*)*2);
ori_data = (char **)malloc(sizeof(char*)*3);
block1 = (char *)malloc(sizeof(char)*blocksize);
block2 = (char *)malloc(sizeof(char)*blocksize);


p_coding = (char **)malloc(sizeof(char*)*8);


trans_block1 = (char *)malloc(sizeof(char)*blocksize);
trans_block2 = (char *)malloc(sizeof(char)*blocksize);
trans_block3 = (char *)malloc(sizeof(char)*blocksize);
trans_block4 = (char *)malloc(sizeof(char)*blocksize);
trans_block5 = (char *)malloc(sizeof(char)*blocksize);
trans_block6 = (char *)malloc(sizeof(char)*blocksize);
trans_block7 = (char *)malloc(sizeof(char)*blocksize);
trans_block8 = (char *)malloc(sizeof(char)*blocksize);


char *dptr , *sptr; 



/* Now, coefficients for 6 new parity m07 m08 m09 m10 m11 m12   */
/* m07 parity */
buf[0][0]=30;
buf[0][1]=48;
buf[0][2]=30;

buf[0][3]=28;
buf[0][4]=29;
buf[0][5]=42;

buf[0][6]=52;
buf[0][7]=45;
buf[0][8]=58;

buf[0][9]=74;
buf[0][10]=61;
buf[0][11]=61;

buf[0][12]=52;
buf[0][13]=55;
buf[0][14]=48;

buf[0][15]=57;
buf[0][16]=54;
buf[0][17]=46;

buf[0][18]=28;
buf[0][19]=38;
buf[0][20]=37;

buf[0][21]=51;
buf[0][22]=51;
buf[0][23]=60;

/* m08 parity */
buf[1][0]=3;
buf[1][1]=5;
buf[1][2]=5;

buf[1][3]=13;
buf[1][4]=11;
buf[1][5]=17;

buf[1][6]=13;
buf[1][7]=11;
buf[1][8]=13;

buf[1][9]=26;
buf[1][10]=16;
buf[1][11]=17;

buf[1][12]=17;
buf[1][13]=18;
buf[1][14]=14;

buf[1][15]=13;
buf[1][16]=9;
buf[1][17]=15;

buf[1][18]=15;
buf[1][19]=14;
buf[1][20]=19;

buf[1][21]=19;
buf[1][22]=19;
buf[1][23]=25;


/* m09 parity */
buf[2][0]=27;
buf[2][1]=43;
buf[2][2]=28;

buf[2][3]=31;
buf[2][4]=31;
buf[2][5]=45;

buf[2][6]=58;
buf[2][7]=57;
buf[2][8]=67;

buf[2][9]=65;
buf[2][10]=55;
buf[2][11]=55;

buf[2][12]=62;
buf[2][13]=75;
buf[2][14]=63;

buf[2][15]=63;
buf[2][16]=58;
buf[2][17]=54;

buf[2][18]=36;
buf[2][19]=42;
buf[2][20]=41;

buf[2][21]=55;
buf[2][22]=57;
buf[2][23]=66;





/* m10 parity */
buf[3][0]=24;
buf[3][1]=37;
buf[3][2]=34;

buf[3][3]=31;
buf[3][4]=35;
buf[3][5]=43;

buf[3][6]=44;
buf[3][7]=29;
buf[3][8]=46;

buf[3][9]=65;
buf[3][10]=53;
buf[3][11]=62;

buf[3][12]=58;
buf[3][13]=70;
buf[3][14]=57;

buf[3][15]=58;
buf[3][16]=52;
buf[3][17]=42;

buf[3][18]=36;
buf[3][19]=54;
buf[3][20]=45;

buf[3][21]=56;
buf[3][22]=48;
buf[3][23]=66;




/* m11 parity */
buf[4][0]=26;
buf[4][1]=40;
buf[4][2]=38;

buf[4][3]=40;
buf[4][4]=47;
buf[4][5]=49;

buf[4][6]=53;
buf[4][7]=41;
buf[4][8]=61;

buf[4][9]=61;
buf[4][10]=51;
buf[4][11]=57;

buf[4][12]=60;
buf[4][13]=75;
buf[4][14]=60;

buf[4][15]=68;
buf[4][16]=56;
buf[4][17]=50;

buf[4][18]=28;
buf[4][19]=44;
buf[4][20]=39;

buf[4][21]=51;
buf[4][22]=44;
buf[4][23]=59;



/* m12 parity */
buf[5][0]=15;
buf[5][1]=23;
buf[5][2]=17;

buf[5][3]=30;
buf[5][4]=34;
buf[5][5]=37;

buf[5][6]=41;
buf[5][7]=34;
buf[5][8]=49;

buf[5][9]=62;
buf[5][10]=51;
buf[5][11]=65;

buf[5][12]=43;
buf[5][13]=54;
buf[5][14]=64;

buf[5][15]=63;
buf[5][16]=51;
buf[5][17]=43;

buf[5][18]=29;
buf[5][19]=41;
buf[5][20]=29;

buf[5][21]=50;
buf[5][22]=45;
buf[5][23]=60;








buf_index=0;





/******************************* start line for Video_m05.mp4, m06 , m07 , m08.mp4 **********************************************************************/





z=0;

integer=0;

tsec = 0.0;
t_total_time =0.0;
t_total_read =0.0;
t_total_write =0.0;
t_total_cal =0.0;

t_whcho =0.0;


//for(z=0;z<4;z++){ /* 4 times loop because there is 4 parity files at /mnt/NewNode1 , /mnt/NewNode2 */
for(z=0;z<6;z++){ /* Now, 6 times loop because there is 6 parity files at /mnt/NewNode1 , /mnt/NewNode2 */

//for(z=0;z<1;z++){

//integer=(z+2)/2;
integer=(z+3)/3;
//for debugging 
//printf("y= %d\n",integer);


// 1. K01 & K02 & k03  

//sprintf(fname1, "%s/Coding/Node1/%s_k01.mp4",curdir,s1);
sprintf(fname1, "/mnt/node1/%s_k01.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node1/%s_k02.mp4",curdir,s1);
sprintf(fname2, "/mnt/node1/%s_k02.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node1/%s_k03.mp4",curdir,s1);
sprintf(fname2, "/mnt/node1/%s_k03.mp4",s1);
f_p3 = fopen(fname3, "rb");


gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
//fclose(f_p1); // causing Segmentation Fault !!

counter += jfread(block2, sizeof(char), blocksize, f_p2);
//fclose(f_p2); // causing Segmentation Fault !!

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);

tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


// do measurement of bus transfer !!
gettimeofday(&t_bus_start, &tz);


ori_data[0] = block1 ;
ori_data[1] = block2 ;
ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}


gettimeofday(&t_bus_end, &tz);
tsec = 0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);



//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_01.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_01_%d.mp4",integer, s1 , z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_01_%d.mp4",integer+10, s1 , z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], buf[z][0] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][1] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][2] , blocksize, p_coding[0], 1);



gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;


fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);


fclose(f_p1); // causing Segmentation Fault !!

fclose(f_p2); // causing Segmentation Fault !!



// 2. K04 & K05 & k06

//sprintf(fname1, "%s/Coding/Node2/%s_k04.mp4",curdir,s1);
sprintf(fname1, "/mnt/node2/%s_k04.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node2/%s_k05.mp4",curdir,s1);
sprintf(fname2, "/mnt/node2/%s_k05.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node2/%s_k06.mp4",curdir,s1);
sprintf(fname3, "/mnt/node2/%s_k06.mp4",s1);
f_p3 = fopen(fname3, "rb");


gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);

counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;



gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;
ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}

gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);



//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_02.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_02.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_02_%d.mp4", integer ,s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_02_%d.mp4", integer+10 ,s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], buf[z][3] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][4] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");

galois_w08_region_multiply(ori_data[2], buf[z][5] , blocksize, p_coding[0], 1);


gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);





// 3. k07 , k08 , k09  

//sprintf(fname1, "%s/Coding/Node3/%s_k07.mp4",curdir,s1);
sprintf(fname1, "/mnt/node3/%s_k07.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node3/%s_k08.mp4",curdir,s1);
sprintf(fname2, "/mnt/node3/%s_k08.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node3/%s_k09.mp4",curdir,s1);
sprintf(fname3, "/mnt/node3/%s_k09.mp4",s1);
f_p3 = fopen(fname3, "rb");



gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);

gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;




gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;
ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}
gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_03.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_03.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_03_%d.mp4",integer, s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_03_%d.mp4",integer+10, s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], buf[z][6] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][7] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");

galois_w08_region_multiply(ori_data[2], buf[z][8] , blocksize, p_coding[0], 1);

gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);




// 4. k10 , k11 , k12  

//sprintf(fname1, "%s/Coding/Node4/%s_k10.mp4",curdir,s1);
sprintf(fname1, "/mnt/node4/%s_k10.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node4/%s_k11.mp4",curdir,s1);
sprintf(fname2, "/mnt/node4/%s_k11.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node4/%s_k12.mp4",curdir,s1);
sprintf(fname3, "/mnt/node4/%s_k12.mp4",s1);
f_p3 = fopen(fname3, "rb");

gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);

gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;




gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;

ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}

gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_04.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_04.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_04_%d.mp4",integer, s1 , z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_04_%d.mp4",integer+10, s1 , z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], buf[z][9] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][10] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][11] , blocksize, p_coding[0], 1);

gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);



// 5. k13 , k14 , k15  

//sprintf(fname1, "%s/Coding/Node5/%s_k13.mp4",curdir,s1);
sprintf(fname1, "/mnt/node5/%s_k13.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node5/%s_k14.mp4",curdir,s1);
sprintf(fname2, "/mnt/node5/%s_k14.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node5/%s_k15.mp4",curdir,s1);
sprintf(fname3, "/mnt/node5/%s_k15.mp4",s1);
f_p3 = fopen(fname3, "rb");




gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;

ori_data[2] = block3 ;


for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}
gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_05.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_05.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_05_%d.mp4",integer, s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_05_%d.mp4",integer+10, s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);
galois_w08_region_multiply(ori_data[0], buf[z][12] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][13] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][14] , blocksize, p_coding[0], 1);


gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);




// 6. k16 , k17 , k18

//sprintf(fname1, "%s/Coding/Node6/%s_k16.mp4",curdir,s1);
sprintf(fname1, "/mnt/node6/%s_k16.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node6/%s_k17.mp4",curdir,s1);
sprintf(fname2, "/mnt/node6/%s_k17.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node6/%s_k18.mp4",curdir,s1);
sprintf(fname3, "/mnt/node6/%s_k18.mp4",s1);
f_p3 = fopen(fname3, "rb");

gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;

gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;

ori_data[2] = block3 ;


for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}

gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_06.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_06.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_06_%d.mp4",integer, s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_06_%d.mp4",integer+10, s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);
galois_w08_region_multiply(ori_data[0], buf[z][15] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][16] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][17] , blocksize, p_coding[0], 1);

gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);



//7. k19 , k20 , k21

//sprintf(fname1, "%s/Coding/Node7/%s_k19.mp4",curdir,s1);
sprintf(fname1, "/mnt/node7/%s_k19.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node7/%s_k20.mp4",curdir,s1);
sprintf(fname2, "/mnt/node7/%s_k20.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node7/%s_k21.mp4",curdir,s1);
sprintf(fname3, "/mnt/node7/%s_k21.mp4",s1);
f_p3 = fopen(fname3, "rb");


gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;

ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}

gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_07.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_07.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_07_%d.mp4",integer, s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_07_%d.mp4",integer+10, s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);
galois_w08_region_multiply(ori_data[0], buf[z][18] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][19] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][20] , blocksize, p_coding[0], 1);


gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);




// 8. k22 , k23 , k24


//sprintf(fname1, "%s/Coding/Node8/%s_k22.mp4",curdir,s1);
sprintf(fname1, "/mnt/node8/%s_k22.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node8/%s_k23.mp4",curdir,s1);
sprintf(fname2, "/mnt/node8/%s_k23.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node8/%s_k24.mp4",curdir,s1);
sprintf(fname3, "/mnt/node8/%s_k24.mp4",s1);
f_p3 = fopen(fname3, "rb");



gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
counter += jfread(block2, sizeof(char), blocksize, f_p2);

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);
tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


gettimeofday(&t_bus_start, &tz);

ori_data[0] = block1 ;
ori_data[1] = block2 ;

ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}

gettimeofday(&t_bus_end, &tz);
tsec=0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);

//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/%s_parity_08.mp4",curdir,s1);
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_08.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_08_%d.mp4",integer, s1, z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_08_%d.mp4",integer+10, s1, z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);
galois_w08_region_multiply(ori_data[0], buf[z][21] , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], buf[z][22] , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], buf[z][23] , blocksize, p_coding[0], 1);



gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;

fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);





// 9. m01 & m02 & m03  

//sprintf(fname1, "%s/Coding/Node9/%s_m01.mp4",curdir,s1);
sprintf(fname1, "/mnt/node9/%s_m01.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node9/%s_m02.mp4",curdir,s1);
sprintf(fname2, "/mnt/node9/%s_m02.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node9/%s_m03.mp4",curdir,s1);
sprintf(fname2, "/mnt/node9/%s_m03.mp4",s1);
f_p3 = fopen(fname3, "rb");


gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
//fclose(f_p1); // causing Segmentation Fault !!

counter += jfread(block2, sizeof(char), blocksize, f_p2);
//fclose(f_p2); // causing Segmentation Fault !!

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);

tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


// do measurement of bus transfer !!
gettimeofday(&t_bus_start, &tz);


ori_data[0] = block1 ;
ori_data[1] = block2 ;
ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}


gettimeofday(&t_bus_end, &tz);
tsec = 0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);



//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_09.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_09_%d.mp4",integer, s1 , z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_09_%d.mp4",integer+10, s1 , z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], 1 , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], 1 , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], 1 , blocksize, p_coding[0], 1);



gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;


fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);


fclose(f_p1); // causing Segmentation Fault !!

fclose(f_p2); // causing Segmentation Fault !!






// 10. m04 & m05 & m06 

//sprintf(fname1, "%s/Coding/Node10/%s_m04.mp4",curdir,s1);
sprintf(fname1, "/mnt/node10/%s_m04.mp4",s1);
f_p1 = fopen(fname1, "rb");

//sprintf(fname2, "%s/Coding/Node10/%s_m05.mp4",curdir,s1);
sprintf(fname2, "/mnt/node10/%s_m05.mp4",s1);
f_p2 = fopen(fname2, "rb");

//sprintf(fname3, "%s/Coding/Node10/%s_m06.mp4",curdir,s1);
sprintf(fname2, "/mnt/node10/%s_m06.mp4",s1);
f_p3 = fopen(fname3, "rb");


gettimeofday(&t_read_start, &tz);

counter += jfread(block1, sizeof(char), blocksize, f_p1);
//fclose(f_p1); // causing Segmentation Fault !!

counter += jfread(block2, sizeof(char), blocksize, f_p2);
//fclose(f_p2); // causing Segmentation Fault !!

counter += jfread(block3, sizeof(char), blocksize, f_p3);


gettimeofday(&t_read_end, &tz);

tsec = 0.0;
tsec += t_read_end.tv_usec;
tsec -= t_read_start.tv_usec;
tsec /= 1000000.0;
tsec += t_read_end.tv_sec;
tsec -= t_read_start.tv_sec;
t_total_read = t_total_read + tsec;


// do measurement of bus transfer !!
gettimeofday(&t_bus_start, &tz);


ori_data[0] = block1 ;
ori_data[1] = block2 ;
ori_data[2] = block3 ;

for (i = 0; i < 2; i++) {
p_coding[i] = (char *)malloc(sizeof(char)*blocksize);
if (p_coding[i] == NULL) { perror("malloc"); exit(1); }
}


gettimeofday(&t_bus_end, &tz);
tsec = 0.0;
tsec += t_bus_end.tv_usec;
tsec -= t_bus_start.tv_usec;
tsec /= 1000000.0;
tsec += t_bus_end.tv_sec;
tsec -= t_bus_start.tv_sec;
//printf("Transfer Time (sec): %0.6f\n", tsec);
t_total_time = t_total_time + tsec;
//printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);



//printf(" counter = %d\n",counter);

// parity file
//sprintf(p_fname, "%s/Coding/NewNode%d/%s_parity_10.mp4",curdir,integer, s1);
//sprintf(p_fname, "/mnt/NewNode%d/%s_parity_10_%d.mp4",integer, s1 , z+1);
sprintf(p_fname, "/mnt/node%d/%s_parity_10_%d.mp4",integer+10, s1 , z+1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_cal_start, &tz);

galois_w08_region_multiply(ori_data[0], 1 , blocksize, p_coding[0], 0);
galois_w08_region_multiply(ori_data[1], 1 , blocksize, p_coding[0], 1);
//printf("galois_w08_region_multiply()!!\n");
galois_w08_region_multiply(ori_data[2], 1 , blocksize, p_coding[0], 1);



gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;


fwrite(p_coding[0], sizeof(char), blocksize, p_fp);
fclose(p_fp);


fclose(f_p1); // causing Segmentation Fault !!

fclose(f_p2); // causing Segmentation Fault !!









///////////////////////////////////////  NOW , 10 new parity made! /////////////////////////////////
////////////////////////////////////// and then, transfered them into NewNode1 , NewNode2 /////////////////


//////////////////////////////////////    NOW, xor between all parities !!!! ////////////////////////////////
/////////////////////////////    parity_01.mp4 ~ parity_08 at once !!    /////////////////////////////////

gettimeofday(&t_whcho_start, &tz);


int bit_counter;
//char *dptr , *sptr; 

counter=0;

par_data = (char **)malloc(sizeof(char*)*8);

//printf(" counter = %d\n",counter);

/////////////////// 1. read 8 parities into memory ////////////////////////
//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_01.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_01_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_01_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block1, sizeof(char), blocksize, f_p1);
par_data[0]= block1;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_02.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_02_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_02_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block2, sizeof(char), blocksize, f_p1);
par_data[1]= block2;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_03.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_03_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_03_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block3, sizeof(char), blocksize, f_p1);
par_data[2]= block3;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_04.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_04_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_04_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block4, sizeof(char), blocksize, f_p1);
par_data[3]= block4;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_05.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_05_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_05_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block5, sizeof(char), blocksize, f_p1);
par_data[4]= block5;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_06.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_06_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_06_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block6, sizeof(char), blocksize, f_p1);
par_data[5]= block6;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_07.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_07_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_07_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block7, sizeof(char), blocksize, f_p1);
par_data[6]= block7;

//sprintf(fname1, "%s/Coding/NewNode%d/%s_parity_08.mp4",curdir,integer,s1);
//sprintf(fname1, "/mnt/NewNode%d/%s_parity_08_%d.mp4",integer,s1,z+1);
sprintf(fname1, "/mnt/node%d/%s_parity_08_%d.mp4",integer+10,s1,z+1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block8, sizeof(char), blocksize, f_p1);
par_data[7]= block8;

// for debugging
//printf(" counter = %d!!!\n",counter);

gettimeofday(&t_whcho_end, &tz);
tsec = 0.0;
tsec += t_whcho_end.tv_usec;
tsec -= t_whcho_start.tv_usec;
tsec /= 1000000.0;
tsec += t_whcho_end.tv_sec;
tsec -= t_whcho_start.tv_sec;
t_whcho = t_whcho + tsec;




bit_counter=0;

dptr=p_coding[0];

gettimeofday(&t_cal_start, &tz);

for(i=0;i<8;i++){

sptr=par_data[i];

if(bit_counter == 0){
memcpy(dptr , sptr , blocksize);
bit_counter=1;
}
else {
galois_region_xor(sptr, dptr, blocksize);
}

}





gettimeofday(&t_cal_end, &tz);
tsec = 0.0;
tsec += t_cal_end.tv_usec;
tsec -= t_cal_start.tv_usec;
tsec /= 1000000.0;
tsec += t_cal_end.tv_sec;
tsec -= t_cal_start.tv_sec;
t_total_cal = t_total_cal + tsec;



//sprintf(p_fname, "%s/Coding/NewNode%d/%s_m%0*d%s", curdir, integer, s1, md, m+1+z, extension);
//sprintf(p_fname, "/mnt/NewNode%d/%s_m%0*d%s", integer, s1, md, m+1+z, extension);
sprintf(p_fname, "/mnt/node%d/%s_m%0*d%s", integer+10, s1, md, m+1+z, extension);

//sprintf(p_fname, "%s/Coding/%s_parity_final.mp4",curdir,s1);
p_fp = fopen(p_fname, "wb");

gettimeofday(&t_write_start, &tz);

fwrite(dptr, sizeof(char), blocksize, p_fp);
fclose(p_fp);

gettimeofday(&t_write_end, &tz);
tsec = 0.0;
tsec += t_write_end.tv_usec;
tsec -= t_write_start.tv_usec;
tsec /= 1000000.0;
tsec += t_write_end.tv_sec;
tsec -= t_write_start.tv_sec;
t_total_write = t_total_write + tsec;


}// for (z=0;z<4;z++) at LOC 900  




//fclose(p_fp);



printf("Total_Transfer Time (sec): %0.6f\n", t_total_time);
//printf("transger_throughput  (MB/sec): %0.6f\n", (((double) size)/1024.0/1024.0)/t_total_time);

printf("Total_Read Time (sec): %0.6f\n", t_total_read);
//printf("read_throughput (MB/sec): %0.6f\n", (((double) size)/1024.0/1024.0)/t_total_read);

printf("Total_Calculation Time (sec): %0.6f\n", t_total_cal);
//printf("cal_throughput (MB/sec): %0.6f\n", (((double) size)/1024.0/1024.0)/t_total_cal);

printf("Total_Write Time (sec): %0.6f\n", t_total_write);
//printf("write_throughput (MB/sec): %0.6f\n", (((double) size)/1024.0/1024.0)/t_total_write);

printf("whcho Time (sec): %0.6f\n", t_whcho);
//printf("whcho (MB/sec): %0.6f\n", (((double) size)/1024.0/1024.0)/t_whcho);







/* Calculate rate in MB/sec and print */



gettimeofday(&t6, &tz);
tsec = 0.0;
tsec += t6.tv_usec;
tsec -= t5.tv_usec;
tsec /= 1000000.0;
tsec += t6.tv_sec;
tsec -= t5.tv_sec;
printf("Calculate_Parity_Total Time (sec): %0.6f\n", tsec);





//fclose(f_p1);
//fclose(f_p2);

//fclose(src_fp);
//fclose(dst_fp);




//whcho added for testing !!!!


/*sprintf(fname1, "/mnt/NewNode1/%s_m05.mp4",s1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block1, sizeof(char), blocksize, f_p1);
par_data[0]= block1;

sprintf(fname1, "/mnt/NewNode1/%s_m06.mp4",s1);
f_p1 = fopen(fname1, "rb");
counter += jfread(block2, sizeof(char), blocksize, f_p1);
par_data[1]= block2;


dptr=par_data[0];
sptr=par_data[1];


galois_region_xor(sptr, dptr, blocksize);


sprintf(p_fname, "/mnt/NewNode1/%s_m_whcho", s1);
p_fp = fopen(p_fname, "wb");
fwrite(dptr, sizeof(char), blocksize, p_fp);
fclose(p_fp);
*/







free(block1);
free(block2);
free(block3);
free(block4);
free(block5);
free(block6);
free(block7);
free(block8);

free(ori_data);
free(p_coding);
free(par_data);




free(trans_block1);
free(trans_block2);
free(trans_block3);
free(trans_block4);
free(trans_block5);
free(trans_block6);
free(trans_block7);
free(trans_block8);


	/* Free allocated memory */
free(s1);
free(fname);
free(block);
free(curdir);



}








/* is_prime returns 1 if number if prime, 0 if not prime */
int is_prime(int w) {
	int prime55[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
	    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,
		    181,191,193,197,199,211,223,227,229,233,239,241,251,257};
	int i;
	for (i = 0; i < 55; i++) {
		if (w%prime55[i] == 0) {
			if (w == prime55[i]) return 1;
			else { return 0; }
		}
	}
}

/* Handles ctrl-\ event */
void ctrl_bs_handler(int dummy) {
	time_t mytime;
	mytime = time(0);
	fprintf(stderr, "\n%s\n", ctime(&mytime));
	fprintf(stderr, "You just typed ctrl-\\ in encoder.c.\n");
	fprintf(stderr, "Total number of read ins = %d\n", readins);
	fprintf(stderr, "Current read in: %d\n", n);
	fprintf(stderr, "Method: %s\n\n", Methods[method]);	
	signal(SIGQUIT, ctrl_bs_handler);
}
