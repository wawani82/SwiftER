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

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gf_complete.h>
#include "galois.h"
#include "jerasure.h"
#include "reed_sol.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

int *reed_sol_r6_coding_matrix(int k, int w)
{
  int *matrix;
  int i, tmp;

  if (w != 8 && w != 16 && w != 32) return NULL;

  matrix = talloc(int, 2*k);
  if (matrix == NULL) return NULL;

  for (i = 0; i < k; i++) matrix[i] = 1;
  matrix[k] = 1;
  tmp = 1;
  for (i = 1; i < k; i++) {
    tmp = galois_single_multiply(tmp, 2, w);
    matrix[k+i] = tmp;
  }
  return matrix;
}


// for encoding 
int *reed_sol_vandermonde_coding_matrix(int k, int m, int w)
{
  int tmp;
  int i, j, index;
  int *vdm, *dist;

//whcho added
  /*vdm = reed_sol_big_vandermonde_distribution_matrix(k+m, k, w);
  if (vdm == NULL) return NULL;
  dist = talloc(int, m*k);
  if (dist == NULL) {
    free(vdm);
    return NULL;
  }*/

//whcho added
  //whcho added
  dist = talloc(int, m*k);
  vdm = talloc(int, (k+m)*k );
  
  //for (j=0; j<k ;j++){
    //vdm[j]=1;
  //}

/* Distibution matrix */
vdm[0]=1;
vdm[24*1]=1;
vdm[24*2]=1;
vdm[24*3]=3;
vdm[24*4]=4;
vdm[24*5]=8;



vdm[1]=1;
vdm[24*1+1]=2;
vdm[24*2+1]=2;
vdm[24*3+1]=2;
vdm[24*4+1]=7;
vdm[24*5+1]=14;



vdm[2]=1;
vdm[24*1+2]=1;
vdm[24*2+2]=3;
vdm[24*3+2]=2;
vdm[24*4+2]=7;
vdm[24*5+2]=8;



vdm[3]=1;
vdm[24*1+3]=2;
vdm[24*2+3]=1;
vdm[24*3+3]=3;
vdm[24*4+3]=7;
vdm[24*5+3]=11;


vdm[4]=1;
vdm[24*1+4]=2;
vdm[24*2+4]=2;
vdm[24*3+4]=3;
vdm[24*4+4]=6;
vdm[24*5+4]=13;


vdm[5]=2;
vdm[24*1+5]=2;
vdm[24*2+5]=4;
vdm[24*3+5]=1;
vdm[24*4+5]=9;
vdm[24*5+5]=21;



vdm[6]=2;
vdm[24*1+6]=3;
vdm[24*2+6]=4;
vdm[24*3+6]=1;
vdm[24*4+6]=12;
vdm[24*5+6]=22;

vdm[7]=1;
vdm[24*1+7]=1;
vdm[24*2+7]=1;
vdm[24*3+7]=2;
vdm[24*4+7]=9;
vdm[24*5+7]=15;

vdm[8]=2;
vdm[24*1+8]=3;
vdm[24*2+8]=2;
vdm[24*3+8]=3;
vdm[24*4+8]=17;
vdm[24*5+8]=19;


vdm[9]=3;
vdm[24*1+9]=4;
vdm[24*2+9]=4;
vdm[24*3+9]=5;
vdm[24*4+9]=14;
vdm[24*5+9]=23;

vdm[10]=3;
vdm[24*1+10]=2;
vdm[24*2+10]=1;
vdm[24*3+10]=5;
vdm[24*4+10]=16;
vdm[24*5+10]=20;

vdm[11]=1;
vdm[24*1+11]=4;
vdm[24*2+11]=2;
vdm[24*3+11]=3;
vdm[24*4+11]=15;
vdm[24*5+11]=22;


vdm[12]=3;
vdm[24*1+12]=4;
vdm[24*2+12]=4;
vdm[24*3+12]=3;
vdm[24*4+12]=12;
vdm[24*5+12]=26;

vdm[13]=2;
vdm[24*1+13]=3;
vdm[24*2+13]=1;
vdm[24*3+13]=5;
vdm[24*4+13]=12;
vdm[24*5+13]=32;

vdm[14]=1;
vdm[24*1+14]=3;
vdm[24*2+14]=1;
vdm[24*3+14]=2;
vdm[24*4+14]=23;
vdm[24*5+14]=18;


vdm[15]=2;
vdm[24*1+15]=4;
vdm[24*2+15]=4;
vdm[24*3+15]=3;
vdm[24*4+15]=14;
vdm[24*5+15]=21;

vdm[16]=3;
vdm[24*1+16]=3;
vdm[24*2+16]=1;
vdm[24*3+16]=2;
vdm[24*4+16]=17;
vdm[24*5+16]=22;

vdm[17]=3;
vdm[24*1+17]=4;
vdm[24*2+17]=4;
vdm[24*3+17]=5;
vdm[24*4+17]=7;
vdm[24*5+17]=11;


vdm[18]=1;
vdm[24*1+18]=1;
vdm[24*2+18]=1;
vdm[24*3+18]=1;
vdm[24*4+18]=9;
vdm[24*5+18]=7;

vdm[19]=2;
vdm[24*1+19]=1;
vdm[24*2+19]=5;
vdm[24*3+19]=3;
vdm[24*4+19]=12;
vdm[24*5+19]=11;

vdm[20]=4;
vdm[24*1+20]=4;
vdm[24*2+20]=5;
vdm[24*3+20]=4;
vdm[24*4+20]=7;
vdm[24*5+20]=9;



vdm[21]=3;
vdm[24*1+21]=4;
vdm[24*2+21]=4;
vdm[24*3+21]=4;
vdm[24*4+21]=17;
vdm[24*5+21]=9;

vdm[22]=3;
vdm[24*1+22]=3;
vdm[24*2+22]=1;
vdm[24*3+22]=3;
vdm[24*4+22]=19;
vdm[24*5+22]=7;

vdm[23]=3;
vdm[24*1+23]=5;
vdm[24*2+23]=5;
vdm[24*3+23]=2;
vdm[24*4+23]=23;
vdm[24*5+23]=7;




  //i = k*k;
  for (j = 0; j < m*k; j++) {
    dist[j] = vdm[j];
    //i++;
  }
  
  /*i = k*k;
  for (j = 0; j < m*k; j++) {
    dist[j] = vdm[i];
    i++;
  }*/

  free(vdm);
  return dist;
}



//for decoding 
int *reed_sol_vandermonde_decoding_matrix(int k, int m, int w)
{
  int tmp;
  int i, j, index;
  int *vdm, *dist;

//whcho added
  /*vdm = reed_sol_big_vandermonde_distribution_matrix(k+m, k, w);
  if (vdm == NULL) return NULL;
  dist = talloc(int, m*k);
  if (dist == NULL) {
    free(vdm);
    return NULL;
  }*/

//whcho added
  //whcho added
  dist = talloc(int, m*k);
  vdm = talloc(int, (k+m)*k );
  
  //for (j=0; j<k ;j++){
    //vdm[j]=1;
  //}

/* Distibution matrix */
vdm[0]=1;
vdm[24*1]=1;
vdm[24*2]=1;
vdm[24*3]=3;
vdm[24*4]=4;
vdm[24*5]=8;
vdm[24*6]=30;
vdm[24*7]=3;
vdm[24*8]=27;
vdm[24*9]=24;
vdm[24*10]=26;
vdm[24*11]=15;



vdm[1]=1;
vdm[24*1+1]=2;
vdm[24*2+1]=2;
vdm[24*3+1]=2;
vdm[24*4+1]=7;
vdm[24*5+1]=14;
vdm[24*6+1]=48;
vdm[24*7+1]=5;
vdm[24*8+1]=43;
vdm[24*9+1]=37;
vdm[24*10+1]=40;
vdm[24*11+1]=23;



vdm[2]=1;
vdm[24*1+2]=1;
vdm[24*2+2]=3;
vdm[24*3+2]=2;
vdm[24*4+2]=7;
vdm[24*5+2]=8;
vdm[24*6+2]=30;
vdm[24*7+2]=5;
vdm[24*8+2]=28;
vdm[24*9+2]=34;
vdm[24*10+2]=38;
vdm[24*11+2]=17;



vdm[3]=1;
vdm[24*1+3]=2;
vdm[24*2+3]=1;
vdm[24*3+3]=3;
vdm[24*4+3]=7;
vdm[24*5+3]=11;
vdm[24*6+3]=28;
vdm[24*7+3]=13;
vdm[24*8+3]=31;
vdm[24*9+3]=31;
vdm[24*10+3]=40;
vdm[24*11+3]=30;


vdm[4]=1;
vdm[24*1+4]=2;
vdm[24*2+4]=2;
vdm[24*3+4]=3;
vdm[24*4+4]=6;
vdm[24*5+4]=13;
vdm[24*6+4]=29;
vdm[24*7+4]=11;
vdm[24*8+4]=31;
vdm[24*9+4]=35;
vdm[24*10+4]=47;
vdm[24*11+4]=34;


vdm[5]=2;
vdm[24*1+5]=2;
vdm[24*2+5]=4;
vdm[24*3+5]=1;
vdm[24*4+5]=9;
vdm[24*5+5]=21;
vdm[24*6+5]=42;
vdm[24*7+5]=17;
vdm[24*8+5]=45;
vdm[24*9+5]=43;
vdm[24*10+5]=49;
vdm[24*11+5]=37;



vdm[6]=2;
vdm[24*1+6]=3;
vdm[24*2+6]=4;
vdm[24*3+6]=1;
vdm[24*4+6]=12;
vdm[24*5+6]=22;
vdm[24*6+6]=52;
vdm[24*7+6]=13;
vdm[24*8+6]=58;
vdm[24*9+6]=44;
vdm[24*10+6]=53;
vdm[24*11+6]=41;


vdm[7]=1;
vdm[24*1+7]=1;
vdm[24*2+7]=1;
vdm[24*3+7]=2;
vdm[24*4+7]=9;
vdm[24*5+7]=15;
vdm[24*6+7]=45;
vdm[24*7+7]=11;
vdm[24*8+7]=57;
vdm[24*9+7]=29;
vdm[24*10+7]=41;
vdm[24*11+7]=34;

vdm[8]=2;
vdm[24*1+8]=3;
vdm[24*2+8]=2;
vdm[24*3+8]=3;
vdm[24*4+8]=17;
vdm[24*5+8]=19;
vdm[24*6+8]=58;
vdm[24*7+8]=13;
vdm[24*8+8]=67;
vdm[24*9+8]=46;
vdm[24*10+8]=61;
vdm[24*11+8]=49;


vdm[9]=3;
vdm[24*1+9]=4;
vdm[24*2+9]=4;
vdm[24*3+9]=5;
vdm[24*4+9]=14;
vdm[24*5+9]=23;
vdm[24*6+9]=74;
vdm[24*7+9]=26;
vdm[24*8+9]=65;
vdm[24*9+9]=65;
vdm[24*10+9]=61;
vdm[24*11+9]=62;


vdm[10]=3;
vdm[24*1+10]=2;
vdm[24*2+10]=1;
vdm[24*3+10]=5;
vdm[24*4+10]=16;
vdm[24*5+10]=20;
vdm[24*6+10]=61;
vdm[24*7+10]=16;
vdm[24*8+10]=55;
vdm[24*9+10]=53;
vdm[24*10+10]=51;
vdm[24*11+10]=51;


vdm[11]=1;
vdm[24*1+11]=4;
vdm[24*2+11]=2;
vdm[24*3+11]=3;
vdm[24*4+11]=15;
vdm[24*5+11]=22;
vdm[24*6+11]=61;
vdm[24*7+11]=17;
vdm[24*8+11]=55;
vdm[24*9+11]=62;
vdm[24*10+11]=57;
vdm[24*11+11]=65;


vdm[12]=3;
vdm[24*1+12]=4;
vdm[24*2+12]=4;
vdm[24*3+12]=3;
vdm[24*4+12]=12;
vdm[24*5+12]=26;
vdm[24*6+12]=52;
vdm[24*7+12]=17;
vdm[24*8+12]=62;
vdm[24*9+12]=58;
vdm[24*10+12]=60;
vdm[24*11+12]=43;


vdm[13]=2;
vdm[24*1+13]=3;
vdm[24*2+13]=1;
vdm[24*3+13]=5;
vdm[24*4+13]=12;
vdm[24*5+13]=32;
vdm[24*6+13]=55;
vdm[24*7+13]=18;
vdm[24*8+13]=75;
vdm[24*9+13]=70;
vdm[24*10+13]=75;
vdm[24*11+13]=54;


vdm[14]=1;
vdm[24*1+14]=3;
vdm[24*2+14]=1;
vdm[24*3+14]=2;
vdm[24*4+14]=23;
vdm[24*5+14]=18;
vdm[24*6+14]=48;
vdm[24*7+14]=14;
vdm[24*8+14]=63;
vdm[24*9+14]=57;
vdm[24*10+14]=60;
vdm[24*11+14]=64;


vdm[15]=2;
vdm[24*1+15]=4;
vdm[24*2+15]=4;
vdm[24*3+15]=3;
vdm[24*4+15]=14;
vdm[24*5+15]=21;
vdm[24*6+15]=57;
vdm[24*7+15]=13;
vdm[24*8+15]=63;
vdm[24*9+15]=58;
vdm[24*10+15]=68;
vdm[24*11+15]=63;


vdm[16]=3;
vdm[24*1+16]=3;
vdm[24*2+16]=1;
vdm[24*3+16]=2;
vdm[24*4+16]=17;
vdm[24*5+16]=22;
vdm[24*6+16]=54;
vdm[24*7+16]=9;
vdm[24*8+16]=58;
vdm[24*9+16]=52;
vdm[24*10+16]=56;
vdm[24*11+16]=51;


vdm[17]=3;
vdm[24*1+17]=4;
vdm[24*2+17]=4;
vdm[24*3+17]=5;
vdm[24*4+17]=7;
vdm[24*5+17]=11;
vdm[24*6+17]=46;
vdm[24*7+17]=15;
vdm[24*8+17]=54;
vdm[24*9+17]=42;
vdm[24*10+17]=50;
vdm[24*11+17]=43;


vdm[18]=1;
vdm[24*1+18]=1;
vdm[24*2+18]=1;
vdm[24*3+18]=1;
vdm[24*4+18]=9;
vdm[24*5+18]=7;
vdm[24*6+18]=28;
vdm[24*7+18]=15;
vdm[24*8+18]=36;
vdm[24*9+18]=36;
vdm[24*10+18]=28;
vdm[24*11+18]=29;


vdm[19]=2;
vdm[24*1+19]=1;
vdm[24*2+19]=5;
vdm[24*3+19]=3;
vdm[24*4+19]=12;
vdm[24*5+19]=11;
vdm[24*6+19]=38;
vdm[24*7+19]=14;
vdm[24*8+19]=42;
vdm[24*9+19]=54;
vdm[24*10+19]=44;
vdm[24*11+19]=41;


vdm[20]=4;
vdm[24*1+20]=4;
vdm[24*2+20]=5;
vdm[24*3+20]=4;
vdm[24*4+20]=7;
vdm[24*5+20]=9;
vdm[24*6+20]=37;
vdm[24*7+20]=19;
vdm[24*8+20]=41;
vdm[24*9+20]=45;
vdm[24*10+20]=39;
vdm[24*11+20]=29;



vdm[21]=3;
vdm[24*1+21]=4;
vdm[24*2+21]=4;
vdm[24*3+21]=4;
vdm[24*4+21]=17;
vdm[24*5+21]=9;
vdm[24*6+21]=51;
vdm[24*7+21]=19;
vdm[24*8+21]=55;
vdm[24*9+21]=56;
vdm[24*10+21]=51;
vdm[24*11+21]=50;


vdm[22]=3;
vdm[24*1+22]=3;
vdm[24*2+22]=1;
vdm[24*3+22]=3;
vdm[24*4+22]=19;
vdm[24*5+22]=7;
vdm[24*6+22]=51;
vdm[24*7+22]=19;
vdm[24*8+22]=57;
vdm[24*9+22]=48;
vdm[24*10+22]=44;
vdm[24*11+22]=45;


vdm[23]=3;
vdm[24*1+23]=5;
vdm[24*2+23]=5;
vdm[24*3+23]=2;
vdm[24*4+23]=23;
vdm[24*5+23]=7;
vdm[24*6+23]=60;
vdm[24*7+23]=25;
vdm[24*8+23]=66;
vdm[24*9+23]=66;
vdm[24*10+23]=59;
vdm[24*11+23]=60;




  //i = k*k;
  for (j = 0; j < m*k; j++) {
    dist[j] = vdm[j];
    //i++;
  }
  
  /*i = k*k;
  for (j = 0; j < m*k; j++) {
    dist[j] = vdm[i];
    i++;
  }*/

  free(vdm);
  return dist;


}


static int prim08 = -1;
static gf_t GF08;

void reed_sol_galois_w08_region_multby_2(char *region, int nbytes)
{
  if (prim08 == -1) {
    prim08 = galois_single_multiply((1 << 7), 2, 8);
    if (!gf_init_hard(&GF08, 8, GF_MULT_BYTWO_b, GF_REGION_DEFAULT, GF_DIVIDE_DEFAULT,
                      prim08, 0, 0, NULL, NULL)) {
      fprintf(stderr, "Error: Can't initialize the GF for reed_sol_galois_w08_region_multby_2\n");
      exit(1);
    }
  }
  GF08.multiply_region.w32(&GF08, region, region, 2, nbytes, 0);
}

static int prim16 = -1;
static gf_t GF16;

void reed_sol_galois_w16_region_multby_2(char *region, int nbytes)
{
  if (prim16 == -1) {
    prim16 = galois_single_multiply((1 << 15), 2, 16);
    if (!gf_init_hard(&GF16, 16, GF_MULT_BYTWO_b, GF_REGION_DEFAULT, GF_DIVIDE_DEFAULT,
                      prim16, 0, 0, NULL, NULL)) {
      fprintf(stderr, "Error: Can't initialize the GF for reed_sol_galois_w16_region_multby_2\n");
      exit(1);
    }
  }
  GF16.multiply_region.w32(&GF16, region, region, 2, nbytes, 0);
}

static int prim32 = -1;
static gf_t GF32;

void reed_sol_galois_w32_region_multby_2(char *region, int nbytes)
{
  if (prim32 == -1) {
    prim32 = galois_single_multiply((1 << 31), 2, 32);
    if (!gf_init_hard(&GF32, 32, GF_MULT_BYTWO_b, GF_REGION_DEFAULT, GF_DIVIDE_DEFAULT,
                      prim32, 0, 0, NULL, NULL)) {
      fprintf(stderr, "Error: Can't initialize the GF for reed_sol_galois_w32_region_multby_2\n");
      exit(1);
    }
  }
  GF32.multiply_region.w32(&GF32, region, region, 2, nbytes, 0);
}

int reed_sol_r6_encode(int k, int w, char **data_ptrs, char **coding_ptrs, int size)
{
  int i;

  /* First, put the XOR into coding region 0 */

  memcpy(coding_ptrs[0], data_ptrs[0], size);

  for (i = 1; i < k; i++) galois_region_xor(data_ptrs[i], coding_ptrs[0], size);

  /* Next, put the sum of (2^j)*Dj into coding region 1 */

  memcpy(coding_ptrs[1], data_ptrs[k-1], size);

  for (i = k-2; i >= 0; i--) {
    switch (w) {
      case 8:  reed_sol_galois_w08_region_multby_2(coding_ptrs[1], size); break;
      case 16: reed_sol_galois_w16_region_multby_2(coding_ptrs[1], size); break;
      case 32: reed_sol_galois_w32_region_multby_2(coding_ptrs[1], size); break;
      default: return 0;
    }

    galois_region_xor(data_ptrs[i], coding_ptrs[1], size);
  }
  return 1;
}

int *reed_sol_extended_vandermonde_matrix(int rows, int cols, int w)
{
  int *vdm;
  int i, j, k;

  if (w < 30 && (1 << w) < rows) return NULL;
  if (w < 30 && (1 << w) < cols) return NULL;

  vdm = talloc(int, rows*cols);
  if (vdm == NULL) { return NULL; }
  
  vdm[0] = 1;
  for (j = 1; j < cols; j++) vdm[j] = 0;
  if (rows == 1) return vdm;

  i=(rows-1)*cols;
  for (j = 0; j < cols-1; j++) vdm[i+j] = 0;
  vdm[i+j] = 1;
  if (rows == 2) return vdm;

  for (i = 1; i < rows-1; i++) {
    k = 1;
    for (j = 0; j < cols; j++) {
      vdm[i*cols+j] = k;
      k = galois_single_multiply(k, i, w);
    }
  }
  return vdm;
}

int *reed_sol_big_vandermonde_distribution_matrix(int rows, int cols, int w)
{
  int *dist;
  int i, j, k;
  int sindex, srindex, siindex, tmp;

  if (cols >= rows) return NULL;
  
  dist = reed_sol_extended_vandermonde_matrix(rows, cols, w);
  if (dist == NULL) return NULL;

  sindex = 0;
  for (i = 1; i < cols; i++) {
    sindex += cols;

    /* Find an appropriate row -- where i,i != 0 */
    srindex = sindex+i;
    for (j = i; j < rows && dist[srindex] == 0; j++) srindex += cols;
    if (j >= rows) {   /* This should never happen if rows/w are correct */
      fprintf(stderr, "reed_sol_big_vandermonde_distribution_matrix(%d,%d,%d) - couldn't make matrix\n", 
             rows, cols, w);
      exit(1);
    }
 
    /* If necessary, swap rows */
    if (j != i) {
      srindex -= i;
      for (k = 0; k < cols; k++) {
        tmp = dist[srindex+k];
        dist[srindex+k] = dist[sindex+k];
        dist[sindex+k] = tmp;
      }
    }
  
    /* If Element i,i is not equal to 1, multiply the column by 1/i */

    if (dist[sindex+i] != 1) {
      tmp = galois_single_divide(1, dist[sindex+i], w);
      srindex = i;
      for (j = 0; j < rows; j++) {
        dist[srindex] = galois_single_multiply(tmp, dist[srindex], w);
        srindex += cols;
      }
    }
 
    /* Now, for each element in row i that is not in column 1, you need
       to make it zero.  Suppose that this is column j, and the element
       at i,j = e.  Then you want to replace all of column j with 
       (col-j + col-i*e).   Note, that in row i, col-i = 1 and col-j = e.
       So (e + 1e) = 0, which is indeed what we want. */

    for (j = 0; j < cols; j++) {
      tmp = dist[sindex+j];
      if (j != i && tmp != 0) {
        srindex = j;
        siindex = i;
        for (k = 0; k < rows; k++) {
          dist[srindex] = dist[srindex] ^ galois_single_multiply(tmp, dist[siindex], w);
          srindex += cols;
          siindex += cols;
        }
      }
    }
  }
  /* We desire to have row k be all ones.  To do that, multiply
     the entire column j by 1/dist[k,j].  Then row j by 1/dist[j,j]. */

  sindex = cols*cols;
  for (j = 0; j < cols; j++) {
    tmp = dist[sindex];
    if (tmp != 1) { 
      tmp = galois_single_divide(1, tmp, w);
      srindex = sindex;
      for (i = cols; i < rows; i++) {
        dist[srindex] = galois_single_multiply(tmp, dist[srindex], w);
        srindex += cols;
      }
    }
    sindex++;
  }

  /* Finally, we'd like the first column of each row to be all ones.  To
     do that, we multiply the row by the inverse of the first element. */

  sindex = cols*(cols+1);
  for (i = cols+1; i < rows; i++) {
    tmp = dist[sindex];
    if (tmp != 1) { 
      tmp = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) dist[sindex+j] = galois_single_multiply(dist[sindex+j], tmp, w);
    }
    sindex += cols;
  }

  return dist;
}

