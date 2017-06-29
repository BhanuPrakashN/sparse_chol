#ifndef _SMATRIX_H_INCLUDED__
#define _SMATRIX_H_INCLUDED__
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <algorithm>  
#include <map>
#include <stdio.h>
#include <stddef.h>
using namespace std;
typedef struct smatrix{
	int nzmax;
	int m;
	int n;
	int *p;
	int *i;
	double *x;
	int nz;
}sm;


typedef struct etree{
	int *par;
	std::map<int,vector<int>* >  *chl;
}et;
et* buildTree(sm *A);




typedef struct Fmatrix{
	int n;
	vector<int> *indices;
	double **x;

}Fm;




sm *smfree (sm *A);
int sm_insert(sm *A,int i,int j,double x);
void smrealloc (sm *A, int nzmax);
sm *smalloc(int m, int n, int nzmax, int values, int triplet);
sm *load (FILE *f);
int print(sm *A);
int cumsum(int *a,int *b,int n);
sm *TripletToCC(sm *T);
map<int,vector<int>* > * col_structure(et *e ,sm *A);
sm *sm_chol(sm *A,et *e,map<int,vector<int>* > * cs);
Fm *fmalloc(int n,vector<int> *indices);
double sm_access(sm *A,int i,int j);
void fconstruct(Fm *F,sm *A,int j);
void printFm(Fm *F);
void Fupdate(Fm *F,Fm *U);


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CSC(A) (A && (A->nz == -1))
#define TRIPLET(A) (A && (A->nz >= 0))
#endif