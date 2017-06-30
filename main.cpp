#include "smatrix.h"
#include <ctime>
using namespace std;


void printet(et *e,int n){
	int s;
	for(int k=0;k<n;k++){
		printf("par[%d] is %d,chl are ",k,e->par[k]);
		s=(int)e->chl->at(k)->size();
		for(int j=0;j<s;j++){
			printf("%d ",e->chl->at(k)->at(j));

		}
		printf("\n");
	}
}
void print_cs(map<int,vector<int>* > *e){
	int s,n=e->size();
	for(int k=0;k<n;k++){
		printf("col_struct[%d] is  ",k);
		s=(int)e->at(k)->size();
		for(int j=0;j<s;j++){
			printf("%d ",e->at(k)->at(j));

		}
		printf("\n");
	}	
}


int main(){
	sm *A;
	FILE *f;
	f=fopen("input.txt","r");
	A=load(f);
	fclose(f);
	sm *B=TripletToCC(A);
	int n=B->n;
	printf("given matrix A:\n");

	

	clock_t t = clock();

	et *e=buildTree(B);

	map<int,vector<int>* > * col_strut=col_structure(e,B);

	sm *L=sm_chol(B,e,col_strut);

	printf("\nTime taken: %.4fs\n", (float)(clock() - t)/CLOCKS_PER_SEC);


	printf("chol factor L:\n");
	for(int j=0;j<n;j++){
		for(int p=L->p[j];p<L->p[j+1];p++){
			printf("%d %d %g\n",L->i[p]+1,j+1,L->x[p]);
		}
	}

	return 0;

}
