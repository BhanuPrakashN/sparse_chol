#include "smatrix.h"


int sm_insert(sm *A,int i,int j,double x){ //insert entry(i,j,x) into triplet form
	if(!TRIPLET(A)||i<0||j<0) return 0;
	if(A->nz >= A->nzmax)smrealloc(A,2*(A->nzmax));
	A->i[A->nz]=i;
	A->p[A->nz]=j;
	A->x[A->nz]=x;
	A->nz++;
	A->m=MAX(A->m,i+1);
	A->n=MAX(A->n,j+1);
	return 1;
}

void smrealloc (sm *A, int nzmax)
{
    A->nzmax = (nzmax <= 0) ? (A->p [A->n]) : nzmax ;
    A->i = (int*)realloc (A->i, MAX(nzmax,1)*sizeof(int)) ;
    if (A->nz >= 0) A->p = (int*)realloc (A->p,MAX(nzmax,1)*sizeof(int)) ;
    if (A->x) A->x = (double*)realloc (A->x, MAX(nzmax,1)*sizeof(double)) ;
}



sm *smalloc(int m, int n, int nzmax, int values, int triplet){
	sm *A = (sm*)calloc(1,sizeof(sm));
	A->m = m ;
	A->n = n ;
	A->nzmax = nzmax;
	A->nz = triplet? 0:-1;
	A->p=(int*)malloc((triplet? nzmax : n+1)*sizeof(int));
	A->i=(int*)malloc(nzmax*sizeof(int));
	A->x=values ? (double*)malloc(nzmax*sizeof(double)) : NULL ;
	return ((!A->p || !A->i || (values && !A->x)) ? smfree (A) : A) ;
}




sm *smfree (sm *A){
    if (!A) return (NULL) ;	
    free (A->p) ;
    free (A->i) ;
    free (A->x) ;
    if(A)free(A);
    return NULL;	
}


sm *load (FILE *f){
    int i, j ;   
    double x ;
    sm *A ;
    if (!f) return (NULL) ;                             
    A = smalloc (0, 0, 1, 1, 1) ;                    
    while (fscanf (f, "%d %d %lg\n", &i, &j, &x) == 3)
    {
        if (!sm_insert (A,  i-1,  j-1, x)) return (smfree (A)) ;
    }
    return (A) ;
}


int print(sm *A){
    printf("triplet of %d by %d\n",A->m,A->n);
    printf("  i      j      x \n");
    for(int p=0;p<A->nz;p++){
        printf("%6d %6d %6g\n",A->i[p],A->p[p],A->x[p]);

    }
    return 0;
}


sm *TripletToCC(sm *T){
    int m,n,nz,*Ti,*Tj,*Cp,*Ci,*temp,p;
    double *Cx,*Tx;
    sm *C;
    if(!TRIPLET(T))return NULL;
    m=T->m;
    n=T->n;
    Ti=T->i;
    Tj=T->p;
    Tx=T->x;
    nz=T->nz;
    C=smalloc(m,n,nz,Tx!=NULL,0);
    temp=(int*)calloc(n,sizeof(int));
    Cp=C->p;
    Ci=C->i;
    Cx=C->x;
    //printf("n=%d\n",C->n);
    for(int k=0;k<nz;k++)temp[Tj[k]]++;
    Cp[n]=cumsum(Cp,temp,n);
    for(int k=0;k<nz;k++){
        p=temp[Tj[k]]++;
        Ci[p]=Ti[k];
        Cx[p]=Tx[k];
    }
    return C;

}

int cumsum(int *a,int *b,int n){
    int nz=0;
    for(int i=0;i<n;i++){
        a[i]=nz;
        nz+=b[i];
        b[i]=a[i];
    }
    return nz;
}



et * buildTree(sm *A){
    et *tree=(et*)malloc(sizeof(et));
    int n=A->n,i=0,next=0;
    int *par,*Ap,*Ai,*ancestor;
    map<int,vector<int>* >  *chl;
    vector<int> *t;
    tree->par=(int*)malloc(n*sizeof(int));
    tree->chl=new map<int,vector<int>* >;
    ancestor=(int*)malloc(n*sizeof(int));
    par=tree->par;
    chl=tree->chl;
    Ap=A->p;
    Ai=A->i;
    for(int k=0;k<n;k++){
        par[k]=-1;
        ancestor[k]=-1;
        t= new vector<int>;
        chl->insert(pair<int,vector<int>* >(k,t));
        for(int p=Ap[k];p<Ap[k+1];p++){
            i=Ai[p];
            while(i!=-1 && i<k){
               // printf("%d\n",i);
                next=ancestor[i];
                ancestor[i]=k;
                if(next==-1){

                    par[i]=k;
                    (chl->at(k))->push_back(i);
                }
                i=next;
            }
        }
    }
    return tree;
}


map<int,vector<int>* > * col_structure(et *e ,sm *A){
    int n=A->n,j;
    map<int,vector<int>* > * col_strut=new map<int,vector<int>* >;
    vector<int> *t;
    int *Ap,*Ai;
    Ap=A->p;
    Ai=A->i;
    int *mark=(int*)malloc(n*sizeof(int));
    for(j=0;j<n;j++){
        t=new vector<int>;
        t->push_back(j);
        col_strut->insert(pair<int,vector<int>* >(j,t));

    }
    for(int i=0;i<n;i++){

        mark[i]=i;
        for(int k=Ap[i];k<Ap[i+1];k++){
            if(Ai[k]<i){
                j=Ai[k];
                while(mark[j]!=i && j!=-1){
                    t=col_strut->at(j);
                    t->push_back(i);
                    mark[j]=i;
                    j=e->par[j];
                }

            }
        }

    }
    return col_strut;
}
double sm_access(sm *A,int i,int j){
    for(int p=A->p[j];p < A->p[j+1];p++){
        if(A->i[p]==i)return A->x[p];
    }
    return 0;
}

Fm *fmalloc(int n,vector<int> *indices){
    Fm *F=(Fm*)malloc(sizeof(Fm));
    F->n=n;
    if(indices!=NULL)sort(indices->begin(),indices->end());
    F->indices=indices;
    F->x=(double **)calloc(n,sizeof(double *));
    for(int i=0;i<n;i++)
        F->x[i]=(double*)calloc(n,sizeof(double));
    return F;
}
void fconstruct(Fm *F,sm *A,int j){
    int n=F->n,i;
    vector<int> *ind=F->indices;
    double **Fx=F->x;
    Fx[0][0]=sm_access(A,j,j);
    for(int p=1;p<n;p++){
        i=ind->at(p);
        Fx[0][p]=sm_access(A,i,j);
        Fx[p][0]=sm_access(A,i,j);
    }


}
void printFm(Fm *F){
    for(int  i=0;i<F->n;i++){
        for(int j=0;j<F->n;j++){
            printf("%6.4g ",F->x[i][j]);
        }
        printf("\n");
    }
}
void Fupdate(Fm *F,Fm *U){
    vector<int> *ind=U->indices;
    int i,j;
    for(int k=0;k<U->n;k++){
        i=ind->at(k);
        for(int p=0;p<U->n;p++){
            j=ind->at(p);
            F->x[i][j]+=U->x[k][p];
        }
    }

}
Fm *UpdateF(Fm *F,double *l_b,vector<int> *Pind){
    Fm *U=(Fm*)malloc(sizeof(Fm));
    U->n=(F->n)-1;
    U->x=(double **)calloc(U->n,sizeof(double *));
    for(int i=0;i<U->n;i++)
        U->x[i]=(double*)calloc(U->n,sizeof(double));
    vector<int> *Gind=F->indices;
    U->indices =new vector<int>;
    int temp=0;
    if(Pind!=NULL)sort(Pind->begin(),Pind->end());
    //printf("Pind:");
    //for(int i=0;i<Pind->size();i++)
      //  printf("%d ",Pind->at(i) );
    //printf("\nGind:");
    //for(int i=0;i<Gind->size();i++)
      //  printf("%d ",Gind->at(i) );  
    

    for(int i=1;i<F->n;i++){
        temp=0;
        while(Pind->at(temp)!=Gind->at(i))temp++;
        U->indices->push_back(temp);
    }
    //printf("\nLind:");
    //for(int i=0;i<U->indices->size();i++)
    //    printf("%d ",U->indices->at(i) );
    for(int i=0;i<U->n;i++){
        for(int j=0;j<U->n;j++){
            U->x[i][j]=F->x[i+1][j+1]-l_b[i]*l_b[j];
        }
    }

    return U;



}
sm *sm_chol(sm *A,et *e,map<int,vector<int>* > * cs){
    sm *L=smalloc (0, 0, 1, 1, 1);
    double l_jj,*l_b;
    int n,c;
    n=A->n;
    vector<int> *temp,*chl;
    Fm *F,*Um;
    map<int,Fm*> *U=new map<int,Fm*>;
    for(int j=0;j<n;j++){
        temp=cs->at(j);
        F=fmalloc(temp->size(),temp);
        fconstruct(F,A,j);
        //printf("\n Frontal[%d]\n",j);
        //for(int k=0;k<F->n;k++)printf("%d ",temp->at(k));
        
        chl=e->chl->at(j);
        for(int k=0;k<chl->size();k++){
            c=chl->at(k);
            Fupdate(F,U->at(c));
            U->erase(c);
        }
        //printFm(F);
        l_b=(double*)malloc(((F->n)-1)*sizeof(double));
        l_jj=sqrt(F->x[0][0]);
        sm_insert (L,  j,  j, l_jj);

        //printf("%g\n",*l_jj);
        for(int k=1;k<F->n;k++){
            l_b[k-1]=F->x[0][k]/(l_jj);
            sm_insert (L,  F->indices->at(k),  j, l_b[k-1]);
            //printf("%g\n",l_b[k-1]);
        }

        if(e->par[j] != -1){
            Um=UpdateF(F,l_b,cs->at(e->par[j]));
            // printf("\n Update[%d]\n",j);
            //printFm(Um);
           
            U->insert(pair<int,Fm*>(j,Um));
        }

    }
    L=TripletToCC(L);
    return L;



}