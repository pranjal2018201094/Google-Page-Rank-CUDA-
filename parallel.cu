#include<stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <cuda.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include "file_read.h"
#include <time.h>

#include "svd.h"

#define debug(x) printf("debug %d \n",x);

double ** A;

double *d_AA;
double * H;
double * Q;
double * q;

double *h_AA;
double *h_H;
double *h_Q;
double *h_q;


long long nodes;
long long edges;

//--------------------------------------------//

__global__ void mal_A_Q_j(double *AA, int nodes, double* z, double* Q, int j , double *sum,int k)
{
    int id = threadIdx.x;
    if(id < nodes){
    	double temp = 0.0;
    	for(int ll =0;ll< nodes ; ll++)
    	{
    		//printf("inside mal_A_Q_j %f %f\n",AA[id*nodes+ll],Q[ll*(k+1)+j]);
    		temp = temp + AA[id*nodes+ll]*Q[ll*(k+1)+j]; 
    	}
        z[id] = temp;
    }
}

__global__ void norm_of_q(double* q, double *sum,int nodes)
{
	sum[0] = 0;
    for(int i=0; i<nodes; i++)
    {
        sum[0] += q[i]*q[i];
    }
	sum[0] = sqrt(sum[0]);
}

__global__ void first_col_Q(double *Q,double* q, double* sum, int nodes,int k)
{
    int id = threadIdx.x;
    if(id<nodes){
    	Q[id*(k+1)+0] = q[id]/sum[0];
    }
}

__global__ void dot_prod_assign_H(double* Q, double*z, int nodes, double* sum, double* H, int i, int j, int k)
{

    double temp=0.0;
    for(int ll=0; ll<nodes; ll++) 
    {
    	temp += Q[ll*(k+1)+i]*z[ll];
    }
    H[i*(k)+j]=temp;

    //printf("temp = %.64f",H[i*(k)+j]);
}

__global__ void update_z(double *z, double *Q, double *H, int nodes,int j,int i,int k)
{

    int id = threadIdx.x;
    if(id<nodes)
    {
    	//printf("\nQ bef 123 z= %.64f ,H= %.64f ,Q= %.64f\n",z[id],H[i*(k)+j],Q[id*(k+1)+i]);
        z[id] = z[id] - H[i*(k)+j]*Q[id*(k+1)+i];
        //printf("\nQ aft 321 z= %.64f ,H= %.64f ,Q= %.64f\n",z[id],H[i*(k)+j],Q[id*(k+1)+i]);
    }

    //for (int ll=0; ll<nodes; ll++)
    //{
    	//printf("\nQ bef 123 z= %.64f ,H= %.64f ,Q= %.64f\n",z[ll],H[i*(k)+j],Q[ll*(k+1)+i]);
      //  z[ll] = z[ll] - H[i*(k)+j]*Q[ll*(k+1)+i];
        //printf("\nQ aft 321 z= %.64f ,H= %.64f ,Q= %.64f\n",z[ll],H[i*(k)+j],Q[ll*(k+1)+i]);
    //}
}

__global__ void assign_norm_in_H(double *z, double *sum, int nodes, double *H, int j,int k)
{
    double temp = 0.0;
    for(int ll=0; ll<nodes; ll++)
    {
        temp = temp + z[ll]*z[ll];
    }
	H[(j+1)*k + j] = sqrt(temp);;
	//printf("norm z = %.64f\n",H[(j+1)*(k)+j] );
}

__global__ void set_Q(double *Q,double *H,int j, int k, int nodes,double *z)
{

    int id = threadIdx.x;
    if(id<nodes)
    {
       Q[id*(k+1)+(j+1)] = z[id]/H[(j+1)*(k)+j];
       //printf("Q he = %.64f / %.64f = %.64f\n",z[id],H[(j+1)*(k)+j],Q[id*(k+1)+(j+1)]);
    }

}

__global__ void new_H(double *H,int k,int nodes)
{
    int id = threadIdx.x;
	if(id<nodes)
	{
		H[id*k+id] -= 1;
	}
}

//--------------------------------------------//

__global__ void update_q(double* q, double* Q, int i, int col, int n)
{
    int id = threadIdx.x;
    if(id<n)
    {
        q[id] = Q[id*col + i];
    }
}




__global__ void norm_q(double* q, double* sum, int n)
{
    int id = threadIdx.x;
    if(id<n){
    	q[id] /= *sum;
    }
}

__global__ void get_v(double* VT, double* v, int k,int n)
{
    int id = threadIdx.x;
    if(id<n)
    {
		v[id] = VT[id*k +k-1+id];    
    }
}

__global__ void mul_Q_v(double *Q, double* v ,double* q, int n, int k)
{
    int id =  threadIdx.x;
    if(id < n){
        double temp = 0;
        for(int i = 0; i<k; i++)
        	temp += Q[id*(k+1)+i]*v[i];
    	
    	q[id] = temp;
    }
}

//----------------------------------------------//

double norm(double *abc,int k){
	// size_t n = sizeof(q) / sizeof(q[0]);
	double sum = 0.0;
	for(int i = 0; i < nodes; i++){
		sum += abc[i] * abc[i];
	}

	return sqrt(sum);
}

//-----------------------------------------------//

void Arnoldi(int k,double *q){

	double *sum ;
	cudaMalloc(&sum, sizeof(double)*1);

	// calculating norm of q
	norm_of_q<<<1, 1>>>(q,sum,nodes);
	cudaDeviceSynchronize();

	// Q[][1] = q/norm
	first_col_Q<<<1, nodes>>>(Q,q, sum, nodes,k);
	cudaDeviceSynchronize();
	
	double *z,*h_z;
	cudaMalloc(&z, sizeof(double) * nodes);
	h_z = (double *)malloc(sizeof(double)*nodes);

	for(int j = 0; j < k; j++)
	{

		// z=Aq[j]
		mal_A_Q_j<<<1, nodes>>>(d_AA, nodes, z,Q, j, sum,k);
		cudaDeviceSynchronize();

		cudaMemcpy(h_z,z, nodes*sizeof(double), cudaMemcpyDeviceToHost);

		//printf("before Z = ");
	    //for (int ll = 0;ll<nodes;ll++)
	    //{
	    //	printf("%.64f ",h_z[ll]);
	    //}
		//printf("\n");
		
		for(int i = 0; i <= j; i++)
		{
			//H[i,j] = Q[][i]_T * z
			dot_prod_assign_H<<< 1, 1>>>(Q, z, nodes, sum, H, i, j, k);
			cudaDeviceSynchronize();

			//z = z − h[i,j] * Q[i]
			update_z<<<1, nodes>>>(z, Q, H, nodes,j,i,k);
			
			cudaDeviceSynchronize();

			cudaMemcpy(h_z,z, nodes*sizeof(double), cudaMemcpyDeviceToHost);

			//printf("\nUP Z = ");
		    //for (int ll = 0;ll<nodes;ll++)
		    //{
		    	//printf("%.64f ",h_z[ll]);
		    //}
			//printf("\n");
		}
		//printf("\n");

		assign_norm_in_H<<<1, 1>>>(z,sum,nodes,H,j,k);
		cudaDeviceSynchronize();

		set_Q<<<1, nodes>>>(Q, H, j, k, nodes,z);
		cudaDeviceSynchronize();

	}
	cudaFree(z);
	cudaFree(sum);
	// debug(2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void call_arnoldi(int k){

	cudaMalloc(&H, sizeof(double) * (k + 1) * k);
	h_H = (double *)malloc(sizeof(double) *(k + 1)*(k));

	cudaMalloc(&Q, sizeof(double) * (nodes) * (k + 1));
	h_Q=(double *)malloc(sizeof(double)*(nodes)*(k+1));

	cudaMalloc(&q, sizeof(double)*nodes);
	h_q=(double *)malloc(sizeof(double)*(nodes));

	for(int i = 0; i < nodes; i++){
		h_q[i] = 1;
	}

	for(int i = 0; i < nodes*(k+1); i++){
		h_Q[i] = 0;
	}

	for(int i = 0; i < k*(k+1); i++){
		h_H[i] = 0;
	}

	cudaMemcpy(q, h_q, nodes*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Q, h_Q, nodes*(k+1)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(H, h_H, k*(k+1)*sizeof(double), cudaMemcpyHostToDevice);

//------------------------------------------------//

    // Instead of Convergence We used 50 iterations (Convergence takes a lot of time)

	for(int pp = 0; pp < 100; pp++)
	{
		
		Arnoldi(k,q);

//----------------- Printing Q H q ----------------//

		cudaMemcpy(h_q , q , sizeof(double)*nodes, cudaMemcpyDeviceToHost);

		cudaMemcpy(h_Q , Q , sizeof(double)*nodes*(k+1), cudaMemcpyDeviceToHost);

		cudaMemcpy(h_H , H , sizeof(double)*(k+1)*k, cudaMemcpyDeviceToHost);

		//printf("\nh_q\n");
		//for (int ll = 0;ll< k;ll++)
		//{
		//	printf("%lf , ",ll,h_q[ll]);
		//}

		//printf("\nQ\n");
		//for (int ll = 0;ll< nodes;ll++)
		//{
		//	for (int mm = 0;mm< k+1;mm++)
		//	{
		//		printf("%lf , ",h_Q[ll*(k + 1) + mm]);
		//	}
		//	printf("\n");
		//}

		//printf("\nH\n");
		//for (int ll = 0;ll< k+1;ll++)
		//{
		//	for (int mm = 0;mm< k;mm++)
		//	{
		//		printf("%lf , ",h_H[ll*k + mm]);
		//	}
		//	printf("\n");
		//}

//-------------------------------------------------//

		new_H<<<1,nodes>>>(H,k,nodes);

//------------------svd initials ----------------//
		
		cusolverDnHandle_t cusolverH = NULL;

	    const int rows = k+1;
	    const int cols = k;
	    const int lda = rows;

	    double *S = NULL;
	    double *U = NULL;
	    double *VT = NULL;
	    double *d_work = NULL;
	    double *d_rwork = NULL;
	    int *devInfo = NULL;
	    int lwork = 0;

	    cudaMalloc((void**)&S,sizeof(double)*cols);
	    cudaMalloc((void**)&U,sizeof(double)*lda*rows);
	    cudaMalloc((void**)&VT,sizeof(double)*lda*cols);
	    cudaMalloc((void**)&devInfo, sizeof(int));

	    cusolverDnCreate(&cusolverH);
	    cusolverDnDgesvd_bufferSize(cusolverH,rows,cols,&lwork);
	    cudaMalloc((void**)&d_work,sizeof(double)*lwork);

//-------------------- SVD ----------------------//
		
		signed char jobu = 'A';
	    signed char jobvt = 'A';
	    cusolverDnDgesvd (
	    cusolverH,
	    jobu,
	    jobvt,
	    rows,
	    cols,
	    H,
	    lda,
	    S,
	    U,
	    lda,
	    VT,
	    lda,
	    d_work,
	    lwork,
	    d_rwork,
	    devInfo);

		cudaDeviceSynchronize();

//------------------------------------------------//

		double * v;
		cudaMalloc(&v , k*sizeof(double));

		get_v<<<1, cols>>>(VT, v, cols, nodes);

		mul_Q_v<<<(nodes)/1024 + 1, 1024>>>(Q, v, q, nodes, cols);

		//debug(16);
		//printf("%d  \n",pp);
		cudaFree(v);
		cudaFree(S);
		cudaFree(U);
		cudaFree(VT);


	}
	cudaFree(H);
	cudaFree(Q);
	cudaFree(q);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int main()
{
	
	nodes = 1000;

	A = (double **)malloc(sizeof(double *) * (nodes));
    for(long long i = 0; i < nodes; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * (nodes));
 	}

//--------------------------------------------------------//

	double damping_factor = 0.85;
	FILE *fptr = fopen("innn.txt", "r");

	if (fptr == NULL){

        printf("Cannot open file \n");

        exit(0);
        debug(1);
    }
    else{

    	char line [128]; /* or other suitable maximum line size */

    	fgets ( line, sizeof line, fptr);
        long long nodes = atoi(line);

        fgets ( line, sizeof line, fptr );
        long long edges = atoi(line);

     	while ( fgets ( line, sizeof line, fptr ) != NULL ) /* read a line */
    	{
            // debug(2);
    	    char* token = strtok(line, " "); 
            
            long long a,b;
            long long k = 0;
            while (token != NULL) { 
                // printf("%s",token);
                if(k==0)
                {
                    a = atoi(token);
                    // printf("%d  ", a); 
                    token = strtok(NULL, " "); 
                }
                else
                {
                    b = atoi(token);
                    // printf("%d  ", b); 
                    token = strtok(NULL, " "); 
                }
                k++;
                
                // count[a]++;
            } 
            A[a][b]=1;
      	}
        double jump = (1 - damping_factor) / (double)nodes;
        //printf("jump = %f\n",jump);
        for(long long i = 0; i < nodes; i++){
            long long val = 0;
            for(long long j = 0; j < nodes; j++){
                if((A[i][j] == 1.0)){
                    val++;
                }
            }
            // long long val = count[i];
            // printf("%d ... ", count);
            // printf("%d\n", val);

            for(long long j = 0; j < nodes ; j++){
                if(val == 0){
                    val = 1;
                }
                A[i][j] = A[i][j]*damping_factor / val + jump;
                // printf("%f ",A[i][j]);
            }
        }
      	fclose (fptr);

    }
    printf("-------End of read------- \n\n");

//--------------------------------------------------------//

    //A transpose

    for (int ll = 0;ll< nodes;ll++)
	{
		for (int mm = 0;mm< ll;mm++)
		{
			double temp = A[ll][mm];
			A[ll][mm] = A[mm][ll];
			A[mm][ll] = temp;
		}
	}

	//printf("A transpose matrix : \n");
	//for (int ll = 0;ll< nodes;ll++)
	//{
	//	for (int mm = 0;mm< nodes;mm++)
	//	{
	//		printf("%lf , ",A[ll][mm]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");

	h_AA = (double *)malloc(sizeof(double) *(nodes*nodes));

	cudaMalloc(&d_AA,sizeof(double) * (nodes) *  nodes);

	int count = -1;
    for(int ll = 0; ll < nodes; ll++){
    	for(int mm = 0; mm < nodes; mm++){
        	h_AA[++count] = A[ll][mm];
    	}
    }

    cudaMemcpy(d_AA, h_AA, nodes*nodes*sizeof(double),cudaMemcpyHostToDevice);

    int k=8;

    double total_time;
	clock_t start, end;
	start = clock();

    call_arnoldi(k);

    end = clock();
	total_time = ((double) (end - start)) / CLOCK_TAI;

    printf("\n\nFinal ANS =\n");
    double tempp = 0.0;
    for(int ppp=0;ppp<nodes;ppp++)
    {
    	tempp = tempp + h_q[ppp]*h_q[ppp];
    }
    tempp = sqrt(tempp);
    for(int ppp=0;ppp<nodes;ppp++)
    {
    	h_q[ppp]=h_q[ppp]/tempp;
    	printf("%d - %f\n",ppp,h_q[ppp]);
    }

	printf("\nTime taken for k = %d is: %f\n",k ,total_time);

    free(h_AA);
    free(A);
    cudaFree(d_AA);

}