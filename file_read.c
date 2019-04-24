#include<stdio.h>
#include <string.h> 

#define debug(x) printf("debug %d \n",x);

void main(){

    double damping_factor = 0.85;
	FILE *fptr = fopen("input2.txt", "r");

	if (fptr == NULL){

        printf("Cannot open file \n");

        exit(0);

    }
    else{

    	char line [ 128 ]; /* or other suitable maximum line size */

    	fgets ( line, sizeof line, fptr );
        long long nodes = atoi(line);
    	printf("%d",nodes);
        printf("\n");
        // debug(1);
        fgets ( line, sizeof line, fptr );
        long long edges = atoi(line);
        printf("%d",edges);
        printf("\n");

        double **A = (double **)malloc(sizeof(double *) * (nodes + 1));
        // long long *count = (long long *)malloc(sizeof(long long) *(nodes + 1));

        // for(long long i = 0; i < nodes + 1; i++){
        //     count[i] = 0;
        // }
        for(long long i = 0; i < nodes + 1; i++){
            A[i] = (double *)malloc(sizeof(double) * (nodes + 1));
        }

        // debug(1);
        // long long i = 0;
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
                    printf("%d  ", a); 
                    token = strtok(NULL, " "); 
                }
                else
                {
                    b = atoi(token);
                    printf("%d  ", b); 
                    token = strtok(NULL, " "); 
                }
                k++;
                A[a][b]=1;
                // count[a]++;
            } 
            printf("\n");
            // i++;
            // fputs ( line, stdout ); /* write the line */
      	}
        double jump = (1 - damping_factor) / (double)nodes;
        printf("%f\n",jump);
        for(long long i = 0; i < 10; i++){
            long long val = 0;
            for(long long j = 0; j < nodes + 1 ; j++){
                if((A[i][j] == 1.0)){
                    val++;
                }
            }
            // long long val = count[i];
            // printf("%d ... ", count);
            printf("%d\n", val);

            for(long long j = 0; j < 10 ; j++){
                if(val == 0){
                    val = 1;
                }
                A[i][j] = A[i][j]*damping_factor / val + jump;
                printf("%f ",A[i][j]);
            }
            printf("\n");
        }
      	fclose (fptr );
    }
}