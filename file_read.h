#include<stdio.h>
#include <string.h> 
#include <stdlib.h>

#define debug(x) printf("debug %d \n",x);

void read_A(double ** A){

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
    	// printf("%lld",nodes);
        // printf("\n");
        // debug(1);
        fgets ( line, sizeof line, fptr );
        long long edges = atoi(line);
        // printf("%lld",edges);
        // printf("\n");

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
            // printf("\n");
            // i++;
            // fputs ( line, stdout ); /* write the line */
      	}
        double jump = (1 - damping_factor) / (double)nodes;
        // printf("jump = %f\n",jump);
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

}
