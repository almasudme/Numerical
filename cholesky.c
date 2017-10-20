
#include<stdio.h>
#include<math.h>
#define size 4

int cholesky(){	
	float A[size][size]= {{16,4,4,-4},{4,10,4,2},{4,4,6,-2},{-4,2,-2,4}};
	float L[size][size]={0};
	float LT[size][size]={0};
	float b[size]={32,26,20,-6};
	int i,j,k;
	float x[size]={0};
	float z[size]={0};
	float temp;
	float norm_x;
	float norm_xbar;
	float x_accurate[size]={0};
/*Ill conditioned matrix*/
	 //Commented out for running the tridiagonal matrix
	for (i=0;i<size;i=i+1){
		for (j=0;j<size;j=j+1){
			if (i+j != 1){
			A[i][j] = 1.0/((double)(i)+(double)(j)-1.0);
			}
		}
	}
	
for (i=0;i<size;i=i+1){
        if (i == 0){
            A[i][i]=4;
            A[i][i+1]=1;
            b[i]=5;
			x_accurate[i]=1;
        } else if (i==size-1){
            A[i][i]=4;
            A[i][i-1]=1;
            b[i]=5;
			x_accurate[i]=1;
        } else {
            A[i][i]=4;
            A[i][i-1]=1;
            A[i][i+1]=1;
            b[i]=6;
			x_accurate[i]=1;
        }


};


/*Implementation of Cholesky algorithm*/
for(i=0;i<size;i++){
	for(j=0;j<=i;j++){ 
		temp=0;
		if(j==i){
			for(k=0;k<j;k++){
				temp = temp +(L[j][k])*(L[j][k]);
			}
			L[j][j]=(sqrt((A[j][j]-temp)));
		}else {
			for(k=0;k<j;k++){
				temp=temp +((L[i][k]*L[j][k]));
			}
			//L[i][j]=(1/(L[j][j]))*((A[i][j])-temp);
			L[i][j]=((A[i][j])-temp)/(L[j][j]);
		}

		if(i<j){
			L[i][j]=0;
		}
	}
}

/*LT = Transpose of L*/
for(i=0;i<size;i++){
	for(j=0;j<size;j++){
		LT[i][j]=L[j][i];
	}
}
temp = 0;

/* forward substitutionLz=b*/
    for(i=0; i<size; i++){
		z[i]=b[i];
		for(j=0; j<i; j++){
				z[i]= z[i]-L[i][j]*z[j];
		}
		z[i] = z[i]/LT[i][i];
    }

/* backward substitution LTx=z*/

for(i=size-1; i>=0; i--)    {
        x[i]= z[i];
        for(j=i+1; j<size; j++)
        {
            x[i]= x[i]-LT[i][j]*x[j];
        }
        x[i] = x[i]/LT[i][i];
    }

/*Norm*/
temp = 0;
for(i=0;i<size;i++){
	temp = temp+(x_accurate[i] -x[i])*(x_accurate[i] -x[i]);	
}
norm_xbar = sqrt(temp);
temp=0;
for(i=0;i<size;i++){
	temp = temp+x[i]*x[i];	
}
norm_x = sqrt(temp);

printf("||x - xbar||2/||x2|| for n = %d : %f \n", size, norm_xbar/norm_x);

return 0;
}

