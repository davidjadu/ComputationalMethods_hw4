#include<stdio.h>
#include<stdlib.h>
#include<math.h>


float *load_matrix(char *filename, int n);
float *multiply(float *m1,float *m2,int n_row_1, int n_col_1, int n_row_2,int n_col_2);
float *traspose(float *A, int n_row,int n_col);
void lu_decomposition(float *A, float *b, float *U, float *L, int n);
void make_data(float *file_data, float *data, float *b,float *t, int n_row);
void solve_upper_triangular(float *matrix, float *b, float *y, float *v, float *g);
void print_in_file(float y, float v, float g);

int main(){
  char *filename = "data_F.dat";
  int n_col=3;
  int n_row = 884;
  float *file_data = load_matrix(filename, n_row);
  float *data = malloc(n_row*n_col*sizeof(float));
  float *b = malloc(n_row*sizeof(float));
  float *theta = malloc(n_row*sizeof(float));
  make_data(file_data, data, b,theta, n_row);
  float *data_traspose = traspose(data,n_row,n_col);
  float *matrix = multiply(data_traspose,data, n_col,n_row,n_row,n_col);
  float *new_b = multiply(data_traspose,b, n_col, n_row, n_row, 1);
  float *U = malloc(n_col*n_col*sizeof(float));
  float *L = malloc(n_col*n_col*sizeof(float));
  lu_decomposition(matrix,new_b,U,L,n_col);
  float a1,a2,a3;
  solve_upper_triangular(U,new_b,&a1,&a2,&a3);
  print_in_file(a1,a2,a3);
  return 0;
}

/* function that puts the data in the desired form
 */
void make_data(float *file_data, float *data, float *b, float *t, int n_row){
  int i,j;
  for (i=0; i<n_row; i++){
      data[i*3]=1;
      t[i]=file_data[2*i];
      data[i*3+1]=t[i];
      data[i*3+2]=t[i]*t[i];
      b[i]=file_data[2*i+1];
  }
}


/* Function that loads a matrix into the program
 */
float *load_matrix(char *filename, int n){
  float *matrix;
  FILE *in;
  int n_row=n, n_cols=2;
  int i;
  int j;

  if(!(in=fopen(filename, "r"))){
    printf("Problem opening file %s\n", filename);
    exit(1);
  }

  matrix = malloc(n_row*n_cols*sizeof(float));

  for (i=0; i<n_row; i++){
    for(j=0;j<n_cols;j++){
      fscanf(in, "%f", &matrix[i*n_cols+j]);
    }
  }

  fclose(in);
  return matrix;
}


/* Function that multiplies two matrices
 */
float *multiply(float *m1,float *m2,int n_row_1, int n_col_1, int n_row_2,int n_col_2){
  if (n_col_1!=n_row_2){
    printf("There was an attempt to multiply two matrices that can't be multiply");
    exit(1);
  }
  float *m;
  int i,j,k, n = n_row_1, p=n_col_2;
  m= malloc(n*p*sizeof(float));
  float actual;
  for(i=0;i<n;i++){
    for (j=0;j<p;j++){
      actual=0;
      for (k=0;k<n_col_1;k++)
	actual+=m1[n_col_1*i+k]*m2[n_col_2*k+j];
    m[p*i+j]=actual; 
    }
  }  
  return m;
}

/* Function that trasposes the matrix A
 */
float *traspose(float *A, int n_row,int n_col){
  float *m = malloc(n_row*n_col*sizeof(float));
  int i, j;
  for (i=0; i<n_row; i++)
    for(j=0; j<n_col; j++)
      m[j*n_row+i]=A[i*n_col+j];
  return m;
}

/* Function that makes the LU decompsoiton, returns the answer by parameter
 */
void lu_decomposition(float *A, float *b,float *U, float *L, int n){
  int i,j,k;
  int found = 0;
  float temp;
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      U[i*n+j]=A[i*n+j];
    }
  }
  for (i=0; i<n; i++){
    L[i*n+i]=1;
    if (U[i*n+i]<pow(1.0/10.0,10)){ // if a(i,i)=0, change lines
      for (j=i+1; j<n && found==0 ; j++){
	if (U[j*n+i] !=0){
	  U[i*n+i]=U[j*n+i];
	  temp = b[i];
	  b[i]=b[j];
	  b[j]=temp;
	  found = 1;
	  for (k=i+1; k<n; k++){ //Swap the two lines
	    temp = U[i*n+k];
	    U[i*n+k]=U[j*n+k];
	    U[j*n+k]=temp;
	  }
	}
      }
      if(found == 0){
	printf("There was an attempt to lu_decompose a singular matrix \n");
	exit(1);
      }
    }
    for(j=i+1;j<n;j++){
	L[j*n+i]=U[j*n+i]/U[i*n+i];
	U[j*n+i]=0;
	b[j]=b[j]-L[j*n+i]*b[i];
	for(k=i+1;k<n; k++){
	  U[j*n+k]=U[j*n+k]-L[j*n+i]*U[i*n+k];
	}
    }
  }
}

/* Función que resuelve un sistema 3x3 que está en forma triangular superior
 */
void solve_upper_triangular(float *matrix, float *b, float *y, float *v, float *g){
  int n=3;
  float *x= malloc(n*sizeof(float));
  x[n-1]=b[n-1]/matrix[n*n-1];

  int i,j;
  float dot;
  for (i=n-2;i>-1;i--){
    dot =0;
    for (j=n-1; j>i;j--){
      dot+=x[j]*matrix[n*i+j];
    }
    x[i]=(b[i]-dot)/(matrix[n*i+i]);
  }
  *y=x[0];
  *v=x[1];
  *g=x[2];
}

/*Función encargada de imprimir los datos en el archivo
 */
void print_in_file(float y, float v, float g){
  char *filename_out="constants_F_data_2.dat";
  FILE *out=fopen(filename_out,"w");
  fprintf(out, "%f %f %f \n", y, v, g);
  fclose(out);
}

