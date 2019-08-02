#include <iostream>
#include <cmath>
using namespace std;

double determinant(double** a, int n){
	double det=0;
	int l, m, temp;
	if(n==1){ return a[1][1];}
	
	if(n==2){
		return a[1][1]*a[2][2] - a[1][2]*a[2][1];
	}

    double **b = new double*[n];
    for(int i=0; i<n; i++)
        b[i] = new double[n];
    
	
	for(int i=1; i<=n; i++){
		l =1; m =1;
		for(int j=2; j<=n; j++){
			m = 1;
			for(int k=1; k<=n; k++){
				if(k!=i) {
					b[l][m] = a[j][k];
					m ++;
				}
			}
			l++;
		}
		
        
		for(int w=1; w<n;w++){
		
			for(int x=1; x<n; x++)
				cout<<b[w][x]<<" ";
				
			cout<<"\n";
		}
		if(i%2) temp = 1;
		else temp = -1;
		
		//cout<<determinant(b, n-1)<<endl;

		det += determinant(b, n-1) * temp * a[1][i];
	}
    for(int i=0; i<n; i++)
        delete [] b[i];
    delete [] b;
	return det;
}

double norm(double *a, int m){
	double temp = 0;
	for(int i=1; i<=m; i++)
		temp += a[m]*a[m];
	return sqrt(temp);
}

void jacobi(double** a, double* b, double* XO, double tol, int maxiter, int n){
	int k=1; double sum;
	double *x = new double[n+1], *y = new double[n+1];
	while(k<=maxiter){
		for(int i=1; i<=n; i++){
			sum = 0;
			for(int j=1; j<=n; j++)
				if(i != j)
					sum -= a[i][j]*XO[j];		
			sum += b[i];
			x[i] = sum/a[i][i];
		}
		
		for(int i=1; i<=n; i++)
			y[i] = x[i] - XO[i];
		
		if(norm(y, n) < tol){
		 	for(int i=1; i<=n; i++)
				XO[i] = x[i];
		 	break;
		}
		
		k++;
		
		for(int i=1; i<=n; i++)
			XO[i] = x[i];	
	}
}

int main(){
	int n;
	cout<<"Enter the dimension of the matrix:";
	cin>>n;
	double **array = new double*[n+1], *b = new double[n+1], *x = new double[n+1];
	for(int i=0; i<n+1; i++)
		array[i] = new double[n+1];
	
	cout<<"Enter the elements of the matrix A row-wise"<<endl;
	for(int i=1; i<n+1; i++)
		for(int j=1; j<n+1; j++)
			cin>>array[i][j];
	
	cout<<"Enter the elements of the vector b row-wise"<<endl;
	for(int j=1; j<n+1; j++)
		cin>>b[j];
	
	for(int i=1; i<n+1; i++)
		x[i] = 0;
	
	if(determinant(array, n) != 0) {
        cout<<"Unique solution exists"<<endl;
        jacobi(array, b, x, 0.01, 100, n);
        for(int i=1; i<=n; i++)
            cout<<x[i]<<endl; 
    
    }
	else cout << "Unique solution doesn't exists"<<endl;
			
	//cout << determinant(array, n);

	return 0;
}