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
		
        
		/**for(int w=1; w<n;w++){
		
			for(int x=1; x<n; x++)
				cout<<b[w][x]<<" ";
				
			cout<<"\n";
		}**/
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

void guass_jordan_elimination(double **A, double *x, double *b, double **A_inv, int n)
//Here A_inv should be intialized as identity matrix
{
    int i, j, k;
    double temp, *pt;
    for(i=1; i<=n; i++){
        x[i] = 0;
        for(j=1; j<=n; j++)
            if(i!=j) A_inv[i][j] = 0;
            else A_inv[i][j] = 1;
    }
    for(i=1; i<=n; i++){
        for(j=1; j<i; j++){
            temp = A[i][j];
            for(k=1; k<=n; k++){

                A[i][k] = A[i][k] - A[j][k]*temp;
                A_inv[i][k] = A_inv[i][k] - A_inv[j][k]*temp;
            }
        }
        temp = A[i][i];
        if(temp!=0)
        for(j=1; j<=n; j++){ 
            A[i][j] = A[i][j]/temp;
            A_inv[i][j] = A_inv[i][j]/temp;
        }
        else{
            pt = A[i];
            A[i] = A[n];
            A[n] = pt;
            i--;
            continue;
        }
        for(j=i-1; j>0; j--){
            temp = A[j][i];
            for(k=1; k<=n; k++){
                A[j][k] = A[j][k] - A[i][k]*temp;
                A_inv[j][k] = A_inv[j][k] - A_inv[i][k]*temp;
            }
        }
    }
    for(i=1; i<=n; i++)
        for(j=1; j<=n; j++)
            x[i] += A_inv[i][j]*b[j];
}

double norm(double *a, int m){
	double temp = 0;
	for(int i=1; i<=m; i++)
		temp += a[m]*a[m];
	return sqrt(temp);
}

int factorial(int n){
	if(n==0) return 1;
	if(n>0){
		int product=1, i=1;
		for(i=1; i<=n; i++)
			product = product*i;
		return product;
	}
	else return 0;
}

int combinations(int n, int m){
	return factorial(n)/(factorial(m)*factorial(n-m));	
}

void jacobi(double** a, double* b, double* XO, double tol, int maxiter, int n){
/**
    This method only works best when the matrix is diagonally dominant.
**/
	int k=1; double sum, *temp_ptr, temp;
	double *x = new double[n+1](), *y = new double[n+1]();
	while(k<=maxiter){
		for(int i=1; i<=n; i++){
			sum = 0;
			for(int j=1; j<=n; j++)
				if(i != j)
					sum -= a[i][j]*XO[j];		
			sum += b[i];
			if(a[i][i]!=0){
				x[i] = sum/a[i][i];
			}
			else{
				temp_ptr=a[i];
				a[i] = a[n-i+1];
				a[n-i+1] = temp_ptr;
				temp = b[i];
				b[i] = b[n-i+1];
				b[n-i+1] = temp;
				k--;
			}
		}
		
		for(int i=1; i<=n; i++)
			y[i] = x[i] - XO[i];
		
		//cout<<"Norm: "<<norm(y, n)<<endl;
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

void generate_combination(int** a, int* b, int reqLen, int start, int currLen, bool check[], int len, int* k){
	if(currLen > reqLen)
		return;
	else if (currLen == reqLen) 
	{
		int e=1;
		for (int i = 1; i <= len; i++) 
		{
			if (check[i] == true) 
			{
				a[*k][e] = b[i];
				e++;
			}
		}
		(*k)++;
		return;
	}

	if (start == len) 
		return;
	
	check[start+1] = true;
	generate_combination(a, b, reqLen, start+1, currLen+1, check, len, k);
	check[start+1] = false;
	generate_combination(a, b, reqLen, start+1, currLen, check, len, k);
	/**int comb1;
	for(int i=w; i<(w+n-m+x); i++){
		comb1 = combinations(n-x-i, m-x-i);
		for(int j=x; j<=comb1; j++)
			a[j][w] = b[k];
		generate_combinations(a, b, n-i, m-i, w+comb1, x+1, k+1);
		k++;
	}**/
}


void print_array(int** a, int m, int n){
	for(int i=1; i<=m; i++){
		for(int j=1; j<=n; j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}
}

void print_array(double** a, int m, int n){
	for(int i=1; i<=m; i++){
		for(int j=1; j<=n; j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}
}

int main(){
	int m, n, number_combinations, i, j, k=1, y, z, w, t, max=1, number_slack;
    char flag;
	cout<<"Enter the dimension of the matrix:";
	cin>>m>>n;
	double **array = new double*[n+1], **array3 = new double*[m+1], **array3_inv = new double*[m+1], *b = new double[m+1], *x = new double[m+1], *c = new double[n+1], a=0, temp_value=0, max_value=0;
	for(i=0; i<n+1; i++)
	{
		array[i] = new double[m+1];
        if(i<m+1){
		  array3[i] = new double[m+1];
		  array3_inv[i] = new double[m+1];
        }
	}
	
	number_combinations = combinations(n,m);
	int **array2 = new int*[number_combinations+1], *b2 = new int[n+1];    
	double **array_sol = new double*[number_combinations+1];
    bool* check = new bool[n+1](), flag_feasible = true;

	for(i=0; i<=n; i++)
		b2[i] = i;
	for(i=0; i<number_combinations+1; i++){
		array2[i] = new int[m+1];
        array_sol[i] = new double[n+1]();
    }
	generate_combination(array2, b2, m, 0, 0, check, n, &k);
	//print_array(array2, number_combinations, m);
	
	cout<<"Enter the elements of the matrix A row-wise"<<endl;
	for(j=1; j<m+1; j++)
		for(i=1; i<n+1; i++)	
			cin>>array[i][j];
	cout<<"Enter the elements of the vector b row-wise"<<endl;
	for(j=1; j<m+1; j++)
		cin>>b[j];
    cout<<"Do you want to maximize some objective function?(Y/N):";	
        cin>>flag;
    if(flag=='Y'){
        cout<<"Enter n cost coeffticients row-wise"<<endl;
        for(i=1; i<n+1; i++)
            cin>>c[i];
        cout<<"Enter number of slack variables: ";
            cin>>number_slack;
    }
    cout<<endl;
	for(i=1; i<n+1; i++)
		x[i] = 0;
	
	for(i=1; i<=number_combinations; i++){
		for(j=1; j<=m; j++)
			for(k=1; k<=m; k++){
				array3[j][k] = array[array2[i][k]][j];
			}

		for(z=1; z<n+1; z++)
			x[z] = 0;

		if(determinant(array3, m) != 0) {
           cout << "Found  Non-Degenerate solution on selecting columns:";
           for(z=1; z<=m; z++){
                cout<<array2[i][z]<<" ";
           }
           cout<<endl;
        	//jacobi(array3, b, x, 0.01, 1000, m);
        	guass_jordan_elimination(array3, x, b, array3_inv, m);

        	t=1;
 
        	cout<<"\"(";
 
        	for(z=1; z<=m; z++){
        		if(z==1){
        			for(y=1; y<array2[i][1]; y++){
        				cout<<"0, ";
                        array_sol[i][t] = 0;
        				t++;
        			}
        		}
        		else if((w=array2[i][z]-array2[i][z-1])>1){
        			for(y=1; y<w; y++){
          				cout<<"0, ";
                        array_sol[i][t] = 0;
        				t++;
        			}
        		}
                //if(t<=n-number_slack)
                //    temp_value += c[t]*x[z];
            	cout<<x[z]<<", ";
                array_sol[i][t] = x[z];
            	t++;
        	}
        	for(;t<=n;t++){
        		cout<<"0, ";
                array_sol[i][t] = 0;
        	}
            cout<<")\n";
            if(flag=='Y'){
            flag_feasible = true; z=1;
            for(t=1; t<=n; t++){
                if(array_sol[i][t]<0){
                    cout<<"But the solution isn't feasible";
                    flag_feasible = false;
                    break;
                }
            }
            if(flag_feasible==true){    
            for(z=1; z<=m; z++){
                temp_value = 0;
                for(t=1; t<=n-number_slack; t++)
                    temp_value += array[t][z]*array_sol[i][t];
                if(temp_value>b[z]+0.02){
                    cout<<temp_value<<" "<<z<<endl;
                    cout<<"But the solution isn't feasible";
                    break;
                }
            }}
            if(z==m+1){
                cout<<"The solution is feasible";
                temp_value=0;
                for(t=1; t<=n-number_slack; t++)
                    temp_value += array_sol[i][t]*c[t];
                if(temp_value>max_value){
                    max_value = temp_value;
                    max = i;
                }
            }
            }
        	cout<<"\n\n";
    	}
		else{
           cout << "Found  Degenerate solution on selecting columns:";
           for(z=1; z<=m; z++){
                cout<<array2[i][z]<<" ";
           }
           cout<<endl<<endl;
        }
	}
    if(flag=='Y'){
        cout<<"\nThe optimal solution is\n(";
        for(t=1; t<=n; t++)
            cout<<array_sol[max][t]<<" ";
        cout<<")\n";
        cout<<"The maximized value of function is:"<<max_value<<endl;
    }
    /**for(i=0; i<n+1; i++)
    {
        delete [] array[i];
        if(i<m+1){
            delete [] array3[i];
            delete [] array3_inv[i];
        }
    }
    for(i=0; i<number_combinations+1; i++)
        delete [] array2[i];
    delete [] array;
    delete [] array2;
    delete [] array3;
    delete [] array3_inv;
    delete [] b;
    delete [] b2;
    delete [] check;
    **/
	//cout << determinant(array, n);
	return 0;
}