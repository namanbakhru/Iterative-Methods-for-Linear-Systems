#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

/**
    In case of memory allocation problem keep re-executing until the memory gets allocated.
**/
void print_simplex_table(double **A, int m, int n, double *b, double *c, int *basic_ind, double *delta){
	int i, j;
    double cost=0;
    cout<<"|-------------------------------------------------------------------------------------------"<<endl;
	cout<<"|B    CB      ";
	for(j=1; j<=m+n; j++)
		cout<<"x"<<j<<"       ";
	cout<<"  b     "<<endl;
    cout<<"|-------------------------------------------------------------------------------------------"<<endl;
    cout.precision(5);
	for(i=1; i<=m; i++){
		cout<<"|"<<setw(5)<<left<<basic_ind[i]<<setw(5)<<left<<c[basic_ind[i]]<<"   ";
		for(j=1; j<=m+n; j++)
			cout<<setw(8)<<A[i][j]<<" ";
		cout<<"  "<<setw(8)<<b[i]<<endl;	
	}
	cout<<"|------------------------------------------------------------------------------------------"<<endl;
	for(j=1; j<=m; j++)
		cost = cost + c[basic_ind[i]]*b[i];
	cout<<"|"<<setw(12)<<" ";
	for(j=1; j<=n+m; j++)
		cout<<setw(8)<<delta[j]<<" ";
    cout<<endl<<"|------------------------------------------------------------------------------------------"<<endl;
    return;
}

int main(){
	int m, n, i, j, number_iterations=0, q, p, *basic_ind = new int[m], *non_basic_ind = new int[n], number_artificial;
	bool opt_flag = false;
    cout<<"Enter the values of m and n:";
	cin>>m>>n;
	double **array = new double*[m], *b = new double[m], *delta = new double[m+n], *c = new double[n+m], temp_cost, min=1, temp=0, big_M;
	array = array - 1;
	b = b - 1;
	delta = delta - 1;
	c = c - 1;
	basic_ind -= 1;
    non_basic_ind -= 1;
	for(i=1; i<=m; i++){
		array[i] = new double[n+m];
		array[i] -= 1;
	}

	cout<<"Enter the elements of the matrix A row-wise"<<endl;
	for(i=1; i<=m; i++){
		for(j=1; j<=n; j++)
			cin>>array[i][j];
		for(; j<=n+m; j++)
			array[i][j] = 0;
		array[i][n+i] = 1;
	}
	for(i=1; i<=m; i++)
		basic_ind[i] = n+i;
    for(i=1; i<=n; i++)
        non_basic_ind[i] = i; 		
	cout<<"Enter the elements of vector b row-wise"<<endl;
	for(i=1; i<=m; i++)
		cin>>b[i];
	
	cout<<"Enter the number of artificial variables"<<endl;
		cin>>number_artificial;
	
	cout<<"Enter the n+m cost coefficients for phase-I: ";
	for(i=1; i<=n+m; i++)
        cin>>c[i];


	//print_array(array, m, m+n);
    while(true){
        for(j=1; j<=m+n; j++){
	        temp_cost=0;
	        for(i=1; i<=m; i++)
	            temp_cost += c[basic_ind[i]]*array[i][j];
	        delta[j] = temp_cost-c[j];			
	    }
        print_simplex_table(array, m, n, b, c, basic_ind, delta);
        q = -1; p =-1;
        min = 0;
        for(i=1; i<=n; i++){
            if(delta[non_basic_ind[i]]<min){
                min = delta[non_basic_ind[i]];
                q = non_basic_ind[i];
            }
        }
        min = 100000;
        if(q == -1){
            cout<<"The solution is optimal"<<endl;
            opt_flag = true;
            break;
        }

        for(i=1; i<=m; i++){
            if(array[i][q]>0){
                if(b[i]/array[i][q]<min){
                    min = b[i]/array[i][q];
                    p = i;
                }
            }
        }
        if(p == -1){
            cout<<"The LPP has unbounded solution"<<endl;
            break;
        }
        cout<<endl<<endl;
        temp = array[p][q];
        for(i=1; i<=n+m; i++)
            array[p][i] = array[p][i]/temp;
        b[p] = b[p]/temp;
        for(j=1; j<=m; j++){
            if(j==p) continue;
            temp = array[j][q];
            for(i=1; i<=m+n; i++)
                array[j][i] = array[j][i] - array[p][i]*temp;
            b[j] = b[j]- b[p]*temp;
        }
	    //basic_ind[p] qNB
	    for(i=1; i<=n; i++)
	        if(non_basic_ind[i] == q)
                break;
        non_basic_ind[i] = basic_ind[p];
	    basic_ind[p] = q;
        number_iterations++;
    }
    //print_simplex_table(array, m, n, b, c, basic_ind, delta);
    cout<<"Number of iterations algorithm took for the given problem:"<<number_iterations<<endl;
    if(opt_flag == true){
		cout<<"Enter the n+m (excluding the number of artificial variables) cost coefficients for phase-II: ";
		for(i=1; i<=n+m-number_artificial; i++)
        	cin>>c[i];
		while(true){
        	for(j=1; j<=m+n-number_artificial; j++){
	        	temp_cost=0;
	        	for(i=1; i<=m; i++)
	            	temp_cost += c[basic_ind[i]]*array[i][j];
	        	delta[j] = temp_cost-c[j];			
	    	}
        	print_simplex_table(array, m, n, b, c, basic_ind, delta);
        	q = -1; p =-1;
        	min = 0;
        	for(i=1; i<=n; i++){
            	if(delta[non_basic_ind[i]]<min){
            	    min = delta[non_basic_ind[i]];
            	    q = non_basic_ind[i];
            	}
        	}
        	min = 100000;
        	if(q == -1){
            	cout<<"The solution is optimal"<<endl;
            	opt_flag = true;
            	break;
        	}

        	for(i=1; i<=m; i++){
            	if(array[i][q]>0){
                	if(b[i]/array[i][q]<min){
                    	min = b[i]/array[i][q];
                    	p = i;
                }
            }
        	}
       	 	if(p == -1){
            	cout<<"The LPP has unbounded solution"<<endl;
            	break;
        	}
        	cout<<endl<<endl;
        	temp = array[p][q];
        	for(i=1; i<=n+m-number_artificial; i++)
            	array[p][i] = array[p][i]/temp;
        	b[p] = b[p]/temp;
        	for(j=1; j<=m; j++){
            	if(j==p) continue;
            	temp = array[j][q];
            	for(i=1; i<=m+n-number_artificial; i++)
            	    array[j][i] = array[j][i] - array[p][i]*temp;
            	b[j] = b[j]- b[p]*temp;
        	}
	    	//basic_ind[p] qNB
	    	for(i=1; i<=n; i++)
	        	if(non_basic_ind[i] == q)
                	break;
        	non_basic_ind[i] = basic_ind[p];
	    	basic_ind[p] = q;
        	number_iterations++;
    	}
	/**
        for(i=1; i<=n; i++)
            delta[i] = 0;
        cout<<endl;
        cout<<"The solution is:";
        temp_cost=0;
        for(i=1; i<=m; i++)
            temp_cost += c[basic_ind[i]]*b[i];            
        for(i=1; i<=m; i++)
            if(basic_ind[i]<=n)
                delta[basic_ind[i]] = b[i];
        for(i=1; i<=n; i++)
            cout<<" x"<<i<<"="<<delta[i];
        cout<<" Z="<<temp_cost<<endl;
        **/
    }
    for(i=1; i<=m; i++){
        array[i] += 1;
        delete [] array[i];
    }
    array = array + 1;
    b = b + 1;
    delta = delta + 1;
    c = c + 1;
    basic_ind += 1;
    non_basic_ind += 1;
    delete [] c;
    delete [] delta;
    delete [] b;
    delete [] array;
    delete [] non_basic_ind;
    delete [] basic_ind;
	return 0;
}
/** Programs to be mailed to:
	operationsresearch18@gmail.com	
**/
