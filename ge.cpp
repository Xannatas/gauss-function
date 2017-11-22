/* Simple linear system solver, from https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/#c-code */

/* Some modifications by me.  Albert Schueller*/

/* Isaiah Banta , Austin Betschart */

/* R-4.7 , R-4.10 , R-4.13 , R-4.21, R-4.24 , R-4.29 */

/* Isaiah Banta and Austin Betschart */

#include <iostream>
#include <fstream>
#include <cmath> // provides pow()
#include <sys/time.h> // provides gettimeofday() (gives milliseconds)
#include <vector> // provides the vector class
#include <cstdlib> // provides the rand() function
#include <ctime> // provides the time() function
#include <algorithm> // provides the sort() function
using namespace std;

/*
My Above and Beyond is my main function, where I put the whole Gauss function in a for loop and then input each of the recieved values into a CSV file which I can then import into Google Sheets to get an immediate reading on my values. This allowed me to make a much more accurate graph because I didn't have to worry about typing so many values in small increments of 50.I have a separate file that also records the log base 2 of the matrix size, and the elapsed time.*/
void print(vector< vector<double> >);
vector<double> gauss(vector< vector<double> >);

int main() {
  int n, i, j, r;
  double initial, end, elapsed;
  struct timeval tp;
  ofstream myfile;
  ofstream myotherfile;
  myfile.open("values.txt");
  myotherfile.open("log_values.txt");
 
 n = 2500;
for(double G=0;G<n+1;G=G+50) { 
  // Number of equations

  vector<double> line(G+1,0);
  vector< vector<double> > A(G,line);
  
  srand(time(0)); // seed the RNG with the current time
  
  for(i=0;i<G;i++) {
   for(j=0;j<G+1;j++) {
   r = rand() % 100;
   A[i][j]=r;
  }
 }

  // Calculate solution
  vector<double> x(G);
  // Call time function here..

  gettimeofday(&tp,NULL);
  initial= tp.tv_sec*1000+tp.tv_usec/1000;

  x = gauss(A);

  gettimeofday(&tp,NULL);
  end=tp.tv_sec*1000+tp.tv_usec/1000;

  elapsed= (end - initial);
  cout << "The total time elapsed for a system of: "<<G<<" equations is: "<<elapsed <<"."<< endl;
 
 
 myfile << G << "," << elapsed << endl;
 myotherfile << log2(G) << "," << log2(elapsed) << endl; 
}
myfile.close();
myotherfile.close();
}
void print(vector< vector<double> > A) {
  int n = A.size();
  for (int i=0; i<n; i++) {
    for (int j=0; j<n+1; j++) {
      cout << A[i][j] << "\t";
      if (j == n-1) {
	cout << "| ";
      }
    }
    cout << "\n";
  }
  cout << endl;
}

vector<double> gauss(vector< vector<double> > A) {
  int n = A.size();

  for (int i=0; i<n; i++) {
    // Search for maximum in this column
    double maxEl = abs(A[i][i]);
    int maxRow = i;
    for (int k=i+1; k<n; k++) {
      if (abs(A[k][i]) > maxEl) {
	maxEl = abs(A[k][i]);
	maxRow = k;
      }
    }

    // Swap maximum row with current row (column by column)
    for (int k=i; k<n+1;k++) {
      double tmp = A[maxRow][k];
      A[maxRow][k] = A[i][k];
      A[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      double c = -A[k][i]/A[i][i];
      for (int j=i; j<n+1; j++) {
	if (i==j) {
	  A[k][j] = 0;
	} else {
	  A[k][j] += c * A[i][j];
	}
      }
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  vector<double> x(n);
  for (int i=n-1; i>=0; i--) {
    x[i] = A[i][n]/A[i][i];
    for (int k=i-1;k>=0; k--) {
      A[k][n] -= A[k][i] * x[i];
    }
  }
  return x;
}


