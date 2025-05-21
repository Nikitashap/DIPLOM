#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include "help.cpp"

using namespace std;


double gamma_function(double x) {
    return tgamma(x);
}

double g_func(double alpha, int k){
    return gamma_function(double(k)-alpha)/gamma_function(-alpha)/gamma_function(double(k)+1.0);
}

double d_plus_func(double x, double alpha){
    return gamma_function(3-alpha)/gamma_function(3)*pow(x,alpha);
}


double d_minus_func(double x, double alpha){
    return gamma_function(3.0-alpha)/gamma_function(3.0)*pow(1.0-x,alpha);
}

double q_x_t(double alpha, double beta, double t, double x, double t_max){
    return (exp(t_max*t)*(3*alpha-2*alpha*x-pow(alpha,2))*(1+beta))/4 + (exp(t_max*t)*(t_max-1)*pow(1-x,2));
}

double F(double x){
    return pow(1-x,2);
}

double F_L(double t, double t_max){
    return exp(t_max*t);
}

double F_R(double t, double t_max){
    return 0;
}



double exact_solution(double x, double t, double t_max) {
    return exp(t_max*t)*pow(1-x,2);  
}



void solve_fractional_diff_eq2(int N, int M, int K,  double dx, double dt, double alpha, function<double(double)> d, function<double(double, double)> q, double beta, double gamma, double t_max ) {
    vector<vector<double>> c(N, vector<double>(M, 0.0));
    vector<vector<double>> A(N, vector<double>(M, 0.0));
    for (int n=0; n<N; n++){
        for(int m=0; m<M; m++){
            if(m==n && m==0){
                A[n][m]=1.0;
            }else if(m==n && m==M-1){
                A[n][m]=1.0;
            }else if(n==0){
                A[n][m]= 0.0;
            }else if(n==N-1){
                A[n][m]=-0.0;
            }else if (m<=n-1){
                A[n][m]=-1.0*g_func(alpha,n-m+1)*d_plus_func(dx*n,alpha)*dt/pow(dx,alpha);;
            }else if(m == n){
                A[n][m]=1-g_func(alpha,1)*d_plus_func(dx*n,alpha)*dt/pow(dx,alpha);
            }else if(m==n+1){
                A[n][m]=-1.0*g_func(alpha,0)*d_plus_func(dx*n,alpha)*dt/pow(dx,alpha);
            }else if(m>=n+2){
                A[n][m]=0.0;
            }
        }
    }

    

    
    for (int i = 0; i < M; ++i) {
        c[0][i] = F(i*dx);  
    }
   
   
    for (int n = 0; n < N-1; n++) {
        c[n+1][0] = F_L((n+1)*dt, t_max); 
        c[n+1][M-1] = F_R((n+1)*dt, t_max);
        
        vector<double> c_n(M,0.0);
        vector<double> q_n(M,0.0);


        for (int i = 0; i < M; i++) {
            c_n[i]=c[n][i];
            q_n[i]=q_x_t(alpha,beta,dt*(n+1),dx*i,t_max)*dt;
        }
        vector<double> res = sumVectors(multiplyMatrixVector(inverseMatrix(A),c_n),multiplyMatrixVector(inverseMatrix(A),q_n));
        for (int i=1; i<M-1; i++){
            c[n+1][i]=res[i];
        }
    }

    
    ofstream outFile("results.txt");
    if (outFile.is_open()) {
        for (int n = 0; n < M; ++n) {
            for (int i = 0; i < N; ++i) {
                outFile << n*dt << "\t" << i*dx << "\t" << c[n][i] << "\t" << exact_solution(i*dx,n*dt,t_max) << endl;
            }
            outFile << endl;
        }
        outFile.close();
        cout << "Results written to file results.txt" << endl;
    } else {
        cout << "Failed to open file for writing!" << endl;
    }

    ofstream outFile2("res2.txt");
    if (outFile2.is_open()) {
            for (int i = 0; i < N; ++i) {
                outFile2 << 1 << "\t" << i*dx << "\t" << c[N-1][i] << "\t" << exact_solution(i*dx,1.0,t_max) << endl;
            }
            
        outFile2.close();
        cout << "Results written to file results.txt" << endl;
    } else {
        cout << "Failed to open file for writing!" << endl;
    }
}

int main() {
    int M = 101;  // Number of spatial points
    int N = 101;  // Number of time steps
    int K = 100;
    double dx = 0.01;  // Spatial step size
    double dt = 0.01;  // Time step size
    double alpha = 1.7;  // Fractional order
    double gamma = 1.0;
    double beta = 1.0;
    double t_max = 1.0;


    solve_fractional_diff_eq2(N, M, K, dx, dt, alpha, d, q,  beta, gamma, t_max);


    return 0;
}
