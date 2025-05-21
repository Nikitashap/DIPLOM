#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include "1.cpp"

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

double d_func(double x){
    return gamma_function(1.2)*pow(x,1.8)/gamma_function(3);
}

double q_x_t(double alpha, double t, double x, double gamma, double beta){
    return gamma_function(4)*pow(x,2)*pow(t,2.2)/gamma_function(3.2) - pow(x,2)*pow(t,3);
}

double F(double x){
    return 0;
}

double F_L(double t){
    return 0;
}

double F_R(double t){
    return pow(t,3);
}



double exact_solution(double x, double t) {
    return t*t*x*x;  // Example exact solution
}



void solve_fractional_diff_eq2(int N, int M, int K,  double dx, double dt, double alpha, double beta, double gamma, double t_max ) {
    vector<vector<double>> c(N, vector<double>(M, 0.0));
    vector<vector<double>> A(N, vector<double>(M, 0.0));
    for (int n=1; n<N-1; n++){
        for(int m=0; m<M; m++){
            if (m<=n-1){
                A[n][m]=-1.0*g_func(alpha,n-m+1)*d_func(dx*n)*pow(dt,gamma)/pow(dx,alpha);
            }else if(m == n){
                A[n][m]=1.0-(g_func(alpha,1)*d_func(dx*n)*pow(dt,gamma)/pow(dx,alpha));
            }else if(m==n+1){
                A[n][m]=-1.0*g_func(alpha,0)*d_func(dx*n)*pow(dt,gamma)/pow(dx,alpha);
            }else if(m>=n+2){
                A[n][m]=0.0;
            }
        }
        A[0][0]=1.0;
        A[N-1][M-1]=1.0;
        for (int j=1; j<M; j++){
            A[0][j]=0.0;
        }
        for (int j=0; j<M-1; j++){
            A[N-1][j]=0.0;
        }
    }
    vector<vector<double>> A_inverse = inverseMatrix(A);


    

    // Set initial condition
    for (int i = 0; i < M; ++i) {
        c[0][i] = F(i*dx);  // Example initial condition
    }
   
    // Time-stepping loop
    for (int n = 0; n < N-1; n++) {
        c[n+1][0] = F_L((n+1)*dt); 
        c[n+1][M-1] = F_R((n+1)*dt);
        
        vector<double> c_n(M,0.0);
        
        vector<double> q_n(M,0.0);
        vector<double> half_res(M,0.0);


        for (int l = 1; l<=n+1; l++){
            for (int i=0; i<M; i++){
                half_res[i]= half_res[i]+c[n-l+1][i]*g_func(gamma,l);
            }

        }
        for (int i=0; i<M; i++){
                half_res[i]=-1.0*half_res[i];
            }
        for (int i = 0; i < M; i++) {
            q_n[i]=q_x_t(alpha,dt*(n+1),dx*i,gamma,beta)*pow(dt,gamma);
        }
        vector<double> res = multiplyMatrixVector(A_inverse,sumVectors(half_res,q_n));
        for (int i=1; i<M-1; i++){
            c[n+1][i]=res[i];
        }
    }

    // Write results to file
    ofstream outFile("time.txt");
    if (outFile.is_open()) {
        for (int n = 0; n < N; ++n) {
            for (int i = 0; i < M; ++i) {
                outFile << n*dt << "\t" << i*dx << "\t" << c[n][i] << "\t" << exact_solution(i*dx,n*dt) << endl;
            }
            outFile << endl;
        }
        outFile.close();
        cout << "Results written to file results.txt" << endl;
    } else {
        cout << "Failed to open file for writing!" << endl;
    }

    ofstream outFile2("time2.txt");
    if (outFile2.is_open()) {
            for (int i = 0; i < M; ++i) {
                outFile2 << 1 << "\t" << i*dx << "\t" << c[N-1][i] << "\t" << exact_solution(i*dx,1.0) << endl;
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
    double alpha = 1.8;  // Fractional order
    double gamma = 0.8;
    double beta = -1.0;
    double t_max = 1.0;



    solve_fractional_diff_eq2(N, M, K, dx, dt, alpha,  beta, gamma, t_max);

    

    return 0;
}
