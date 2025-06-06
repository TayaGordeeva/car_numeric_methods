#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

const int n = 200;
const int m = 200;
const int n2 = 2*n, m2=2*m;
const double a = 0.0, b = 2.0, c = 0.0, d = 1.0;
const double h = (b - a) / n, k = (d - c) / m;
const double h1 = (b - a) / n2, k1 = (d - c) / m2;

//const double tau = 1.59999616*pow(10, -6);
const double tau = (h*h*k*k)/(2*(h*h+k*k));
const double tau1 = (h1*h1*k1*k1)/(2*(h1*h1+k1*k1));


vector<double> v((n - 1) * (m - 1), 0.0);
vector<double> u((n - 1) * (m - 1), 0.0);


vector<double> b_r((n - 1) * (m - 1), 0.0);
vector<double> b_r_main1((n - 1) * (m - 1), 0.0);
vector<double> b_r_main2((n2 - 1) * (m2 - 1), 0.0);


//for test task

double f_test(double x, double y) {
    return -(pow(M_PI, 2) * (pow(x, 2) + pow(y, 2)) * exp(pow(sin(M_PI * x * y), 2)) * (pow(sin(2 * M_PI * x * y), 2) + 2 * cos(2 * M_PI * x * y)));
}

double u_test(double x, double y) {
    return exp(pow(sin(M_PI * x * y), 2));
}

void b_test(vector<double> &b_r_) {
    int ind = 0;
    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            ind = (j - 1) * (n - 1) + i - 1;
            if (i == 1) b_r_[ind] += u_test(0, k * j) / (h * h);
            if (i == n - 1) b_r_[ind] += u_test(2.0, k * j) / (h * h);
            if (j == 1) b_r_[ind] += u_test(i * h, 0) / (k * k);
            if (j == m - 1) b_r_[ind] += u_test(i * h, 1.0) / (k * k);

            b_r_[ind] += f_test(i * h, j * k);
        }
    }
}

void r_test(vector<double> &x, vector<double> &r) {
    double a_ = 2 * (1 / (h * h) + 1 / (k * k));
    vector<double> tmp((n - 1) * (m - 1), 0.0);
    for (int j = 1; j < m; j++) {
        for (int i = 1; i < n; i++) {
            int ind = (j - 1) * (n - 1) + i - 1;
            tmp[ind] = a_ * x[ind];

            if (1 < i) tmp[ind] -= x[ind - 1] / (h * h);
            if (i < n - 1) tmp[ind] -= x[ind + 1] / (h * h);
            if (j > 1) tmp[ind] -= x[ind - (n - 1)] / (k * k);
            if (j < m - 1) tmp[ind] -= x[ind + (n - 1)] / (k * k);
        }
    }

    for (int i = 0; i < (n - 1) * (m - 1); i++)
        r[i] = tmp[i] - b_r[i];
}

void x_s(vector<double> &x, vector<double> r) {
    for (int i = 0; i < (n - 1) * (m - 1); i++)
        x[i] = x[i] - tau * r[i];
}

//for main task
double f_main(double x, double y) {
    return -(abs(x*x - 2*y));
}

double mu1 (double y) {
    return pow(sin(M_PI*y), 2);
}

double mu2 (double y) {
    return pow(sin(2*M_PI*y), 2);
}

double mu3 (double x) {
    return pow(sin(M_PI*x), 2);
}

double mu4 (double x) {
    return pow(sin(2*M_PI*x), 2);
}

void b_(vector<double> &b_r_, int n_, int m_, double h_, double k_) {
    int ind = 0;
    for (int j = 1; j < m_; j++) {
        for (int i = 1; i < n_; i++) {
            ind = (j - 1) * (n_ - 1) + i - 1;
            if (i == 1) b_r_[ind] += mu1(k_ * j) / (h_ * h_);
            if (i == n_ - 1) b_r_[ind] += mu2(k_ * j) / (h_ * h_);
            if (j == 1) b_r_[ind] += mu3(i * h_) / (k_ * k_);
            if (j == m_ - 1) b_r_[ind] += mu4(i * h_) / (k_ * k_);

            b_r_[ind] += f_main(i * h_, j * k_);
        }
    }
}

void r_main_(vector<double> &x, vector<double> &r, vector<double> &b_r1, int n_, int m_, double h_, double k_) {
    double a_ = 2 * (1 / (h_ * h_) + 1 / (k_ * k_));
    vector<double> tmp((n_ - 1) * (m_ - 1), 0.0);
    for (int j = 1; j < m_; j++) {
        for (int i = 1; i < n_; i++) {
            int ind = (j - 1) * (n_ - 1) + i - 1;
            tmp[ind] = a_ * x[ind];

            if (1 < i) tmp[ind] -= x[ind - 1] / (h_ * h_);
            if (i < n_ - 1) tmp[ind] -= x[ind + 1] / (h_ * h_);
            if (j > 1) tmp[ind] -= x[ind - (n_ - 1)] / (k_ * k_);
            if (j < m_ - 1) tmp[ind] -= x[ind + (n_ - 1)] / (k_ * k_);
        }
    }

    for (int i = 0; i < (n_ - 1) * (m_ - 1); i++)
        r[i] = tmp[i] - b_r1[i];
}

void x_s_main(vector<double> &x, vector<double> r, int n_, int m_, double t) {
    for (int i = 0; i < (n_ - 1) * (m_ - 1); i++)
        x[i] = x[i] - t * r[i];
}

//write result
void writeParametersToCSV(int N, double eps) {
    ofstream file("parameters.txt");
    file << n << "\n";
    file << m << "\n";
    file << eps << "\n";
    file << N << "\n";
    file.close();
}

void writeTestResultsToCSV(int iterations, double final_norm, double max_error, double norm_0, double err, double x, double y) {
    ofstream file("test_results.txt");
    file << iterations << "\n"; //Iterations (l)
    file << norm_0 << "\n"; //start norm
    file << scientific << err << "\n";
    file << final_norm << "\n"; //Final norm
    file << max_error << "\n"; //Max error
    file << x <<"\n";
    file << y <<"\n";

    file.close();
}

void writeMainResultsToCSV(int coarse_iterations, double coarse_norm,  int fine_iterations, double fine_norm, double max_error, double xmax, double ymax, double n1, double n2, double r_1, double r_2) {
    ofstream file("main_results.txt");
    file << coarse_iterations << "\n"; //Coarse iterations (l)
    file << scientific << setprecision(10) << coarse_norm << "\n"; //Coarse final norm
    file << fine_iterations << "\n"; //Fine iterations (l1)
    file << scientific << setprecision(10) << fine_norm << "\n"; //Fine final norm
    file << scientific << setprecision(10) << max_error << "\n"; //Max error between grids
    file << fixed << setprecision(10) << xmax << "\n"; //Error location x
    file << fixed << setprecision(10) << ymax << "\n"; //Error location y
    file << fixed << setprecision(10) << n1 << "\n"; //start norm1
    file << fixed << setprecision(10) << n2 << "\n"; //start norm2
    file << fixed << setprecision(10) << r_1 << "\n"; //start r1
    file << fixed << setprecision(10) << r_2 << "\n"; //start r2
    file.close();
}


int main() {
    vector<double> x_i((n - 1) * (m - 1), 0.0), tmp_i, x_main1((n - 1) * (m - 1), 0.0), x_main2((n2 - 1) * (m2 - 1), 0.0);
    const int N = 2000000;
    const double eps = pow(10, -12), eps1 = pow(10, -13);
    double norm;
    int l = 0, l1 = 0;

    writeParametersToCSV(N, eps);

    cout<<tau1<<endl;

    //for test task
    //ofstream file_test_u("u_test.csv");
//
//
    //for (int i = 0; i < n; i++)
    //    file_test_u << u_test(a+i*h, 0.0) << ",";
    //file_test_u << u_test(2.0, 0.0) << "\n";
//
    //for (int i = 0; i < n; i++)
    //    file_test_u << u_test(a+i*h, 0.0) << ",";
    //file_test_u << u_test(2.0, 0.0) << "\n";
//
//
    //for (int j = 1; j < m; j++) {
    //    for (int i = 1; i < n; i++) {
    //        if (i==1) {
    //            file_test_u <<u_test(0.0, c + j * k)<<",";
    //        }
//
    //        int ind = (j - 1) * (n - 1) + (i - 1);
    //        double xi = a + i * h;
    //        double yj = c + j * k;
    //        u[ind] = u_test(xi, yj);
    //        file_test_u << u[ind] << ",";
//
    //        if (i == n-1) {
    //            file_test_u <<u_test(2.0, c + j * k);
    //        }
    //    }
    //    file_test_u << "\n";
    //}
//
    //for (int i = 0; i < n; i++)
    //    file_test_u << u_test(a+i*h, 1.0) << ",";
    //file_test_u << u_test(2.0, 1.0) << "";
//
//
    //file_test_u.close();
    //vector<double> r((n - 1) * (m - 1), 0.0);
    //b_test(b_r);
    //double norm_0; //start nevyazka for test
//
    //for (l = 0; l < N; l++) {
    //    norm = 0.0;
    //    tmp_i = x_i;
    //    r_test(x_i, r);
    //    x_s(x_i, r);
//
    //    for (int i = 0; i < (n - 1) * (m - 1); i++)
    //        norm = max(norm, abs(x_i[i] - tmp_i[i]));
//
    //    if (l==0) {
    //        norm_0=0.0;
    //        for (int i=0; i<(n - 1) * (m - 1); i++)
    //        norm_0 = max(norm_0, r[i]);
    //    }
//
    //    if (norm <= eps) {
    //        break;
    //    }
    //}
//
    double r_1=0.0;
    //for (int i=0; i<(n - 1) * (m - 1); i++)
    //    r_1 = max(r_1, abs(r[i]));
//
    ////v_test
//
    //ofstream file_test_v("v_test.csv");
//
    //for (int i = 0; i < n; i++)
    //    file_test_v << u_test(a+i*h, 0.0) << ",";
    //file_test_v << u_test(2.0, 0.0) << "\n";
//
    //for (int i = 0; i < n; i++)
    //    file_test_v << u_test(a+i*h, 0.0) << ",";
    //file_test_v << u_test(2.0, 0.0) << "\n";
//
    //for (int j = 1; j < m; ++j) {
    //    for (int i = 1; i < n; ++i) {
    //        if (i==1) {
    //            file_test_v <<u_test(0.0, c + j * k)<<",";
    //        }
    //        
    //        int ind = (j - 1) * (n - 1) + (i - 1);
    //        file_test_v << x_i[ind]<<",";
//
    //        if (i == n-1) {
    //            file_test_v <<u_test(2.0, c + j * k);
    //        }
    //    }
    //    file_test_v <<"\n";
    //}
//
    //for (int i = 0; i < n; i++)
    //    file_test_v << u_test(a+i*h, 1.0) << ",";
    //file_test_v << u_test(2.0, 1.0);
//
//
    //file_test_v.close();
//
    ////|u-v| file
//
    //ofstream file_test_uv("uv_test.csv");
//
    //for (int i = 0; i < n; i++)
    //    file_test_uv << "0.0,";
    //file_test_uv << "0.0\n";
//
    //for (int i = 0; i < n; i++)
    //    file_test_uv  << "0.0,";
    //file_test_uv << "0.0\n";
//
    //for (int j = 1; j < m; ++j) {
    //    for (int i = 1; i < n; ++i) {
    //        if (i==1) {
    //            file_test_uv <<"0.0,";
    //        }
    //        
    //        int ind = (j - 1) * (n - 1) + (i - 1);
    //        file_test_uv << abs(x_i[ind] - u[ind])<<",";
//
    //        if (i == n-1) {
    //            file_test_uv <<"0.0";
    //        }
    //    }
    //    file_test_uv <<"\n";
    //}
//
    //for (int i = 0; i < n; i++)
    //    file_test_uv << "0.0,";
    //file_test_uv << "0.0\n";
//
//
    //file_test_uv.close();
//
//
    double xmax = 0.0;
    double ymax = 0.0;
    double max_err = 0.0;
//
    //for (int j = 1; j < m; ++j) {
    //    for (int i = 1; i < n; ++i) {
    //        int ind_coarse = (j - 1) * (n - 1) + (i - 1);
    //        if (max_err < abs(x_i[ind_coarse] - u[ind_coarse])) {
    //            max_err = abs(x_i[ind_coarse] - u[ind_coarse]);
    //            ymax = j * k;
    //            xmax = i * h;
    //        }
    //    }
    //}
//
//
    //cout << max_err << " "<<l<<" "<<norm<<" " << 1/(4/(h*h)*pow(sin(M_PI/(2*n)), 2) + 4/(k*k)*pow(sin(M_PI/(2*m)), 2)) * r_1 << endl;
    //writeTestResultsToCSV(l, r_1, max_err, norm_0, norm, ymax, xmax);
//
//

    //for main task
    vector<double> r_main((n - 1) * (m - 1), 0.0);
    vector<double> r_main_0((n - 1) * (m - 1), 0.0);

    b_(b_r_main1, n, m, h, k);
    double norm_0_1, norm_0_2;

    for (l = 0; l < N; l++) {
        norm = 0.0;
        tmp_i = x_main1;
        r_main_(x_main1, r_main, b_r_main1, n, m, h, k);
        x_s_main(x_main1, r_main, n, m, tau);

        for (int i = 0; i < (n - 1) * (m - 1); i++)
            norm = max(norm, abs(x_main1[i] - tmp_i[i]));

        if (l==0) norm_0_1 = norm;


        if (norm <= eps) {
            break;
        }
    } // result: x_main1 
    double normm = norm;

    r_1=0.0;
    for (int i=0; i<(n - 1) * (m - 1); i++)
        r_1 = max(r_1, abs(r_main[i]));


    cout << "main task: (n, m) - " <<l<<" "<<norm<<endl;

    ofstream file_main_u("v1_main.csv");

    for (int i = 0; i < n; i++)
        file_main_u << mu3(a+i*h) << ",";
    file_main_u << mu3(2.0) << "\n";

    for (int i = 0; i < n; i++)
        file_main_u << mu3(a+i*h) << ",";
    file_main_u << mu3(2.0) << "\n";

   
    for (int j = 1; j < m; ++j) {
        for (int i = 1; i < n; ++i) {
            if (i==1) {
                file_main_u <<mu1(c + j * k)<<",";
            }
            
            int ind = (j - 1) * (n - 1) + (i - 1);
            file_main_u << x_main1[ind]<<",";

            if (i == n-1) {
                file_main_u <<mu2(c + j * k);
            }
        }
        file_main_u <<"\n";
    }

    for (int i = 0; i < n; i++)
        file_main_u << mu3(a+i*h) << ",";
    file_main_u << mu3(2.0);


    file_main_u.close();


    vector<double> r_main1((n2 - 1) * (m2 - 1), 0.0);
    b_(b_r_main2, n2, m2, h1, k1);

    for (l1 = 0; l1 < 2*N; l1++) {
        norm = 0.0;
        tmp_i = x_main2;
        r_main_(x_main2, r_main1, b_r_main2, n2, m2, h1, k1);
        x_s_main(x_main2, r_main1, n2, m2, tau1);

        for (int i = 0; i < (n2 - 1) * (m2 - 1); i++)
            norm = max(norm, abs(x_main2[i] - tmp_i[i]));

        if (l1==0) norm_0_2 = norm;

        if (norm <= eps1) {
            break;
        }
    } // result: x_main2

    double r_2=0.0;
    for (int i=0; i<(n2 - 1) * (m2 - 1); i++)
        r_2 = max(r_2, abs(r_main1[i]));

    cout << "           (2n, 2m) - " <<l1<<" "<<norm<<endl;

    ofstream file_main_v("v2_main.csv");

    for (int i = 0; i < n2; i++)
        file_main_v << mu3(a+i*h1) << ",";
    file_main_v << mu3(2.0) << "\n";

    for (int i = 0; i < n2; i++)
        file_main_v << mu3(a+i*h1) << ",";
    file_main_v << mu3(2.0) << "\n";

   
    for (int j = 1; j < m2; ++j) {
        for (int i = 1; i < n2; ++i) {
            if (i==1) {
                file_main_v <<mu1(c + j * k1)<<",";
            }
            
            int ind = (j - 1) * (n2 - 1) + (i - 1);
            file_main_v << x_main2[ind]<<",";

            if (i == n2-1) {
                file_main_v <<mu2(c + j * k1);
            }
        }
        file_main_v <<"\n";
    }

    for (int i = 0; i < n2; i++)
        file_main_v << mu3(a+i*h1) << ",";
    file_main_v << mu3(2.0);


    file_main_v.close();



    xmax = 0.0;
    ymax = 0.0;
    double max_err_main = 0.0;
    for (int j = 1; j < m; ++j) {
        for (int i = 1; i < n; ++i) {
            int ind_coarse = (j - 1) * (n - 1) + (i - 1);
            int ind_fine = (2*j - 1) * (n2 - 1) + (2*i - 1);
            if (max_err_main < abs(x_main1[ind_coarse] - x_main2[ind_fine])) {
                max_err_main = abs(x_main1[ind_coarse] - x_main2[ind_fine]);
                ymax = j * k;
                xmax = i * h;
            }
        }
    }

     cout << max_err_main << " " << xmax << " " << ymax;

    ofstream file_main_uv("uv_main.csv");

    for (int i = 0; i < n; i++)
        file_main_uv << "0.0,";
    file_main_uv << "0.0\n";

    for (int i = 0; i < n; i++)
        file_main_uv << "0.0,";
    file_main_uv << "0.0\n";

   
    for (int j = 1; j < m; ++j) {
        for (int i = 1; i < n; ++i) {
            if (i==1) {
                file_main_uv << "0.0,";
            }
            
            int ind = (j - 1) * (n - 1) + (i - 1);
            int ind_fine = (2*j - 1) * (n2 - 1) + (2*i - 1);

            file_main_uv << abs(x_main1[ind] - x_main2[ind_fine])<<",";

            if (i == n-1) {
                file_main_uv <<"0.0";
            }
        }
        file_main_uv <<"\n";
    }

    for (int i = 0; i < n; i++)
        file_main_uv << "0.0,";
    file_main_uv << "0.0";


    file_main_uv.close();

    writeMainResultsToCSV(l, normm, l1, norm, max_err_main, xmax, ymax, norm_0_1, norm_0_2, r_1, r_2);

    return 0;
}