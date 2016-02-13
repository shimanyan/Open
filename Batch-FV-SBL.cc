#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include "math.hpp"
#include <boost/random.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include "./func_co_variance_update.cc"
#include "./func_noise_generator.cc"
#include "./func_mean_update.cc"
#include "./func_d_cap_update.cc"
#include "./func_k_d_gaussian_kernel.cc"
#include "./func_k_d_online_prediction.cc"
#include "./func_delete_Check_old_l.cc"
#include "./func_k_d_design_matrix_update.cc"

using namespace std;

#if 0
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
文献
Dmitriy Shutin, Thomas Buchgraber, Sanjeev R. Kulkarni, and H. Vincent Poor.
Fast variational sparse bayesian learning with automatic relevance determination for superimposed signals. IEEE Transactions on Signal Processing, 59(12):6257– 6261, 2011.
で導入されているアルゴリズムを実装した。
FV-SBLのバッチ学習版となっている。
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
#endif

//g++ -Wall -std=c++11 -I/usr/local/include Batch-FV-SBL.cc
//valgrind --leak-check=full --show-leak-kinds=all ./a.out
int total_num = 0;
int dimension = 0;
FILE *gp;

typedef boost::numeric::ublas::matrix<double> type_Matrix;
typedef boost::numeric::ublas::vector<double> type_Vector;

//std::default_random_engine generatorは"./func_noise_generator.cc"内で宣言されている
std::uniform_real_distribution<double> my_uniform_distribution(-10,10);

int main(){
	generator.seed(1410);
    type_Matrix input_data_x(1,1);
    type_Matrix target_t(1,1);

    int data_num = 0;
    bool check = 0;
    string line;
    while(getline(cin,line)){
        if(line[0] == '#'){
            continue;
        }
        istringstream ss(line);

        if (check == 0){
            ss >> total_num;
            ss >> dimension;
            input_data_x.resize(total_num,dimension);
            target_t.resize(total_num,1);
            cout << "total_num:" << total_num << endl;
            cout << "dimension:" << dimension << endl;
            check = 1;
        }

        else{
          ss >> input_data_x(data_num,0);
          ss >> target_t(data_num,0);
          data_num++;
        }

        if(data_num == total_num){
            break;
        }
    }

	// INITIALIZE
	type_Matrix basic_func_check(total_num,1);
	for (int i = 0; i <= total_num-1; i++){
		basic_func_check(i,0) = 1;
	}
	int input_data_num = total_num;
	int basic_func_num = total_num;

    double noise_precision = 10.0;
    double c = pow(10,-6);
    double d = pow(10,-6);
    double c_cap,d_cap;
    
    type_Matrix design_matrix(total_num,total_num);
    design_matrix = k_d_design_matrix_update(input_data_num, basic_func_num, basic_func_check, input_data_x);

	type_Matrix alpha(total_num,1);
	for (int i = 0; i <= total_num-1; i++){
		alpha(i,0) = 10.0;
	}
	//alphaの初期値には0を与える

	type_Matrix co_variance(basic_func_num,basic_func_num);
	co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
	type_Matrix mean(basic_func_num,1);
	mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);


    for (int i = 0; i < 30; i++){

        double error = 0;
        cout << "残差を計算:" << endl;
        for (int i = 0; i <= 30; i++){
          double x = -10.0 + 2.0*(10.0/30.0)*i;
          type_Matrix x_vector(1,1);
          x_vector(0,0) = x;
          error += abs(sin(x) - k_d_online_predicton(x_vector, basic_func_num, basic_func_check, mean,input_data_x));
        }
        error = error/31.0;
        cout << i << "回目のイテレーションの開始時点での残差は:" << error << endl;

       int l = 0;
       while(1){

        	if (delete_Check_old_l(l, noise_precision ,basic_func_num, input_data_num, design_matrix, alpha,
                  input_data_x, input_data_x, target_t) > 0){
        		alpha(l,0) = delete_Check_old_l(l, noise_precision ,basic_func_num, input_data_num, design_matrix, alpha,
                  input_data_x, input_data_x, target_t);

        	    co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
        	    l++;
        	    if (l == (basic_func_num) ){
        		  break;
                }
        	}

        	
        	else{
        		int count = -1;
        		int basic_func_pos_2 = -1;
        		while(1){
        			basic_func_pos_2++;
        			if (basic_func_check(basic_func_pos_2,0) == 1){
        				count++;
        				if (count == l){
        					break;
        				}
        			}
        		}

        		basic_func_check(basic_func_pos_2,0) = 0;

        		type_Matrix alpha_without_l(basic_func_num-1,1);

        		for (int j = 0; j <= l-1; j++){
        			alpha_without_l(j,0) = alpha(j,0);
        		}
        		for (int j = l+1; j <= basic_func_num-1; j++){
        			alpha_without_l(j-1,0) = alpha(j,0);
        		}
        		alpha.resize(basic_func_num-1,1);
        		alpha = alpha_without_l;
        		basic_func_num--;

        		design_matrix = k_d_design_matrix_update(input_data_num, basic_func_num, basic_func_check, input_data_x);
        		co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
                mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);

                if (l == (basic_func_num) ){
        		  break;
                }
        	}
        }
   
    mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);
    c_cap = c + 0.5*input_data_num;
    d_cap = d_cap_update(d, input_data_num, basic_func_num, target_t, design_matrix, co_variance, mean);
    noise_precision = c_cap/d_cap;

   }

        
    gp = popen("gnuplot -persist","w");
	fprintf(gp, "set multiplot\n"); 
    fprintf(gp, "set xrange [-10:10]\n");
    fprintf(gp, "set yrange [-5.0:5.0]\n");
    fprintf(gp, "set xlabel \"time\"\n");	// ラベル表示
    fprintf(gp, "set ylabel \"error\"\n");
    fprintf(gp, "plot '-' with lines linestyle 6 linetype 6 linewidth 1 dt '_.-.' title \"y(x)\"\n");

    for (int i = 0; i <= 1000; i++){
        double x = -10.0 + 0.02*i;
        type_Matrix x_vector(1,1);
        x_vector(0,0) = x;
        fprintf(gp,"%f\t%f\n", x, k_d_online_predicton(x_vector, basic_func_num, basic_func_check, mean,input_data_x) );    // データの書き込み
    }
    fprintf(gp, "e\n");
    fflush(gp);
    pclose(gp);
}