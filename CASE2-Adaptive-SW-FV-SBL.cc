#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
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
#include "./func_k_d_add_Check_new_l.cc"
#include "./func_k_d_online_prediction.cc"
#include "./func_k_d_online_prediction_SW.cc"
#include "./func_delete_Check_old_l.cc"
#include "./func_k_d_design_matrix_update.cc"
#include "./func_k_d_design_matrix_update_SW.cc"
#include "./func_SEGMENTATION.cc"

using namespace std;
#if 0
___________________________________＿＿＿＿＿
卒論のCASE2の実験のためのAdaptive-SW-FV-SBLのプログラム。
CASE2はSudden Concept Driftの発生を想定
___________________________________＿＿＿＿＿
#endif

//g++ -Wall -std=c++11 -I/usr/local/include CASE2-Adaptive-SW-FV-SBL.cc
//valgrind --leak-check=full --show-leak-kinds=all ./a.out < CASE2-dataset.txt
int total_num = 0;
int dimension = 0;
FILE *gp;

typedef boost::numeric::ublas::matrix<double> type_Matrix;
typedef boost::numeric::ublas::vector<double> type_Vector;

//std::default_random_engine generatorは"./func_noise_generator.cc"内で宣言されている
std::uniform_real_distribution<double> my_uniform_distribution(-10,10);

int main(){
	generator.seed(1510);
    type_Matrix input_data_container(1,1);
    type_Matrix target_t_container(1,1);

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
            input_data_container.resize(total_num,dimension);
            target_t_container.resize(total_num,1);
            check = 1;
        }

        else{
          ss >> input_data_container(data_num,0);
          ss >> target_t_container(data_num,0);
          data_num++;
        }

        if(data_num == total_num){
            break;
        }
    }

	// INITIALIZE
    double offset = 0.0;

    type_Matrix gap(total_num,1);
    double prediction_error[total_num];
    for (int i = 0; i < total_num; i++){
        prediction_error[i] = 0.0;
    }

    int tau = 0; //時間
    int basic_func_num = 1;
    int input_data_num = 1;
    int k_max = 300; //sliding-windowブロックの長さの最大値
    int k_min = 30;
	type_Matrix input_data_x(1,dimension);
    type_Matrix target_t(1,1);

    for (int i = 0; i <= dimension-1; i++){
        input_data_x(0,i) = input_data_container(0,i);
    }
    target_t(0,0) = target_t_container(0,0);

	type_Matrix basic_func_check(total_num,1);
	for (int i = 0; i <= total_num-1; i++){
		basic_func_check(i,0) = 0;
	}
    basic_func_check(0,0) = 1;

    double noise_precision = 10.0;
    double c = 0;
    double d = 0;
    double c_cap,d_cap;
    
    type_Matrix design_matrix(1,1);
    design_matrix = k_d_design_matrix_update_SW(input_data_num, basic_func_num, basic_func_check, input_data_x, input_data_container);

	type_Matrix alpha(1,1);
	alpha(0,0) = 0.0;
	//alphaの初期値には0を与える

	type_Matrix co_variance(basic_func_num,basic_func_num);
	co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
	type_Matrix mean(basic_func_num,1);
	mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);

    type_Matrix vector_for_prediction(1,1);
    vector_for_prediction(0,0) = input_data_x(0,0);
    gap(0,0) = abs(target_t(0,0) - k_d_online_predicton_SW(vector_for_prediction, basic_func_num, basic_func_check, mean,input_data_x,input_data_container));

    cout << "時刻:" << tau << "\n";
    type_Matrix x_vector(1,1);
    for (int i = 0; i < 100; i++){
        x_vector(0,0) = my_uniform_distribution(generator);
        prediction_error[0] += 100.0*abs(sin(x_vector(0,0)) + offset - k_d_online_predicton_SW(x_vector, basic_func_num, basic_func_check, mean,input_data_x,input_data_container));
    }
    prediction_error[0] = 0.01*(0.01*prediction_error[0]);

    for (tau = 1; tau <= total_num-1; tau++){
        if (tau == 500){
            offset = 2.0;
        }

        c_cap = c + 0.5*input_data_num;
        d_cap = d_cap_update(d, input_data_num, basic_func_num, target_t, design_matrix, co_variance, mean);
        noise_precision = c_cap/d_cap;

        if (input_data_num < k_max){
            input_data_num++;
        }
        input_data_x.resize(input_data_num,dimension);
        target_t.resize(input_data_num,1);


        type_Matrix new_input_data_x(input_data_num,dimension);
        type_Matrix new_target_t(input_data_num,1);
        for (int num = 0; num < input_data_num; num++){
          for (int i = 0; i <= dimension-1; i++){
            new_input_data_x(num,i) = input_data_container(tau-input_data_num+1 + num,i);
          }
            new_target_t(num,0) = target_t_container(tau-input_data_num+1 + num,0);
        }

        input_data_x = new_input_data_x;
        target_t = new_target_t;

        vector_for_prediction(0,0) = input_data_x(input_data_num-1,0);
        gap(tau,0) = abs(target_t(input_data_num-1,0) - k_d_online_predicton_SW(vector_for_prediction, basic_func_num, basic_func_check, mean,input_data_x,input_data_container));
        if (input_data_num > k_min){
          input_data_num = SEGMENTATION(input_data_num, tau, gap);
        }
        cout << "input_data_num:" << input_data_num << endl;

        input_data_x.resize(input_data_num,dimension);
        target_t.resize(input_data_num,1);


        new_input_data_x.resize(input_data_num,dimension);
        new_target_t.resize(input_data_num,1);
        for (int num = 0; num < input_data_num; num++){
          for (int i = 0; i <= dimension-1; i++){
            new_input_data_x(num,i) = input_data_container(tau-input_data_num+1 + num,i);
          }
            new_target_t(num,0) = target_t_container(tau-input_data_num+1 + num,0);
        }

        input_data_x = new_input_data_x;
        target_t = new_target_t;

        design_matrix = k_d_design_matrix_update_SW(input_data_num, basic_func_num, basic_func_check, input_data_x, input_data_container);
        mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);
        co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);

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
                design_matrix = k_d_design_matrix_update_SW(input_data_num, basic_func_num, basic_func_check, input_data_x, input_data_container);
        		co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
                mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);

                if (l == (basic_func_num) ){
        		  break;
                }
        	}
        }

       if (k_d_add_Check_new_l(tau, noise_precision ,basic_func_num, input_data_num, design_matrix, alpha,
         input_data_container, input_data_x, target_t) > 0){
         double new_alpha = k_d_add_Check_new_l(tau, noise_precision ,basic_func_num, input_data_num, design_matrix, alpha,
         input_data_container, input_data_x, target_t);
         basic_func_num++;
         basic_func_check(tau,0) = 1;
         design_matrix = k_d_design_matrix_update_SW(input_data_num, basic_func_num, basic_func_check, input_data_x, input_data_container);
         alpha.resize(basic_func_num,1);
         alpha(basic_func_num-1,0) = new_alpha;
         co_variance = co_variance_update(noise_precision, basic_func_num, design_matrix, alpha);
       }
   
    mean = mean_update(noise_precision, basic_func_num, design_matrix, co_variance, target_t);

    cout << "時刻:" << tau << "\n";
    for (int i = 0; i < 100; i++){
        x_vector(0,0) = my_uniform_distribution(generator);
        prediction_error[tau] += 100.0*abs(sin(x_vector(0,0)) + offset - k_d_online_predicton_SW(x_vector, basic_func_num, basic_func_check, mean,input_data_x,input_data_container));
    }
    prediction_error[tau] = 0.01*(0.01*prediction_error[tau]);
}

for (int i = 0; i < total_num; i++){
        cout << setprecision(15);
        cout << prediction_error[i] << endl;
    }

    gp = popen("gnuplot -persist","w");
    fprintf(gp, "set multiplot\n"); 
    fprintf(gp, "set xrange [1:1000]\n");
    fprintf(gp, "set yrange [0:2.5]\n");
    fprintf(gp, "set xlabel \"Time\"\n");   // ラベル表示
    fprintf(gp, "set ylabel \"predicrion error\"\n");
    fprintf(gp, "plot '-' with lines linestyle 6 linetype 6 linewidth 1 dt '_.-.' title \"Adaptive-SW-FV-SBL\"\n");

    for (int i = 1; i <= 1000; i++){
        fprintf(gp,"%f\t%f\n", (double)i, prediction_error[i-1]);    // データの書き込み
    }
    fprintf(gp, "e\n");
    fflush(gp);
    pclose(gp);
}