#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <random>
#include "math.hpp"
#include <iomanip>
#include <boost/random.hpp>
#include "./func_noise_generator.cc"

//g++ -std=c++11 -Wall -I/usr/local/include CASE2-dataset-generator.cc
//valgrind --leak-check=full --show-leak-kinds=all ./a.out

using namespace std;
#if 0
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
卒業研究の実験"CASE2"で用いたsynthetic dataのストリームデータのデータセット"CASE2-dataset.txt"を生成する関数。
Sudden Concept Driftの発生を想定。
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
#endif


//std::default_random_engine generatorは"./func_noise_generator.cc"内で宣言されている
std::uniform_real_distribution<double> my_uniform_distribution(-10,10);

int main(){
	generator.seed(1410);
	double CASE2_input[1000];
	double CASE2_target[1000];
	double parameter[1000];
	for (int i = 1; i <= 500; i++){
		parameter[i-1] = 0.0;
	}

	for (int i = 501; i <= 1000; i++){
		parameter[i-1] = 2.0;
	}


	for (int tau = 1; tau <= 1000; tau++){
		CASE2_input[tau-1] = my_uniform_distribution(generator);
		CASE2_target[tau-1] = sin(CASE2_input[tau-1]) +parameter[tau-1] + noise_generator();
	}

	for (int i = 0; i < 1000; i++){
		cout << setprecision(15);
		cout << CASE2_input[i] << " " << CASE2_target[i] << "\n";
	}
}