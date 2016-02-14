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

//g++ -std=c++11 -Wall -I/usr/local/include CASE1-dataset-generator.cc
//valgrind --leak-check=full --show-leak-kinds=all ./a.out

using namespace std;
#if 0
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
卒業研究の実験"CASE1"で用いたsynthetic dataのストリームデータのデータセット"CASE1-dataset.txt"を生成する関数。
Concept Driftは発生せず，時刻1~時刻1000までデータ源は定常
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
#endif

//std::default_random_engine generatorは"./func_noise_generator.cc"内で宣言されている
std::uniform_real_distribution<double> my_uniform_distribution(-10,10);

int main(){
	generator.seed(1410);
	double CASE1_input[1000];
	double CASE1_target[1000];

	for (int tau = 1; tau <= 1000; tau++){
		CASE1_input[tau-1] = my_uniform_distribution(generator);
		CASE1_target[tau-1] = sin(CASE1_input[tau-1]) + noise_generator();
	}

	for (int i = 0; i < 1000; i++){
		cout << setprecision(15);
		cout << CASE1_input[i] << " " << CASE1_target[i] << "\n";
	}
}