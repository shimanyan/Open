#include <cassert>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <algorithm>
typedef boost::numeric::ublas::matrix<double> type_Matrix;

using namespace std;
#if 0
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
文献
M. Drozdzal, J. Vitria, S. Segui, C. Malagelada, F. Azpiroz, and P. Radeva. Intestinal event segmentation for endoluminal video analysis. 
In Image Processing (ICIP),2014 IEEE International Conference on, pages 3592–3596, Oct 2014.
で導入されているストリームデータに対するSegmentation Algorithmを用いて
sliding-windowブロックの長さの最適化を行なう関数
＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿
#endif

int SEGMENTATION(int k, int tau, type_Matrix &gap){
  int width = k;
  double max = -pow(10,6);
  double min = pow(10,6);
  //この範囲に収まらないデータは必ず正規化する必要がある。
  type_Matrix data_sequence(gap.size1(),1);
  data_sequence = gap;
  for (int i = tau-width+1; i <= tau; i++){
    if (data_sequence(i,0) > max){
      max = data_sequence(i,0);
    }
    if (data_sequence(i,0) < min){
      min = data_sequence(i,0);
    }
  }

  for (int i = tau-width+1; i <= tau; i++){
    data_sequence(i,0) = abs(data_sequence(i,0) - min) / abs(max-min) ;
    #if 0
    SEGMENTATION ALGORITHMを適用するデータ列に含まれるすべての要素はその値を[0,1]にリスケーリングする必要があるため、
    このようにリスケーリングする

    きちんとリスケーリングされているかのテスト:
    assert((data_sequence(i,0) >= 0)&&(data_sequence(i,0) <= 1.0));
    =>通過！
    #endif
  }

  double delta = 0.1;

  for (int i = tau-width+1; i <= tau-1; i++){
    //n_0 >= 1 かつ n_1 >= 1 となる範囲で動かす
    int n_0 = i - (tau-width+1) + 1;
    //1番古い要素からi番目までのデータ数
    int n_1 = tau - (i+1) + 1;
    //cout << n_0 << "+" << n_1 << "=" << width << endl;
    //i+1番目から一番新しい要素までのデータ数
    #if 0
    assert((n_0 >= 1)&&(n_1 >= 1));
    =>通過！
    #endif

    double u_W0 = 0;
    double u_W1 = 0;



    double H = 1.0 / ((1.0/n_0) + (1.0/n_1));
    double delta_2 = delta/width;
    double eps_cut = sqrt((log(4.0/delta_2)) / (2.0*H));

    for (int j = tau-width+1; j <= i; j++){
      u_W0 += data_sequence(j,0);
    }
    u_W0 = u_W0/n_0;

    for (int j = i+1; j <= tau; j++){
      u_W1 += data_sequence(j,0);
    }
    u_W1 = u_W1/n_1;

    if (abs(u_W0 - u_W1) >= eps_cut){
      width = n_1;
      break;
    }
  }

  return std::max(width,30);
}