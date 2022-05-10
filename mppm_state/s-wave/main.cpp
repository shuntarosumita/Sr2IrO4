//-*- coding: utf-8 -*-
// Time-stamp: <2017-10-18 18:05:01 shunta>
// ----------------------------------------
// Sr2IrO4の超伝導計算
// ----------------------------------------
#include "hamiltonian.h"
typedef std::complex<double> Complex;
using namespace std;
using namespace Hamiltonian;

// メインプログラム
int main(int argc, char* argv[]) {
  // 物理量
  Dvector t(3);                          // 飛び移り積分
  //Cvector order(8);                      // 秩序パラメータ
  boost::numeric::ublas::vector<Cvector> order(3, Cvector(8));
  double h;                              // 磁気モーメントの大きさ
  Dvector theta(8);                      // 磁気モーメントの角度
  Dvector alpha(4);                      // 反対称スピン軌道相互作用の強さ
  double mu;                             // 化学ポテンシャル
  Dvector q = zero_vector<double>(3);    // 重心運動量

  // 固定パラメータの設定(-++-状態)
  t[0] = 1.00; t[1] = 0.26; t[2] = 0.10;
  theta[0] = 348.0 * M_PI / 180.0;
  theta[1] = 192.0 * M_PI / 180.0;
  theta[2] = 168.0 * M_PI / 180.0;
  theta[3] = 12.0  * M_PI / 180.0;
  theta[4] = 168.0 * M_PI / 180.0;
  theta[5] = 12.0  * M_PI / 180.0;
  theta[6] = 348.0 * M_PI / 180.0;
  theta[7] = 192.0 * M_PI / 180.0;

  // 秩序パラメータの初期化
  Complex Delta = 0.02;
  Cvector order_init = Cvector(8, Delta);

  // 可変パラメータを設定
  if(argc > 5) {
    h = atof(argv[1]);
    if(h < 0.0 || h > 1.0) {
      cerr << "h range: [0.0, 1.0]" << endl;
      return 0;
    }

    alpha[0] = atof(argv[2]);
    if(alpha[0] < 0.0 || alpha[0] > 1.0) {
      cerr << "alpha1 range: [0.0, 1.0]" << endl;
      return 0;
    }

    alpha[1] = atof(argv[3]);
    if(alpha[1] < 0.0 || alpha[1] > 1.0) {
      cerr << "alpha2 range: [0.0, 1.0]" << endl;
      return 0;
    }

    alpha[2] = atof(argv[4]);
    alpha[3] = 0.00;
    if(alpha[2] < 0.0 || alpha[2] > 1.0) {
      cerr << "alpha2 range: [0.0, 1.0]" << endl;
      return 0;
    }

    mu = atof(argv[5]);
    if(mu <= - 3 || mu > 3) {
      cerr << "mu range: (-3, 3]" << endl;
      return 0;
    }
  }
  else {
    cerr << "Usage: " << argv[0] << " h alpha1 alpha2 alpha3 mu" << endl;
    return 0;
  }

  // // 秩序パラメータを計算
  // double T = 0.001;                      // 温度
  // double beta = 1.0 / T;                 // 逆温度
  // long loopcount, i;
  // loopcount = 0;
  // order[0] = order_init;
  // order[1] = calculate_order_parameter(t, order[0], h, theta, alpha, mu, q, beta);
  // order[2] = calculate_order_parameter(t, order[1], h, theta, alpha, mu, q, beta);
  // while(abs((order[0][0] - order[1][0]) / order[0][0]) > delta) {
  //   cout << fixed << setprecision(8) << loopcount << "  " << order[0][0] << endl;

  //   loopcount += 1;
  //   for(i = 0; i < 8; i++) {
  //     order[0][i] = abs(order[0][i] - pow(order[1][i] - order[0][i], 2.0) / (order[0][i] - 2.0 * order[1][i] + order[2][i]));
  //   }

  //   order[1] = calculate_order_parameter(t, order[0], h, theta, alpha, mu, q, beta);
  //   order[2] = calculate_order_parameter(t, order[1], h, theta, alpha, mu, q, beta);

  //   if(loopcount > 200) {
  //     return 0;
  //   }
  // }

  // output_eigen_normal(t, h, theta, alpha, mu);

  order[0] = order_init;
  // output_eigen_kxky(t, order[0], h, theta, alpha, mu, q);
  // output_eigen_kykz(t, order[0], h, theta, alpha, mu, q);
  // output_eigen_kzkx(t, order[0], h, theta, alpha, mu, q);
  // output_eigen2(t, order[0], h, theta, alpha, mu, q);
  output_bandgap(t, order[0], h, theta, alpha, mu);
  // output_gap_angle(t, order[0], h, theta, alpha, mu, q);
  // test_kx(t, order[0], h, theta, alpha, mu, q);
  // test_ky(t, order[0], h, theta, alpha, mu, q);
  // test_kz(t, order[0], h, theta, alpha, mu, q);

  return 0;
}
