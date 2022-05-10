//-*- coding: utf-8 -*-
// Time-stamp: <2016-12-08 11:19:11 shunta>
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
  double h;                              // 磁気モーメントの大きさ
  Dvector theta(8);                      // 磁気モーメントの角度
  Dvector alpha(4);                      // 反対称スピン軌道相互作用の強さ
  double mu;                             // 化学ポテンシャル
  Dvector q = zero_vector<double>(3);    // 重心運動量

  // 固定パラメータの設定(-+-+状態)
  t[0] = 1.00; t[1] = 0.26; t[2] = 0.10;
  theta[0] = 348.0 * M_PI / 180.0;
  theta[1] = 192.0 * M_PI / 180.0;
  theta[2] = 168.0 * M_PI / 180.0;
  theta[3] = 12.0  * M_PI / 180.0;
  theta[4] = 348.0 * M_PI / 180.0;
  theta[5] = 192.0 * M_PI / 180.0;
  theta[6] = 168.0 * M_PI / 180.0;
  theta[7] = 12.0  * M_PI / 180.0;

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

  // electron_density(t, h, theta, alpha, mu);
  // output_eigen_normal_kx(t, h, theta, alpha, mu);
  output_eigen_normal2(t, h, theta, alpha, mu);
  // output_asym(t, h, theta, alpha, mu);

  return 0;
}
