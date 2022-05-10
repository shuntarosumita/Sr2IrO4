//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-11-30 11:23:10 shunta>

//
// Sr2IrO4のハミルトニアンライブラリヘッダファイル
//

// 二重インクルードを防止する。
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
//#include <vector>
#include <sys/stat.h>

// Boost
#include <boost/numeric/ublas/vector.hpp>     // ベクトル用のヘッダ
#include <boost/numeric/ublas/matrix.hpp>     // 普通の行列用のヘッダ
#include <boost/numeric/ublas/triangular.hpp> // 三角行列用のヘッダ．前進消去，後退代入に必要
#include <boost/numeric/ublas/lu.hpp>         // LU分解，前進消去，後退代入用のヘッダ
#include <boost/numeric/bindings/lapack/heev.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>         // ストリーム入出力用のヘッダ

// OPENMP
#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::complex<double> Complex;
typedef boost::numeric::ublas::matrix<Complex> Cmatrix;
typedef boost::numeric::ublas::vector<Complex> Cvector;
typedef boost::numeric::ublas::vector<double> Dvector;
using namespace std;
using namespace boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace Hamiltonian {

  // 各種固定パラメータ
  constexpr long N[3] = {256, 256, 8};               // サイト数
  constexpr double V = (double)(N[0] * N[1] * N[2]); // 体積
  constexpr long DIM = 32;                           // Hamiltonianの次元
  constexpr double delta = 1.0e-2;                   // 収束条件
  constexpr Complex imag_unit(0.0, 1.0);             // 虚数単位

  //
  // 複素数行列クラス
  //
  class Zmatrix {
  public:
    Cmatrix hamil;                              // Hamiltonian
    Cmatrix hamil_normal;                       // Hamiltonian (normal part)
    Dvector eigen;                              // 固有値を入れる配列
    Dvector eigen_normal;                       // 固有値を入れる配列 (normal part)

    //
    // コンストラクタ
    //
    Zmatrix() {
      hamil = zero_matrix<Complex>(DIM, DIM);
      eigen = Dvector(DIM, 0.0);
    }

    //
    // コンストラクタ
    //   k: 波数
    //   t: 飛び移り積分
    //   order: 秩序パラメータ
    //   h: 磁気モーメントの大きさ
    //   theta: 磁気モーメントの角度
    //   alpha: 反対称スピン軌道相互作用の強さ
    //   mu: 化学ポテンシャル
    //   q: Cooper対の重心運動量
    //
    Zmatrix(Dvector &k, Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
      init_zmatrix(k, t, order, h, theta, alpha, mu, q);
    }

    //
    // 初期化
    //
    void init_zmatrix(Dvector &k, Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
      unsigned int i, j;
      Dvector kplus = k;
      Dvector kminus = - k + q;
      Cmatrix diagonal_kplus = normal_part(kplus, t, alpha, mu) + moment_part(h, theta);
      Cmatrix diagonal_kminus = - trans(normal_part(kminus, t, alpha, mu) + moment_part(h, theta));
      Cmatrix offdiagonal = order_part(k, order);

      hamil = zero_matrix<Complex>(DIM, DIM);
      eigen = Dvector(DIM, 0.0);

      // Hamiltonianの左上
      for(i = 0; i < DIM / 2; i++) {
        for(j = 0; j < DIM / 2; j++) {
          hamil(i, j) = diagonal_kplus(i, j);
        }
      }      

      // Hamiltonianの右上
      for(i = 0; i < DIM / 2; i++) {
        for(j = 0; j < DIM / 2; j++) {
          hamil(i, j + DIM / 2) = offdiagonal(i, j);
        }
      }      

      // Hamiltonianの右下
      for(i = 0; i < DIM / 2; i++) {
        for(j = 0; j < DIM / 2; j++) {
          hamil(i + DIM / 2, j + DIM / 2) = diagonal_kminus(i, j);
        }
      }      

      for(i = 0; i < hamil.size1(); i++) {
        for(j = 0; j < i; j++) {
          hamil(i, j) = conj(hamil(j, i));
        }
      }

    }

    //
    // Normal part Hamiltonianの初期化
    //
    void init_zmatrix_normal(Dvector &k, Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
      hamil_normal = normal_part(k, t, alpha, mu) + moment_part(h, theta);
      eigen_normal = Dvector(DIM / 2, 0.0);
    }

    //
    // Hamiltonianを対角化する
    //
    long zheev(void) {
      Dvector e_copy(hamil.size1());
      matrix<Complex, column_major> h_copy(hamil.size1(), hamil.size2());
      long info;

      // Hamiltonianのコピーを作る
      for(size_t i = 0; i < hamil.size1(); i++) {
        for(size_t j = i; j < hamil.size2(); j++) {
          h_copy(i, j) = hamil(i, j);
        }
      }

      // lapackを用いて対角化
      info = lapack::heev('V', 'U', h_copy, e_copy, lapack::optimal_workspace());
      BOOST_UBLAS_CHECK(info == 0, internal_logic());

      // 固有値と固有ベクトルを代入
      hamil = h_copy;
      eigen.swap(e_copy);

      return info;
    }

    //
    // Normal part Hamiltonianを対角化する
    //
    long zheev_normal(void) {
      Dvector e_copy(hamil_normal.size1());
      matrix<Complex, column_major> h_copy(hamil_normal.size1(), hamil_normal.size2());
      long info;

      // Hamiltonianのコピーを作る
      for(size_t i = 0; i < hamil_normal.size1(); i++) {
        for(size_t j = i; j < hamil_normal.size2(); j++) {
          h_copy(i, j) = hamil_normal(i, j);
        }
      }

      // lapackを用いて対角化
      info = lapack::heev('V', 'U', h_copy, e_copy, lapack::optimal_workspace());
      BOOST_UBLAS_CHECK(info == 0, internal_logic());

      // 固有値と固有ベクトルを代入
      hamil_normal = h_copy;
      eigen_normal.swap(e_copy);

      return info;
    }

    //
    // Nearest neighbour hopping
    //
    double nn(Dvector &k, Dvector &t) {
      return - 4.0 * t[0] * cos(k[0] * 0.5) * cos(k[1] * 0.5);
    }

    //
    // Next nearest neighbour hopping
    //
    double nnn(Dvector &k, Dvector &t) {
      return - 2.0 * t[1] * (cos(k[0]) + cos(k[1]));
    }

    //
    // Third nearest neighbour hopping x
    //
    Complex tnnx(Dvector &k, Dvector &t) {
      return - t[2] * cos(k[0] * 0.5) * exp(- imag_unit * k[2] * 0.25);
    }

    //
    // Third nearest neighbour hopping y
    //
    Complex tnny(Dvector &k, Dvector &t) {
      return - t[2] * cos(k[1] * 0.5) * exp(- imag_unit * k[2] * 0.25);
    }

    //
    // パウリ行列
    //
    Cmatrix Sigma(long component) {
      Cmatrix mat = zero_matrix<Complex>(2, 2);
      if (component == 0) {
        // sigma_x
        mat(0, 1) = 1.0;
        mat(1, 0) = 1.0;
      }
      else if (component == 1) {
        // sigma_y
        mat(0, 1) = - imag_unit;
        mat(1, 0) = imag_unit;
      }
      else if (component == 2) {
        // sigma_z
        mat(0, 0) = 1.0;
        mat(1, 1) = - 1.0;
      }
      else {
        // それ以外は単位行列を返す
        mat = identity_matrix<Complex>(2);
      }

      return mat;
    }

    //
    // Order parameter part Hamiltonian
    //
    Cmatrix order_part(Dvector &k, Cvector &order) {
      Cmatrix mat = zero_matrix<Complex>(DIM / 2, DIM / 2);
      long i, spin;

      for (i = 0; i < 8; i++) {
        for (spin = 0; spin < 4; spin++) {
          mat(i * 2 + (spin / 2), i * 2 + (spin % 2)) = order[i] * imag_unit * Sigma(1)(spin / 2, spin % 2);
        }
      }

      return mat;
    }

  private:
    //
    // Normal part Hamiltonian
    //
    Cmatrix normal_part(Dvector &k, Dvector &t, Dvector &alpha, double mu) {
      long pm, layer, sl, spin;

      // 層内を表す部分
      // (副格子の足 -> スピンの足)
      Complex H_intra_layer[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_intra_layer[layer][sl][spin] =
              ((nnn(k, t) - mu) * Sigma(-1)(spin / 2, spin % 2) * Sigma(-1)(sl / 2, sl % 2)
              + nn(k, t) * Sigma(-1)(spin / 2, spin % 2) * Sigma(0)(sl / 2, sl % 2)) * Sigma(-1)(layer / 2, layer % 2);
          }
        }
      }

      // 層間を表す部分(その1)
      // (副格子の足 -> スピンの足)
      Complex H_inter_layer1[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_inter_layer1[layer][sl][spin] =
              (real(tnnx(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(-1)(sl / 2, sl % 2)
               + real(tnny(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(0)(sl / 2, sl % 2)) * Sigma(0)(layer / 2, layer % 2)
              - (imag(tnnx(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(-1)(sl / 2, sl % 2)
                 + imag(tnny(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(0)(sl / 2, sl % 2)) * Sigma(1)(layer / 2, layer % 2);
          }
        }
      }

      // 層間を表す部分(その2)
      // (副格子の足 -> スピンの足)
      Complex H_inter_layer2[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_inter_layer2[layer][sl][spin] =
              (real(tnny(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(-1)(sl / 2, sl % 2)
               + real(tnnx(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(0)(sl / 2, sl % 2)) * Sigma(0)(layer / 2, layer % 2)
              + (imag(tnny(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(-1)(sl / 2, sl % 2)
                 + imag(tnnx(k, t)) * Sigma(-1)(spin / 2, spin % 2) * Sigma(0)(sl / 2, sl % 2)) * Sigma(1)(layer / 2, layer % 2);
          }
        }
      }

      // 層内の反対称スピン軌道相互作用を表す部分(その1)
      // (副格子の足 -> スピンの足)
      Complex H_intra1_ASOC[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_intra1_ASOC[layer][sl][spin] =
              - alpha[0] * cos(k[0] * 0.5) * cos(k[1] * 0.5) * Sigma(2)(spin / 2, spin % 2) * Sigma(1)(sl / 2, sl % 2) * Sigma(-1)(layer / 2, layer % 2);
          }
        }
      }

      // 層内の反対称スピン軌道相互作用を表す部分(その2)
      // (副格子の足 -> スピンの足)
      Complex H_intra2_ASOC[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_intra2_ASOC[layer][sl][spin] =
              alpha[1] * ( sin(k[0]) * cos(k[1]) * Sigma(0)(spin / 2, spin % 2) - sin(k[1]) * cos(k[0]) * Sigma(1)(spin / 2, spin % 2) ) * Sigma(2)(sl / 2, sl % 2) * Sigma(2)(layer / 2, layer % 2);
          }
        }
      }

      // 層間の反対称スピン軌道相互作用を表す部分(その1・x方向)
      // (副格子の足 -> スピンの足)
      Complex H_inter1x_ASOC[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_inter1x_ASOC[layer][sl][spin] =
              alpha[2] * ( cos(k[2] * 0.25) * sin(k[0] * 0.5) * Sigma(0)(spin / 2, spin % 2) - 2.0 * sin(k[2] * 0.25) * cos(k[0] * 0.5) * Sigma(2)(spin / 2, spin % 2) ) * Sigma(1)(sl / 2, sl % 2) * Sigma(1)(layer / 2, layer % 2)
              + alpha[2] * ( sin(k[2] * 0.25) * sin(k[0] * 0.5) * Sigma(0)(spin / 2, spin % 2) + 2.0 * cos(k[2] * 0.25) * cos(k[0] * 0.5) * Sigma(2)(spin / 2, spin % 2) ) * Sigma(1)(sl / 2, sl % 2) * Sigma(0)(layer / 2, layer % 2);
          }
        }
      }

      // 層間の反対称スピン軌道相互作用を表す部分(その1・y方向)
      // (副格子の足 -> スピンの足)
      Complex H_inter1y_ASOC[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_inter1y_ASOC[layer][sl][spin] =
              - alpha[2] * ( cos(k[2] * 0.25) * sin(k[1] * 0.5) * Sigma(1)(spin / 2, spin % 2) - 2.0 * sin(k[2] * 0.25) * cos(k[1] * 0.5) * Sigma(2)(spin / 2, spin % 2) ) * Sigma(1)(sl / 2, sl % 2) * Sigma(1)(layer / 2, layer % 2)
              + alpha[2] * ( sin(k[2] * 0.25) * sin(k[1] * 0.5) * Sigma(1)(spin / 2, spin % 2) + 2.0 * cos(k[2] * 0.25) * cos(k[1] * 0.5) * Sigma(2)(spin / 2, spin % 2) ) * Sigma(1)(sl / 2, sl % 2) * Sigma(0)(layer / 2, layer % 2);
          }
        }
      }

      // 層間の反対称スピン軌道相互作用を表す部分(その2)
      // (副格子の足 -> スピンの足)
      Complex H_inter2_ASOC[4][4][4];
      for (layer = 0; layer < 4; layer++) {
        for (sl = 0; sl < 4; sl++) {
          for (spin = 0; spin < 4; spin++) {
            H_inter2_ASOC[layer][sl][spin] =
              alpha[3] * cos(k[2] * 0.5) * ( sin(k[0] * 0.5) * cos(k[1] * 0.5) * Sigma(0)(spin / 2, spin % 2) - sin(k[1] * 0.5) * cos(k[0] * 0.5) * Sigma(1)(spin / 2, spin % 2) ) * Sigma(2)(sl / 2, sl % 2) * Sigma(2)(layer / 2, layer % 2);
          }
        }
      }

      // 行列を一度配列を用いて作る
      // (+-の足 -> 層の足 -> 副格子の足 -> スピンの足)
      Complex mat_copy[4][4][4][4];
      for (pm = 0; pm < 4; pm++) {

        if (pm == 0 or pm == 3) {
          for (layer = 0; layer < 4; layer++) {
            for (sl = 0; sl < 4; sl++) {
              for (spin = 0; spin < 4; spin++) {
                mat_copy[pm][layer][sl][spin] =
                  H_intra_layer[layer][sl][spin]
                  + H_inter_layer1[layer][sl][spin]
                  + H_intra1_ASOC[layer][sl][spin]
                  + H_intra2_ASOC[layer][sl][spin]
                  + H_inter1y_ASOC[layer][sl][spin];
              }
            }
          }
        }

        else if (pm == 1) {
          for (layer = 0; layer < 4; layer++) {
            for (sl = 0; sl < 4; sl++) {
              for (spin = 0; spin < 4; spin++) {
                mat_copy[pm][layer][sl][spin] =
                  H_inter_layer2[layer][sl][spin]
                  + H_inter1x_ASOC[layer][sl][spin]
                  + H_inter2_ASOC[layer][sl][spin];
              }
            }
          }
        }

        else {
          for (layer = 0; layer < 4; layer++) {
            for (sl = 0; sl < 4; sl++) {
              for (spin = 0; spin < 4; spin++) {
                mat_copy[pm][layer][sl][spin] = 0.0;
              }
            }
          }
        }

      }

      // 配列をCmatrixへ移す
      Cmatrix mat(DIM / 2, DIM / 2);
      for (pm = 0; pm < 4; pm++) {
        for (layer = 0; layer < 4; layer++) {
          for (sl = 0; sl < 4; sl++) {
            for (spin = 0; spin < 4; spin++) {
              mat((pm / 2) * 8 + (layer / 2) * 4 + (sl / 2) * 2 + (spin / 2), (pm % 2) * 8 + (layer % 2) * 4 + (sl % 2) * 2 + (spin % 2)) = mat_copy[pm][layer][sl][spin];
            }
          }
        }
      }

      for(unsigned int i = 0; i < mat.size1(); i++) {
        for(unsigned int j = 0; j < i; j++) {
          mat(i, j) = conj(mat(j, i));
        }
      }

      return mat;

    }

    //
    // Moment part Hamiltonian
    //
    Cmatrix moment_part(double h, Dvector &theta) {
      Cmatrix mat = zero_matrix<Complex>(DIM / 2, DIM / 2);
      unsigned int i, spin;

      for (i = 0; i < 8; i++) {
        for (spin = 0; spin < 4; spin++) {
          mat(i * 2 + (spin / 2), i * 2 + (spin % 2)) = - h * ( cos(theta[i]) * Sigma(0)(spin / 2, spin % 2) + sin(theta[i]) * Sigma(1)(spin / 2, spin % 2) );
        }
      }

      return mat;
    }

  };

  // フェルミ分布関数
  inline double fermi(double e, double beta);

  // 秩序パラメータを計算する
  Cvector calculate_order_parameter(Dvector &t, Cvector &order_in, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q, double beta);

  // 自由エネルギーを計算する
  double calculate_free_energy(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q, double beta);

  // 超伝導感受率の固有値を計算する
  void calculate_susceptibility(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu, double T, double U);

  // 固有値をファイルに出力する
  void output_eigen(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q);

  void test(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q);

  // ギャップの角度依存性を出力する
  void output_gap_angle(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q);

  // Normal part Hamiltonianの固有値をファイルに出力する
  void output_eigen_normal(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu);
}

#endif // HAMILTONIAN_H
