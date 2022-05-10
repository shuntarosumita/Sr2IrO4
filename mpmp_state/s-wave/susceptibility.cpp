//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2017-11-21 16:32:24 shunta>

#include "hamiltonian.h"

//
// フェルミ分布関数
//
inline double Hamiltonian::fermi(double e, double beta) {
  return 1.0 / (exp(beta * e) + 1.0);
}

//
// 超伝導感受率の固有値を計算する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   T: 温度
//
void Hamiltonian::calculate_susceptibility(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu, double T, double U) {
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  Dvector k(3), q(3);                                 // 波数
  int qnx, nx, ny, nz, m[2], i[2];                    // ループ用変数
  long info;                                          // 対角化zheevの関数値
  double beta = 1.0 / T;                              // 逆温度

  // 配列の動的確保
  Complex *****Unitary = new Complex****[N[0] + 1];     // Normal part Hamiltonianを対角化するユニタリ行列
  double ****Eigen = new double***[N[0] + 1];           // Normal part Hamiltonianを対角化したときの固有値
  for(nx = 0; nx <= N[0]; nx++) {
    Unitary[nx] = new Complex***[N[1] + 1];
    Eigen[nx] = new double**[N[1] + 1];
    for(ny = 0; ny <= N[1]; ny++) {
      Unitary[nx][ny] = new Complex**[N[2] + 1];
      Eigen[nx][ny] = new double*[N[2] + 1];
      for(nz = 0; nz <= N[2]; nz++) {
        Unitary[nx][ny][nz] = new Complex*[DIM / 2];
        Eigen[nx][ny][nz] = new double[DIM / 2];
        for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
          Unitary[nx][ny][nz][i[0]] = new Complex[DIM / 2];
        }
      }
    }
  }

  // kx, ky, kz について，-pi から pi まで波数を変えて各点でGreen関数を計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[2]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // Normal part Hamiltonianを定義・対角化する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
        info = zmatrix.zheev_normal();

        for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
          for(i[1] = 0; i[1] < DIM / 2; i[1]++) {
            Unitary[nx][ny][nz][i[0]][i[1]] = zmatrix.hamil_normal(i[0], i[1]);
          }
          Eigen[nx][ny][nz][i[0]] = zmatrix.eigen_normal[i[0]];
        }

      }
    }
  }

  // ファイル出力関連
  stringstream filename;
  ofstream output;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(3) << "./data/sus_eigen_T" << T << setprecision(2) << "_U" << U
           << "_h" << h << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
  output.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  q[1] = 0.0;
  q[2] = 0.0;

  // qx について，波数を変えて各点で感受率を計算
  for(qnx = - N[0] / 8; qnx <= N[0] / 8; qnx++) {
    q[0] = 2.0 * M_PI * (double)qnx / (double)N[0];
    Cmatrix chi = zero_matrix<Complex>(8, 8);         // 超伝導感受率

    for(nx = 0; nx <= N[0]; nx++) {
      int nx2 = N[0] - nx + qnx;
      if(nx2 < 0) {
        nx2 = nx2 + N[0];
      }
      else if(nx2 > N[0]) {
        nx2 = nx2 - N[0];
      }

      for(ny = 0; ny <= N[1]; ny++) {
        int ny2 = N[1] - ny;

        for(nz = 0; nz <= N[2]; nz++) {
          int nz2 = N[2] - nz;

          for(m[0] = 0; m[0] < 8; m[0]++) {
            for(m[1] = 0; m[1] < 8; m[1]++) {
              for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
                for(i[1] = 0; i[1] < DIM / 2; i[1]++) {
                  Complex calc = 0.0;
                  calc += Unitary[nx][ny][nz][2 * m[0]][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1]][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0] + 1][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1] + 1][i[1]]); // (0 1 0 1)
                  calc += - Unitary[nx][ny][nz][2 * m[0]][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1] + 1][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0] + 1][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1]][i[1]]); // (0 1 1 0)
                  calc += - Unitary[nx][ny][nz][2 * m[0] + 1][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1]][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0]][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1] + 1][i[1]]); // (1 0 0 1)
                  calc += Unitary[nx][ny][nz][2 * m[0] + 1][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1] + 1][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0]][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1]][i[1]]); // (1 0 1 0)

                  // ゼロ割りを防止する
                  if(abs(Eigen[nx][ny][nz][i[0]] + Eigen[nx2][ny2][nz2][i[1]]) > 1.0e-8) {
                    calc *= (Hamiltonian::fermi(- Eigen[nx2][ny2][nz2][i[1]], beta) - Hamiltonian::fermi(Eigen[nx][ny][nz][i[0]], beta)) / (Eigen[nx][ny][nz][i[0]] + Eigen[nx2][ny2][nz2][i[1]]);
                  }
                  else {
                    calc *= beta * exp(beta * Eigen[nx][ny][nz][i[0]]) / pow(exp(beta * Eigen[nx][ny][nz][i[0]]) + 1.0, 2.0);
                  }

                  chi(m[0], m[1]) += calc;
                }
              }
            }
          }
        }

      }
    }

    chi *= U / V;

    // 固有値を求める
    Dvector eigen(chi.size1());
    matrix<Complex, column_major> chi_copy(chi.size1(), chi.size2());

    // コピーを作る
    for(size_t i = 0; i < chi.size1(); i++) {
      for(size_t j = i; j < chi.size2(); j++) {
        chi_copy(i, j) = chi(i, j);
      }
    }

    // lapackを用いて対角化
    info = lapack::heev('V', 'U', chi_copy, eigen, lapack::optimal_workspace());
    BOOST_UBLAS_CHECK(info == 0, internal_logic());

    // 超伝導感受率の最大固有値を出力
    output << fixed << setprecision(8) << q[0] << "  " << eigen[7] << endl;
    cout << q[0] << "  " << eigen << endl;

  }

  // メモリ解放
  for(nx = 0; nx <= N[0]; nx++) {
    for(ny = 0; ny <= N[1]; ny++) {
      for(nz = 0; nz <= N[2]; nz++) {
        for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
          delete[] Unitary[nx][ny][nz][i[0]];
        }
        delete[] Unitary[nx][ny][nz];
        delete[] Eigen[nx][ny][nz];
      }
      delete[] Unitary[nx][ny];
      delete[] Eigen[nx][ny];
    }
    delete[] Unitary[nx];
    delete[] Eigen[nx];
  }
  delete[] Unitary;
  delete[] Eigen;

}

//
// 超伝導感受率の固有値を計算する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   T: 温度
//
// void Hamiltonian::calculate_susceptibility(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu, double T, double U) {
//   Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
//   Dvector k(3), q(3);                                 // 波数
//   int qnx, nx, ny, nz, m[2], i[2];                    // ループ用変数
//   long info;                                          // 対角化zheevの関数値
//   double beta = 1.0 / T;                              // 逆温度

//   // 配列の動的確保
//   Complex *****Unitary = new Complex****[N[0] + 1];     // Normal part Hamiltonianを対角化するユニタリ行列
//   double ****Eigen = new double***[N[0] + 1];           // Normal part Hamiltonianを対角化したときの固有値
//   for(nx = 0; nx <= N[0]; nx++) {
//     Unitary[nx] = new Complex***[N[1] + 1];
//     Eigen[nx] = new double**[N[1] + 1];
//     for(ny = 0; ny <= N[1]; ny++) {
//       Unitary[nx][ny] = new Complex**[N[2] + 1];
//       Eigen[nx][ny] = new double*[N[2] + 1];
//       for(nz = 0; nz <= N[2]; nz++) {
//         Unitary[nx][ny][nz] = new Complex*[DIM / 2];
//         Eigen[nx][ny][nz] = new double[DIM / 2];
//         for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
//           Unitary[nx][ny][nz][i[0]] = new Complex[DIM / 2];
//         }
//       }
//     }
//   }

//   // kx, ky, kz について，-pi から pi まで波数を変えて各点でGreen関数を計算
//   for(nx = 0; nx <= N[0]; nx++) {
//     k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

//     for(ny = 0; ny <= N[1]; ny++) {
//       k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

//       for(nz = 0; nz <= N[2]; nz++) {
//         k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

//         // Normal part Hamiltonianを定義・対角化する
//         zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
//         info = zmatrix.zheev_normal();

//         for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
//           for(i[1] = 0; i[1] < DIM / 2; i[1]++) {
//             Unitary[nx][ny][nz][i[0]][i[1]] = zmatrix.hamil_normal(i[0], i[1]);
//           }
//           Eigen[nx][ny][nz][i[0]] = zmatrix.eigen_normal[i[0]];
//         }

//       }
//     }
//   }

//   // ファイル出力関連
//   stringstream filename;
//   ofstream output;

//   // データを出力するためのフォルダを作成
//   mkdir("./data", 0755);

//   // 出力ファイルを作成
//   filename << fixed << setprecision(3) << "./data/sus_eigen_T" << T << setprecision(2) << "_U" << U
//            << "_h" << h << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
//   output.open(filename.str().c_str());
//   filename.str("");                      // バッファのクリア
//   filename.clear(stringstream::goodbit); // フラグのクリア

//   q[1] = 0.0;
//   q[2] = 0.0;

//   // qx について，波数を変えて各点で感受率を計算
//   for(qnx = - N[0] / 8; qnx <= N[0] / 8; qnx++) {
//     q[0] = 2.0 * M_PI * (double)qnx / (double)N[0];
//     Cmatrix chi = zero_matrix<Complex>(8, 8);         // 超伝導感受率

//     for(nx = 0; nx <= N[0]; nx++) {
//       int nx2 = N[0] - nx + qnx;
//       if(nx2 < 0) {
//         nx2 = nx2 + N[0];
//       }
//       else if(nx2 > N[0]) {
//         nx2 = nx2 - N[0];
//       }

//       for(ny = 0; ny <= N[1]; ny++) {
//         int ny2 = N[1] - ny;

//         for(nz = 0; nz <= N[2]; nz++) {
//           int nz2 = N[2] - nz;

//           for(m[0] = 0; m[0] < 8; m[0]++) {
//             for(m[1] = 0; m[1] < 8; m[1]++) {
//               for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
//                 for(i[1] = 0; i[1] < DIM / 2; i[1]++) {
//                   Complex calc = 0.0;
//                   calc += Unitary[nx][ny][nz][2 * m[0]][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1]][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0] + 1][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1] + 1][i[1]]); // (0 1 0 1)
//                   calc += - Unitary[nx][ny][nz][2 * m[0]][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1] + 1][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0] + 1][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1]][i[1]]); // (0 1 1 0)
//                   calc += - Unitary[nx][ny][nz][2 * m[0] + 1][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1]][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0]][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1] + 1][i[1]]); // (1 0 0 1)
//                   calc += Unitary[nx][ny][nz][2 * m[0] + 1][i[0]] * conj(Unitary[nx][ny][nz][2 * m[1] + 1][i[0]]) * Unitary[nx2][ny2][nz2][2 * m[0]][i[1]] * conj(Unitary[nx2][ny2][nz2][2 * m[1]][i[1]]); // (1 0 1 0)

//                   // ゼロ割りを防止する
//                   if(abs(Eigen[nx][ny][nz][i[0]] + Eigen[nx2][ny2][nz2][i[1]]) > 1.0e-8) {
//                     calc *= (Hamiltonian::fermi(- Eigen[nx2][ny2][nz2][i[1]], beta) - Hamiltonian::fermi(Eigen[nx][ny][nz][i[0]], beta)) / (Eigen[nx][ny][nz][i[0]] + Eigen[nx2][ny2][nz2][i[1]]);
//                   }
//                   else {
//                     calc *= beta * exp(beta * Eigen[nx][ny][nz][i[0]]) / pow(exp(beta * Eigen[nx][ny][nz][i[0]]) + 1.0, 2.0);
//                   }

//                   chi(m[0], m[1]) += calc;
//                 }
//               }
//             }
//           }
//         }

//       }
//     }

//     chi *= U / V;

//     // 固有値を求める
//     Dvector eigen(chi.size1());
//     matrix<Complex, column_major> chi_copy(chi.size1(), chi.size2());

//     // コピーを作る
//     for(size_t i = 0; i < chi.size1(); i++) {
//       for(size_t j = i; j < chi.size2(); j++) {
//         chi_copy(i, j) = chi(i, j);
//       }
//     }

//     // lapackを用いて対角化
//     info = lapack::heev('V', 'U', chi_copy, eigen, lapack::optimal_workspace());
//     BOOST_UBLAS_CHECK(info == 0, internal_logic());

//     // 超伝導感受率の最大固有値を出力
//     output << fixed << setprecision(8) << q[0] << "  " << eigen[7] << endl;
//     cout << q[0] << "  " << eigen << endl;

//   }

//   // メモリ解放
//   for(nx = 0; nx <= N[0]; nx++) {
//     for(ny = 0; ny <= N[1]; ny++) {
//       for(nz = 0; nz <= N[2]; nz++) {
//         for(i[0] = 0; i[0] < DIM / 2; i[0]++) {
//           delete[] Unitary[nx][ny][nz][i[0]];
//         }
//         delete[] Unitary[nx][ny][nz];
//         delete[] Eigen[nx][ny][nz];
//       }
//       delete[] Unitary[nx][ny];
//       delete[] Eigen[nx][ny];
//     }
//     delete[] Unitary[nx];
//     delete[] Eigen[nx];
//   }
//   delete[] Unitary;
//   delete[] Eigen;

// }
