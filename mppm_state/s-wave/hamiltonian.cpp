//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2017-02-09 15:25:02 shunta>

//
// Sr2IrO4のハミルトニアンライブラリ
//

#include "hamiltonian.h"

//
// フェルミ分布関数
//
inline double Hamiltonian::fermi(double e, double beta) {
  return 1.0 / (exp(beta * e) + 1.0);
}

//
// 秩序パラメータを計算する
//   t: 飛び移り積分
//   order_in: 秩序パラメータ(入力)
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//   beta: 逆温度
//
Cvector Hamiltonian::calculate_order_parameter(Dvector &t, Cvector &order_in, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q, double beta) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  Cvector order_out(8);                               // 秩序パラメータ(出力)
  int nx, ny, nz, i, j;                               // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // kx, ky, kz について，-pi から pi まで波数を変えて各点で固有ベクトルを計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[2]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order_in, h, theta, alpha, mu, q);
    
        // 対角化する
        info = zmatrix.zheev();

        // 固有ベクトルをもとに秩序パラメータを計算する
        for(i = 0; i < DIM; i++) {
          for(j = 0; j < 8; j++) {
            order_out[j] += - U / V * conj(zmatrix.hamil(i, 2 * j + 17)) * zmatrix.hamil(i, 2 * j) * fermi(zmatrix.eigen[i], beta);
          }
        }
      }
    }
  }

  // 計算した秩序パラメータを返す
  return order_out;
}

//
// 固有値をファイルに出力する(E = 0に近いところのみ)
//   t: 飛び移り積分
//   order: 秩序パラメータ
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//
void Hamiltonian::output_eigen_kxky(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3), wave_num(3);                          // 波数(wave_numは実際に出力する波数)
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, ny, nz, i;                                  // ループ用変数
  int nMax = 1600;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kx, ky について，0 から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= nMax; nx++) {
      k[0] = 2.0 * M_PI * (double)nx / (double)nMax - M_PI;

      for(ny = 0; ny <= nMax; ny++) {
        k[1] = 2.0 * M_PI * (double)ny / (double)nMax - M_PI;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

        // 対角化する
        info = zmatrix.zheev();

        // if(zmatrix.eigen(DIM / 2) <= 10.0 * real(order[0])) {
          // 波数を設定する(piを超えたときは折り返す)
          for (i = 0; i < 3; i++) {
            if(k(i) + q(i) * 0.5 <= M_PI) {
              wave_num(i) = k(i) + q(i) * 0.5;
            }
            else {
              wave_num(i) = k(i) + q(i) * 0.5 - 2.0 * M_PI;
            }
          }

          // 波数を出力する
          eigen_values[nz] << fixed << setprecision(8) << wave_num[0] << "  " << wave_num[1] << "  " << wave_num[2];

          // エネルギー固有値を出力する
          eigen_values[nz] << fixed << setprecision(8) << "  " << zmatrix.eigen(DIM / 2) << endl;
        // }
      }
      eigen_values[nz] << endl;
    }

    eigen_values[nz].close();
  }

}

void Hamiltonian::output_eigen_kykz(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3), wave_num(3);                          // 波数(wave_numは実際に出力する波数)
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, ny, nz, i;                                  // ループ用変数
  int nMax = 800;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kxzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kxpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(nx = 0; nx < 2; nx++) {
    k[0] = M_PI * nx;

    // ky, kz について，0 から pi まで波数を変えて各点で固有値を計算
    for(ny = 0; ny <= nMax; ny++) {
      k[1] = M_PI * (double)ny / (double)nMax;

      for(nz = 0; nz <= nMax; nz++) {
        k[2] = M_PI * (double)nz / (double)nMax;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

        // 対角化する
        info = zmatrix.zheev();

        // if(zmatrix.eigen(DIM / 2) <= 10.0 * real(order[0])) {
          // 波数を設定する(piを超えたときは折り返す)
          for (i = 0; i < 3; i++) {
            if(k(i) + q(i) * 0.5 <= M_PI) {
              wave_num(i) = k(i) + q(i) * 0.5;
            }
            else {
              wave_num(i) = k(i) + q(i) * 0.5 - 2.0 * M_PI;
            }
          }

          // 波数を出力する
          eigen_values[nx] << fixed << setprecision(8) << wave_num[0] << "  " << wave_num[1] << "  " << wave_num[2];

          // エネルギー固有値を出力する
          eigen_values[nx] << fixed << setprecision(8) << "  " << zmatrix.eigen(DIM / 2) << endl;
        // }
      }
      eigen_values[nx] << endl;
    }

    eigen_values[nx].close();
  }

}

void Hamiltonian::output_eigen_kzkx(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3), wave_num(3);                          // 波数(wave_numは実際に出力する波数)
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nz, nx, ny, i;                                  // ループ用変数
  int nMax = 800;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kyzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kypi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(ny = 0; ny < 2; ny++) {
    k[1] = M_PI * ny;

    // kz, kx について，0 から pi まで波数を変えて各点で固有値を計算
    for(nz = 0; nz <= nMax; nz++) {
      k[2] = M_PI * (double)nz / (double)nMax;

      for(nx = 0; nx <= nMax; nx++) {
        k[0] = M_PI * (double)nx / (double)nMax;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

        // 対角化する
        info = zmatrix.zheev();

        // if(zmatrix.eigen(DIM / 2) <= 10.0 * real(order[0])) {
          // 波数を設定する(piを超えたときは折り返す)
          for (i = 0; i < 3; i++) {
            if(k(i) + q(i) * 0.5 <= M_PI) {
              wave_num(i) = k(i) + q(i) * 0.5;
            }
            else {
              wave_num(i) = k(i) + q(i) * 0.5 - 2.0 * M_PI;
            }
          }

          // 波数を出力する
          eigen_values[ny] << fixed << setprecision(8) << wave_num[0] << "  " << wave_num[1] << "  " << wave_num[2];

          // エネルギー固有値を出力する
          eigen_values[ny] << fixed << setprecision(8) << "  " << zmatrix.eigen(DIM / 2) << endl;
        // }
      }
      eigen_values[ny] << endl;
    }

    eigen_values[ny].close();
  }

}

//
// 固有値をファイルに出力する(E = 0に近いところのみ)
//   t: 飛び移り積分
//   order: 秩序パラメータ
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//
void Hamiltonian::output_eigen2(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3), wave_num(3);                          // 波数(wave_numは実際に出力する波数)
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nz, nx, ny, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-s_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
  eigen_values.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // kx, ky, kz について，0 から pi まで波数を変えて各点で固有値を計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[2]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

        // 対角化する
        info = zmatrix.zheev();

        // 波数を設定する(piを超えたときは折り返す)
        for (i = 0; i < 3; i++) {
          if(k(i) + q(i) * 0.5 <= M_PI) {
            wave_num(i) = k(i) + q(i) * 0.5;
          }
          else {
            wave_num(i) = k(i) + q(i) * 0.5 - 2.0 * M_PI;
          }
        }

        // 波数を出力する
        eigen_values << fixed << setprecision(8) << wave_num[0] << "  " << wave_num[1] << "  " << wave_num[2];

        // エネルギー固有値を出力する
        eigen_values << fixed << setprecision(8) << "  " << zmatrix.eigen(DIM / 2) << endl;
      }
    }
    eigen_values << endl;
  }

  eigen_values.close();
}

//
// バンド基底で表示したギャップの大きさを出力する
//   t: 飛び移り積分
//   order: 秩序パラメータ
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//
void Hamiltonian::output_bandgap(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nz, nx, ny, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値
  Cmatrix calc =
    zero_matrix<Complex>(DIM / 2, DIM / 2);           // 計算用
  Cmatrix calc2 = zero_matrix<Complex>(2, 2);         // 計算用(2*2行列)

  // ファイル出力関連
  stringstream filename;
  ofstream gap;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/bandgap_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
  gap.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // kx, ky, kz について，-pi から pi まで波数を変えて各点で固有値を計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[2]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // 波数を出力する
        gap << fixed << setprecision(8) << k(0) << "  " << k(1) << "  " << k(2);

        // Normal part Hamiltonianを定義し対角化する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
        info = zmatrix.zheev_normal();

        // バンド基底で表示した秩序変数を計算する
        calc = prod(herm(zmatrix.hamil_normal), zmatrix.order_part(k, order));
        k = - k;
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
        info = zmatrix.zheev_normal();
        calc = prod(calc, conj(zmatrix.hamil_normal));

        // Kramersの縮退を考慮して2*2ずつ行列をとり，ゲージ不変な量を計算・出力
        for(i = 0; i < DIM / 4; i++) {
          calc2(0, 0) = calc(2 * i, 2 * i);
          calc2(0, 1) = calc(2 * i, 2 * i + 1);
          calc2(1, 0) = calc(2 * i + 1, 2 * i);
          calc2(1, 1) = calc(2 * i + 1, 2 * i + 1);
          calc2 = prod(herm(calc2), calc2);
          gap << fixed << setprecision(8) << "  " << sqrt(0.5 * abs(calc2(0, 0) + calc2(1, 1)));
        }
        gap << endl;

        // 波数を戻しておく
        k = - k;
      }
    }
    gap << endl;
  }

  gap.close();

}

void Hamiltonian::test_kx(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, ny, i;                                      // ループ用変数
  int nyMax = 4000;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kxzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kxpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  k[2] = 0.0;
  for(nx = 0; nx < 2; nx++) {
    k[0] = M_PI * nx;

    // ky について，0 から pi まで波数を変えて各点で固有値を計算
    for(ny = 0; ny <= nyMax; ny++) {
      k[1] = M_PI * (double)ny / (double)nyMax;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

      // 対角化する
      info = zmatrix.zheev();

      // 波数を出力する
      eigen_values[nx] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

      // エネルギー固有値を出力する
      for(i = 0; i < DIM; i++) {
        eigen_values[nx] << fixed << setprecision(8) << "  " << zmatrix.eigen(i);
      }

      // // Normal part Hamiltonianを定義・対角化する
      // zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
      // info = zmatrix.zheev_normal();

      // // エネルギー固有値を出力する
      // for(i = 0; i < DIM / 2; i++) {
      //   eigen_values[nx] << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
      // }

      eigen_values[nx] << endl;
    }

    eigen_values[nx].close();
  }

}

void Hamiltonian::test_ky(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, ny, i;                                      // ループ用変数
  int nxMax = 4000;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kyzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kypi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  k[2] = 0.0;
  for(ny = 0; ny < 2; ny++) {
    k[1] = M_PI * ny;

    // kx について，0 から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= nxMax; nx++) {
      k[0] = M_PI * (double)nx / (double)nxMax;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

      // 対角化する
      info = zmatrix.zheev();

      // 波数を出力する
      eigen_values[ny] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

      // エネルギー固有値を出力する
      for(i = 0; i < DIM; i++) {
        eigen_values[ny] << fixed << setprecision(8) << "  " << zmatrix.eigen(i);
      }

      // // Normal part Hamiltonianを定義・対角化する
      // zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
      // info = zmatrix.zheev_normal();

      // // エネルギー固有値を出力する
      // for(i = 0; i < DIM / 2; i++) {
      //   eigen_values[ny] << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
      // }

      eigen_values[ny] << endl;
    }

    eigen_values[ny].close();
  }

}

void Hamiltonian::test_kz(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, nz, i;                                      // ループ用変数
  int nxMax = 4000;
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/test_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  k[1] = 0.5 * M_PI;
  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kx について，0 から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= nxMax; nx++) {
      k[0] = M_PI * (double)nx / (double)nxMax;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

      // 対角化する
      info = zmatrix.zheev();

      // 波数を出力する
      eigen_values[nz] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

      // エネルギー固有値を出力する
      for(i = 0; i < DIM; i++) {
        eigen_values[nz] << fixed << setprecision(8) << "  " << zmatrix.eigen(i);
      }

      // // Normal part Hamiltonianを定義・対角化する
      // zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
      // info = zmatrix.zheev_normal();

      // // エネルギー固有値を出力する
      // for(i = 0; i < DIM / 2; i++) {
      //   eigen_values[nz] << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
      // }

      eigen_values[nz] << endl;
    }

    eigen_values[nz].close();
  }

}

//
// ギャップの角度依存性を出力する
//   t: 飛び移り積分
//   order: 秩序パラメータ
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//   q: Cooper対の重心運動量
//
void Hamiltonian::output_gap_angle(Dvector &t, Cvector &order, double h, Dvector &theta, Dvector &alpha, double mu, Dvector &q) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nz, nr, nt, i;                                  // ループ用変数
  long nrMax = 4000;                                  // nrの最大値
  long ntMax = 100;                                   // ntの最大値
  double Emin;                                        // エネルギーの最小値
  double kr, kt;                                      // 波数の動径成分と角度成分
  long info;                                          // 対角化zheevの関数値
  Cmatrix calc =
    zero_matrix<Complex>(DIM / 2, DIM / 2);           // 計算用
  Cmatrix calc2 = zero_matrix<Complex>(2, 2);         // 計算用(2*2行列)

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/gap_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  filename << fixed << setprecision(2) << "./data/gap_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kr, kt について，波数を変えて各点で固有値を計算
    for(nt = 0; nt <= ntMax; nt++) {
      kt = M_PI * (double)nt / (double)(2 * ntMax);
      Emin = 100.0;

      for(nr = 0; nr <= nrMax; nr++) {
        k[0] = M_PI * (double)nr / (double)nrMax * cos(kt);
        k[1] = M_PI * (double)nr / (double)nrMax * sin(kt);

        // Hamiltonianを定義する
        zmatrix.init_zmatrix(k, t, order, h, theta, alpha, mu, q);

        // 対角化する
        info = zmatrix.zheev();

        // 最小のエネルギー固有値(ギャップの大きさ)のみを抽出
        if(zmatrix.eigen(DIM / 2) < Emin) {
          Emin = zmatrix.eigen(DIM / 2);
          kr = M_PI * (double)nr / (double)nrMax;
        }
      }

      // 波数とエネルギー固有値を出力する
      eigen_values[nz] << fixed << setprecision(8) << kt << "  " << kr << "  " << Emin;

      // Normal part Hamiltonianを定義し対角化する
      k[0] = kr * cos(kt);
      k[1] = kr * sin(kt);
      zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
      info = zmatrix.zheev_normal();

      // バンド基底で表示した秩序変数を計算する
      calc = prod(herm(zmatrix.hamil_normal), zmatrix.order_part(k, order));
      k = - k;
      zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
      info = zmatrix.zheev_normal();
      calc = prod(calc, conj(zmatrix.hamil_normal));

      // Kramersの縮退を考慮して2*2ずつ行列をとり，ゲージ不変な量を計算・出力
      for(i = 0; i < DIM / 4; i++) {
        calc2(0, 0) = calc(2 * i, 2 * i);
        calc2(0, 1) = calc(2 * i, 2 * i + 1);
        calc2(1, 0) = calc(2 * i + 1, 2 * i);
        calc2(1, 1) = calc(2 * i + 1, 2 * i + 1);
        calc2 = prod(herm(calc2), calc2);
        eigen_values[nz] << fixed << setprecision(8) << "  " << sqrt(0.5 * abs(calc2(0, 0) + calc2(1, 1)));
      }
      eigen_values[nz] << endl;
    }

    eigen_values[nz].close();
  }

}

//
// Normal part Hamiltonianの固有値をファイルに出力する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//
void Hamiltonian::output_eigen_normal(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
  Dvector k(3);                                       // 波数
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  int nx, ny, nz, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
  eigen_values.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // kx, ky, kz について，-pi から pi まで波数を変えて各点で固有値を計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[2]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // Hamiltonianを定義・対角化する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);
        info = zmatrix.zheev_normal();

        // 波数を出力する
        eigen_values << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

        // エネルギー固有値を出力する
        for(i = 0; i < DIM / 4; i++) {
          eigen_values << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(2 * i);
        }
        eigen_values << endl;

      }
    }
    eigen_values << endl;

  }
  eigen_values.close();

}
