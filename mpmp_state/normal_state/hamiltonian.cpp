//-*- coding: utf-8; mode: c++ -*-
// Time-stamp: <2016-11-26 17:17:29 shunta>

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
// 電子密度を計算する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//
void Hamiltonian::electron_density(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
  Dvector k(3);                                       // 波数
#ifdef _OPENMP
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
#else
  static Hamiltonian::Zmatrix zmatrix;                // 複素数行列
#endif
  int nz, nx, ny, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値
  double density = 0.0;                               // 

  // kx, ky, kz について，-pi から pi まで波数を変えて各点で固有値を計算
  for(nx = 0; nx <= N[0]; nx++) {
    k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

    for(ny = 0; ny <= N[1]; ny++) {
      k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

      for(nz = 0; nz <= N[1]; nz++) {
        k[2] = 2.0 * M_PI * (double)nz / (double)N[2] - M_PI;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);

        // 対角化する
        info = zmatrix.zheev_normal();

        // 電子密度
        for(i = 0; i < DIM / 2; i++) {
          if(zmatrix.eigen_normal(i) <= 0.0) {
            density += 1.0;
          }
        }
      }
    }
  }

  density /= (double)DIM / 4.0 * V;
  cout << "n = " << density << endl;
  
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
#ifdef _OPENMP
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
#else
  static Hamiltonian::Zmatrix zmatrix;                // 複素数行列
#endif
  int nz, nx, ny, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzero.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kpi.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kx, ky について，-pi から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= N[0]; nx++) {
      k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

      for(ny = 0; ny <= N[1]; ny++) {
        k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

        // Hamiltonianを定義する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);

        // 対角化する
        info = zmatrix.zheev_normal();

        // 波数を出力する
        eigen_values[nz] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

        // エネルギー固有値を出力する
        for(i = 0; i < DIM / 2; i++) {
          eigen_values[nz] << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
        }
        eigen_values[nz] << endl;

      }
      eigen_values[nz] << endl;
    }

    eigen_values[nz].close();
  }

}

//
// Normal part Hamiltonianの(kxだけ動かした)固有値をファイルに出力する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//   mu: 化学ポテンシャル
//
void Hamiltonian::output_eigen_normal_kx(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
  Dvector k(3);                                       // 波数
#ifdef _OPENMP
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
#else
  static Hamiltonian::Zmatrix zmatrix;                // 複素数行列
#endif
  int nz, nx, i;                                      // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream eigen_values[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzero2.d";
  eigen_values[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kpi2.d";
  eigen_values[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  k[1] = 0.0;
  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kx について，-pi から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= N[0]; nx++) {
      k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

      // Hamiltonianを定義する
      zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);

      // 対角化する
      info = zmatrix.zheev_normal();

      // 波数を出力する
      eigen_values[nz] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

      // エネルギー固有値を出力する
      for(i = 0; i < DIM / 2; i++) {
        eigen_values[nz] << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
      }
      eigen_values[nz] << endl;
    }

    eigen_values[nz].close();
  }

}

//
// Normal part Hamiltonianの代表的な固有値をファイルに出力する
// (Γ -> X -> S -> -> Y -> Γ -> Z -> U -> R -> T -> Z -> U -> X -> S -> R -> T -> Y)
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//
void Hamiltonian::output_eigen_normal2(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
  static Dvector k = zero_vector<double>(3);          // 波数
  static double knum = 1.0;

  // ファイル出力関連
  stringstream filename;
  static ofstream eigen_values;

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/ev-n_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << ".d";
  eigen_values.open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 繰り返し使うローカル関数を定義
  struct {
    void operator()(int coord, bool pm, Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu){
#ifdef _OPENMP
      Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
#else
      static Hamiltonian::Zmatrix zmatrix;                // 複素数行列
#endif
      int n, i;                                           // ループ用変数
      long info;                                          // 対角化zheevの関数値

      for(n = 1; n <= N[0]; n++) {
        if (pm) {
          k[coord] = M_PI * (double)n / (double)N[0];
        }
        else {
          k[coord] = M_PI * (double)(N[0] - n) / (double)N[0];
        }

        // Hamiltonianを定義する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);

        // 対角化する
        info = zmatrix.zheev_normal();

        // 波数を出力する
        eigen_values << fixed << setprecision(8) << knum;

        // エネルギー固有値を出力する
        for(i = 0; i < DIM / 2; i++) {
          eigen_values << fixed << setprecision(8) << "  " << zmatrix.eigen_normal(i);
        }
        eigen_values << endl;

        knum += 1.0;
      }
    }
  } output;


  cout << "Gamma: " << k << endl;

  // Γ点 -> X点
  output(0, true, t, h, theta, alpha, mu);
  cout << "X: " << k << endl;

  // X点 -> S点
  output(1, true, t, h, theta, alpha, mu);
  cout << "S: " << k << endl;

  // S点 -> Y点
  output(0, false, t, h, theta, alpha, mu);
  cout << "Y: " << k << endl;

  // Y点 -> Γ点
  output(1, false, t, h, theta, alpha, mu);
  cout << "Gamma: " << k << endl;

  // Γ点 -> Z点
  output(2, true, t, h, theta, alpha, mu);
  cout << "Z: " << k << endl;

  // Z点 -> U点
  output(0, true, t, h, theta, alpha, mu);
  cout << "U: " << k << endl;

  // U点 -> R点
  output(1, true, t, h, theta, alpha, mu);
  cout << "R: " << k << endl;

  // R点 -> T点
  output(0, false, t, h, theta, alpha, mu);
  cout << "T: " << k << endl;

  // T点 -> Z点
  output(1, false, t, h, theta, alpha, mu);
  cout << "Z: " << k << endl;

  // Z点 -> U点
  output(0, true, t, h, theta, alpha, mu);
  cout << "U: " << k << endl;

  // U点 -> X点
  output(2, false, t, h, theta, alpha, mu);
  cout << "X: " << k << endl;

  // X点 -> S点
  output(1, true, t, h, theta, alpha, mu);
  cout << "S: " << k << endl;

  // S点 -> R点
  output(2, true, t, h, theta, alpha, mu);
  cout << "R: " << k << endl;

  // R点 -> T点
  output(0, false, t, h, theta, alpha, mu);
  cout << "T: " << k << endl;

  // T点 -> Y点
  output(2, false, t, h, theta, alpha, mu);
  cout << "Y: " << k << endl << endl;

  eigen_values.close();

}

//
// Normal part energyの非対称部分をファイルに出力する
//   t: 飛び移り積分
//   h: 磁気モーメントの大きさ
//   theta: 磁気モーメントの角度
//   alpha: 反対称スピン軌道相互作用の強さ
//
void Hamiltonian::output_asym(Dvector &t, double h, Dvector &theta, Dvector &alpha, double mu) {
  Dvector k(3);                                       // 波数
#ifdef _OPENMP
  Hamiltonian::Zmatrix zmatrix;                       // 複素数行列
  Hamiltonian::Zmatrix zmatrix_op;                    // 複素数行列
#else
  static Hamiltonian::Zmatrix zmatrix;                // 複素数行列
  static Hamiltonian::Zmatrix zmatrix_op;             // 複素数行列
#endif
  int nz, nx, ny, i;                                  // ループ用変数
  long info;                                          // 対角化zheevの関数値

  // ファイル出力関連
  stringstream filename;
  ofstream asym[2];

  // データを出力するためのフォルダを作成
  mkdir("./data", 0755);

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/asym_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kzero.d";
  asym[0].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  // 出力ファイルを作成
  filename << fixed << setprecision(2) << "./data/asym_h" << h
           << "_alpha" << alpha[0] << "_" << alpha[1] << "_" << alpha[2] << "_mu" << mu << "_kpi.d";
  asym[1].open(filename.str().c_str());
  filename.str("");                      // バッファのクリア
  filename.clear(stringstream::goodbit); // フラグのクリア

  for(nz = 0; nz < 2; nz++) {
    k[2] = M_PI * nz;

    // kx, ky について，-pi から pi まで波数を変えて各点で固有値を計算
    for(nx = 0; nx <= N[0]; nx++) {
      k[0] = 2.0 * M_PI * (double)nx / (double)N[0] - M_PI;

      for(ny = 0; ny <= N[1]; ny++) {
        k[1] = 2.0 * M_PI * (double)ny / (double)N[1] - M_PI;

        // 波数を出力する
        asym[nz] << fixed << setprecision(8) << k[0] << "  " << k[1] << "  " << k[2];

        // Hamiltonianを定義する
        zmatrix.init_zmatrix_normal(k, t, h, theta, alpha, mu);

        // 対角化する
        info = zmatrix.zheev_normal();

        // kxの符号を反転する
        k[0] = - k[0];

        // Hamiltonianを定義する
        zmatrix_op.init_zmatrix_normal(k, t, h, theta, alpha, mu);

        // 対角化する
        info = zmatrix_op.zheev_normal();

        // エネルギー固有値の差を出力する
        for(i = 0; i < DIM / 2; i++) {
          asym[nz] << fixed << setprecision(8) << "  " << (zmatrix.eigen_normal(i) - zmatrix_op.eigen_normal(i));
        }
        asym[nz] << endl;

        // kxの符号を戻しておく
        k[0] = - k[0];
      }
      asym[nz] << endl;
    }

    asym[nz].close();
  }

}
