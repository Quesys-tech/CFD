#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include <cmath>
using namespace std;

#include "QLA/include/qla.hpp"
using namespace QLA;

const int n = 50;
constexpr double h = 1 / (double)n;
const double dt = 1e-4;
const double re = 40;

//Poisson 方程式のソルバー
//http://nalab.mind.meiji.ac.jp/~mk/labo/text/poisson.pdf
// http://www.damtp.cam.ac.uk/user/reh10/lectures/nst-mmii-chapter2.pdf p.23
//のアルゴリズムを利用
void Poisson_solver_SQlattice(
    Matrix &u,         //離散化された支配方程式が成立する場
    Matrix const &f,   //poisson方程式の離散化した右辺
    double _h,         //格子の感覚
    double omg = 1.2,  // 加速係数 1 omg<2
    double eps = 1e-10 //終了条件
)
{
    unsigned int m = u.getCols();
    assert(m == u.getRows());
    assert(m == f.getRows());

    const double h2 = _h * _h; //h^2
    unsigned int count_iter = 0;
    while (true)
    {
        double err = 0;

        for (int i = 1; i < m - 1; i++)
        {
            for (int j = 1; j < m - 1; j++)
            {
                double y = (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1) - h2 * f(i, j)) / 4.0;
                double err_ij = fabs(omg * (y - u(i, j)));
                u(i, j) = u(i, j) + omg * (y - u(i, j));

                err = max(err_ij, err);
            }
        }

        if (err <= eps)
        {
            return;
        }
        count_iter++;

        if (count_iter > 1e9)
        {
            cerr << "failed to solve" << endl;
            return;
        }
    }
}


void showProgress(float progress)
{
    int width = 70;
    std::cout << "[";

    int pos = (int)round(progress * (float)width);
    for (size_t i = 0; i < width; i++)
    {
        if (i < pos)
        {
            std::cout << "#";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "]" << (int)round(progress * 100.0) << "%\r";
    std::cout.flush();
}

int main(int argc, char const *argv[])
{
    Matrix psi(n + 1, n + 1, 0);  //流れ関数\psi
    Matrix omg(n + 1, n + 1);     //渦度 omg
    Matrix omg_new(n + 1, n + 1); //次の時刻の omg

    //流れ関数の境界条件
    //すべての要素が0で初期化されているので不要

    double max_t = 10; //計算する時刻の最大値
    double t = 0;      //現在時刻
    for (int n = 0; n * dt <= max_t; n++)
    {

        //渦度の境界条件の設定
        for (size_t i = 0; i < omg.getCols(); i++) //左右
        {
            omg(0, i) = -2.0 * psi(1, i) / (h * h);
            omg(omg.getRows() - 1, i) = -2.0 * psi(omg.getRows() - 1, i) / (h * h);
            //CD上の境界条件
        }
        for (size_t i = 0; i < omg.getCols(); i++) //上下
        {
            omg(i, 0) = -2.0 * psi(i, 1) / (h * h);
            omg(i, omg.getCols() - 1) = -2.0 * (psi(i, omg.getCols() - 1) + h) / (h * h);
        }

//渦度の計算
#ifdef WITH_OMP
#pragma omp parallel for
#endif
        for (int i = 1; i < omg.getRows() - 1; i++)
        {
            for (int j = 1; j < omg.getCols() - 1; j++)
            {
                double h2 = h * h; //h^2
                double rhs = ((psi(i + 1, j) - psi(i - 1, j)) * (omg(i, j + 1) - omg(i, j - 1)) / 4.0 - (psi(i, j + 1) - psi(i, j - 1)) * (omg(i + 1, j) - omg(i - 1, j))) / 4.0 + (omg(i - 1, j) + omg(i + 1, j) + omg(i, j + 1) + omg(i, j - 1) - 4.0 * omg(i, j)) / re;
                omg_new(i, j) = omg(i, j) + dt / h2 * rhs;
            }
        }
        omg.copy(omg_new);

        Poisson_solver_SQlattice(psi, -1.0 * omg, h, 1.2); //流れ関数の計算

        t = (double)(n + 1) * dt; //時刻を1ステップ進める

        showProgress(t / max_t); //進捗の表示

        if (n % 10000 == 0)
        {
            ofstream ofs(string("results/psi_") + to_string(t) + string(".txt"));
            if (ofs.fail())
            {
                cerr << "failed to write at" << string("results/psi_") + to_string(t) + string(".txt") << endl;
            }
            for (int i = 1; i < psi.getRows() - 1; i++)
            {
                for (int j = 1; j < psi.getCols() - 1; j++)
                {
                    ofs << (double)i * h << " " << (double)j * h << " " << psi(i, j) << endl;
                }
                ofs << endl;
            }
        }
    }
    cout << endl;
    return 0;
}