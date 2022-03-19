#pragma once
#include<vector>
#include<iostream>
#include "Matrix.h"
using namespace std;


void print(Matrix<double> matrix);
void print(const vector<double>& matrix, int n, int m);
void print(vector<int> vec);

class SimplexMethod {
private:
    Matrix<double> A;
    Matrix<double> b;
    Matrix<double> c;
    Matrix<double> B_Nk;
    Matrix<double> d_Nk;
    vector<int> Nk_ind;
    vector<int> Lk_ind;
    vector<int> M_ind;
    vector<int> N_ind;
    vector<bool> verify;
    int index_jk;
    int index_i;
    Matrix<double> referenceVector;
    void add();
    bool calcD_1();
    bool calcD_2();
    void set_new_B();
    double theta(Matrix<double> a, Matrix<double>b);
    bool degenerate(Matrix<double> vec, Matrix<double> u);
public:
    SimplexMethod(const vector<double> A1, const  vector<double> b1, const vector<double>c1, const  int M1, const  int N1);
    void run();
};
void SimplexMethod::run() {
    add();
    bool a, b, c;
    int iter = 0;
    while (true) {
        a = calcD_1();
        if (a == false) {
            cout << "solve" << endl;;
            print(referenceVector);
            break;
        }
        c = calcD_2();
        set_new_B();
        iter++;
    }
}

void print(Matrix<double> matrix) {
    for (int j = 0; j < matrix.getLines(); j++) {
        for (int i = 0; i < matrix.getColumns(); i++) {
            cout << matrix.getValue(i, j) << ' ';
        }
        cout << endl;
    }
}
void print(const vector<double>& matrix, int n, int m)
 {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i * n + j] << ' ';
        }
        cout << endl;
    }
}
bool degenerate(Matrix<double> vec, int m) {
    int sum_ = 0;
    for (int i = 0; i < vec.getLines(); i++) {
        if (vec.getValue(0, i) > 0) {
            sum_ += 1;
        }
    }
    return !(sum_ == m);
}
void sort_ind(vector<int>& vec, int k) {
    int c = vec[k];
    vec.erase(vec.begin() + k);
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] > c) {
            vec.insert(vec.begin() + i, c);
            return;
        }
    }
    vec.push_back(c);
}
int search(vector<int>& vec, int value) {
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] == value) {
            return i;
        }
    }
}
void print(vector<int> vec) {
    for (int i : vec) {
        cout << i << ' ';
    }
    cout << endl;
}



void SimplexMethod::set_new_B() {
    Matrix<double> un = B_Nk * A.Samp({ index_jk }, M_ind);
    Matrix<double> F(Nk_ind.size(), Nk_ind.size());
    for (int i = 0; i < Nk_ind.size(); i++) {
        F.setValue(1, i, i);
    }
    int index = 0;
    for (int i = 0; i < Nk_ind.size(); i++) {
        if (referenceVector.getValue(0, Nk_ind[i]) == 0) {
            index = i;
            break;
        }
    }
    for (int i = 0; i < Nk_ind.size(); i++) {
        F.setValue(un.getValue(0, i), index, i);
    }
    double num = F.getValue(index, index);
    for (int i = 0; i < Nk_ind.size(); i++) {
        if (i != index) {
            F.setValue(-F.getValue(index, i) / num, index, i);
        }
        else {
            F.setValue(1 / num, index, index);
        }
    }
    int swp = Nk_ind[index];
    Nk_ind[index] = index_jk;
    int ind = search(Lk_ind, index_jk);
    Lk_ind[ind] = swp;
    //sort_ind(Lk_ind, ind);
    //sort_ind(Nk_ind, index);
    B_Nk = F * B_Nk;
}
bool SimplexMethod::calcD_2() {
    Matrix<double> u_n(N_ind.size(), 1);
    Matrix<double> u_nk = B_Nk * A.Samp({ index_jk }, M_ind);
    for (int i = 0; i < Nk_ind.size(); i++) {
        u_n.setValue(u_nk.getValue(0, i), 0, Nk_ind[i]);
    }
    u_n.setValue(-1, 0, index_jk);
    if (!degenerate(referenceVector, u_n)) {
        return false;
    }
    double th = theta(referenceVector, u_n);
    referenceVector = referenceVector - u_n.multiply(th);
    return true;
}
bool SimplexMethod::calcD_1() {
    Matrix<double> c_Nk = c.Samp({ 0 }, Nk_ind);
    Matrix<double> d_N = (c.Tr() - c_Nk.Tr() * B_Nk * A).Tr();
    d_Nk = d_N;

    index_jk = -1;
    for (int j = 0; j < Lk_ind.size(); j++) {
        if (d_Nk.getValue(0, Lk_ind[j]) < 0) {
            index_jk = Lk_ind[j];
            return true;
        }
    }
    return false;
}
void SimplexMethod::add() {
    int N = A.getColumns();
    int M = A.getLines();
    vector<double> new_A(M * (N + M));
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N + M; j++) {
            if (j < N) {
                new_A[i * (N + M) + j] = A.getValue(j, i);
            }
            else {
                new_A[i * (N + M) + j] = 0;
            }
        }
    }
    for (int i = 0; i < M; i++) {
        new_A[(N + M) * i + N + i] = 1;
    }
    A.setMatrix(new_A, M, N + M);
    vector<double> refVec(N + M, 0);

    for (int i = N; i < N + M; i++) {
        refVec[i] = b.getValue(0, i - N);
    }
    referenceVector.setMatrix(refVec, N + M, 1);
    vector<double> new_vec_c(M + N, 1);
    for (int i = 0; i < N; i++) {
        new_vec_c[i] = c.getValue(0, i);
    }
    c.setMatrix(new_vec_c, N + M, 1);
    vector<int> new_NK(M);
    for (int i = N; i < N + M; i++) {
        new_NK[i - N] = i;
    }
    vector<int> new_LK(N);
    for (int i = 0; i < N; i++) {
        new_LK[i] = i;
    }

    vector<int> M_(A.getLines());
    for (int i = 0; i < A.getLines(); i++) {
        M_[i] = i;
    }
    M_ind = M_;
    vector<int> N_(A.getColumns());
    for (int i = 0; i < A.getColumns(); i++) {
        N_[i] = i;
    }
    N_ind = N_;

    Nk_ind = new_NK;
    Lk_ind = new_LK;
    B_Nk = A.Samp(Nk_ind, M_ind);
}
double SimplexMethod::theta(Matrix<double> a, Matrix<double>b) {
    double min_ = INT_MAX;
    double min_2;
    for (int i = 0; i < Nk_ind.size(); i++) {
        for (int j = 0; j < Nk_ind.size(); j++) {
            if (b.getValue(0, Nk_ind[i]) > 0) {
                min_ = min(min_, a.getValue(0, Nk_ind[i]) / b.getValue(0, Nk_ind[i]));
            }
        }
    }
    return min_;
}
bool SimplexMethod::degenerate(Matrix<double> vec, Matrix<double> u) {
    int sum_ = 0;
    for (int i = 0; i < Nk_ind.size(); i++) {
        if (vec.getValue(0, Nk_ind[i]) == 0 && u.getValue(0, Nk_ind[i]) > 0) {
            return false;
        }
    }
    return true;
}
SimplexMethod::SimplexMethod(const vector<double> A1, const  vector<double> b1, const vector<double>c1, const  int M1, const  int N1) {
    A.setMatrix(A1, M1, N1);
    b.setMatrix(b1, M1, 1);
    c.setMatrix(c1, M1, 1);
    verify = vector<bool>(N1, 0);
}