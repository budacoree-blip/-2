#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <cmath>
#include <windows.h>

using namespace std;

const int N = 1024;
const long long TOTAL_OPS = 2LL * N * N * N;

struct Matrix {
    float* data;
    int n;

    Matrix(int size) : n(size) {
        data = new float[size * size];
    }

    ~Matrix() {
        delete[] data;
    }

    float& operator()(int i, int j) {
        return data[i * n + j];
    }

    const float& operator()(int i, int j) const {
        return data[i * n + j];
    }

    void fillRandom() {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<float> dist(-1.0f, 1.0f);
        for (int i = 0; i < n * n; ++i) {
            data[i] = dist(gen);
        }
    }

    void fillZero() {
        for (int i = 0; i < n * n; ++i) {
            data[i] = 0.0f;
        }
    }
};

void multiplyNaive(Matrix& C, const Matrix& A, const Matrix& B) {
    C.fillZero();
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < N; ++k) {
            float aik = A(i, k);
            for (int j = 0; j < N; ++j) {
                C(i, j) += aik * B(k, j);
            }
        }
    }
}

void multiplyBLAS(Matrix& C, const Matrix& A, const Matrix& B) {
    C.fillZero();
    const int BLOCK = 32;

    for (int i = 0; i < N; i += BLOCK) {
        int i_max = min(i + BLOCK, N);
        for (int j = 0; j < N; j += BLOCK) {
            int j_max = min(j + BLOCK, N);
            for (int k = 0; k < N; k += BLOCK) {
                int k_max = min(k + BLOCK, N);

                for (int ii = i; ii < i_max; ++ii) {
                    for (int kk = k; kk < k_max; ++kk) {
                        float aik = A(ii, kk);
                        for (int jj = j; jj < j_max; ++jj) {
                            C(ii, jj) += aik * B(kk, jj);
                        }
                    }
                }
            }
        }
    }
}

void multiplyOptimized(Matrix& C, const Matrix& A, const Matrix& B) {
    C.fillZero();
    const int BLOCK = 64;

    for (int i = 0; i < N; i += BLOCK) {
        int i_end = min(i + BLOCK, N);
        for (int j = 0; j < N; j += BLOCK) {
            int j_end = min(j + BLOCK, N);
            for (int k = 0; k < N; k += BLOCK) {
                int k_end = min(k + BLOCK, N);

                for (int ii = i; ii < i_end; ++ii) {
                    for (int kk = k; kk < k_end; ++kk) {
                        float aik = A(ii, kk);
                        for (int jj = j; jj < j_end; ++jj) {
                            C(ii, jj) += aik * B(kk, jj);
                        }
                    }
                }
            }
        }
    }
}

double measureTime(void (*func)(Matrix&, const Matrix&, const Matrix&),
    Matrix& C, const Matrix& A, const Matrix& B, const string& name) {

    cout << "Выполняется: " << name << "... " << flush;

    LARGE_INTEGER frequency, start, end;
    QueryPerformanceFrequency(&frequency);

    QueryPerformanceCounter(&start);
    func(C, A, B);
    QueryPerformanceCounter(&end);

    double duration = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;
    double mflops = (TOTAL_OPS / duration) * 1e-6;

    cout << "время = " << fixed << setprecision(2) << duration << " с, "
        << "производительность = " << setprecision(0) << mflops << " MFLOPS" << endl;

    return mflops;
}

void checkResult(const Matrix& C1, const Matrix& C2, const string& name1, const string& name2) {
    float maxDiff = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float diff = fabs(C1(i, j) - C2(i, j));
            if (diff > maxDiff) maxDiff = diff;
        }
    }

    if (maxDiff < 0.01f) {
        cout << "  [OK] " << name1 << " и " << name2 << " совпадают" << endl;
    }
    else {
        cout << "  [WARNING] Макс. отклонение: " << maxDiff << endl;
    }
}

int main() {
    system("chcp 1251 > nul");

    cout << "========================================" << endl;
    cout << "     ПЕРЕМНОЖЕНИЕ МАТРИЦ " << N << "x" << N << endl;
    cout << "========================================" << endl;
    cout << "Число операций: " << TOTAL_OPS << endl;
    cout << "========================================" << endl << endl;

    cout << "Инициализация матриц..." << endl;
    Matrix A(N), B(N);
    A.fillRandom();
    B.fillRandom();
    cout << "Готово!" << endl << endl;

    Matrix C1(N), C2(N), C3(N);

    cout << "Запуск тестов:" << endl;
    cout << "----------------------------------------" << endl;

    double mflops1 = measureTime(multiplyNaive, C1, A, B, "1. Классический");
    double mflops2 = measureTime(multiplyBLAS, C2, A, B, "2. BLAS-подобный");
    double mflops3 = measureTime(multiplyOptimized, C3, A, B, "3. Оптимизированный");

    cout << "\n========================================" << endl;
    cout << "           РЕЗУЛЬТАТЫ" << endl;
    cout << "========================================" << endl;
    cout << " 1. Классический:         " << fixed << setprecision(0) << mflops1 << " MFLOPS" << endl;
    cout << " 2. BLAS-подобный:        " << mflops2 << " MFLOPS" << endl;
    cout << " 3. Оптимизированный:     " << mflops3 << " MFLOPS" << endl;
    cout << "----------------------------------------" << endl;

    if (mflops1 > 0) {
        cout << " Ускорение (BLAS/классический):     x" << (mflops2 / mflops1) << endl;
        cout << " Ускорение (оптимиз./классический): x" << (mflops3 / mflops1) << endl;

        double efficiency = (mflops3 / mflops2) * 100.0;
        cout << " Эффективность от BLAS:              " << efficiency << "%" << endl;

        if (efficiency >= 30.0) {
            cout << " [OK] Цель достигнута (>30% от BLAS)" << endl;
        }
    }
    cout << "========================================" << endl;

    cout << "\nПроверка корректности:" << endl;
    cout << "----------------------------------------" << endl;
    checkResult(C2, C1, "BLAS", "классическим");
    checkResult(C3, C1, "Оптимизированный", "классическим");

    cout << "\nПрограмма успешно завершена!" << endl;
    cout << "\nНажмите Enter для выхода...";
    cin.get();

    return 0;
}
