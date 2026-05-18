#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace chrono;

const int N = 1024;

class Matrix {
public:
    vector<float> data;
    int n;

    Matrix(int size) : n(size), data(size* size) {}

    float& operator()(int i, int j) { return data[i * n + j]; }
    const float& operator()(int i, int j) const { return data[i * n + j]; }

    void fillRandom() {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<float> dist(-10.0f, 10.0f);
        for (auto& val : data) {
            val = dist(gen);
        }
    }

    void fillZero() {
        fill(data.begin(), data.end(), 0.0f);
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

void multiplyOptimized(Matrix& C, const Matrix& A, const Matrix& B) {
    C.fillZero();
    const int BLOCK = 64;

    for (int i = 0; i < N; i += BLOCK) {
        int i_max = min(i + BLOCK, N);
        for (int j = 0; j < N; j += BLOCK) {
            int j_max = min(j + BLOCK, N);
            for (int k = 0; k < N; k += BLOCK) {
                int k_max = min(k + BLOCK, N);

                for (int ii = i; ii < i_max; ++ii) {
                    for (int kk = k; kk < k_max; ++kk) {
                        float aik = A(ii, kk);
                        int c_off = ii * N + j;
                        int b_off = kk * N + j;

                        for (int jj = j; jj < j_max; ++jj) {
                            C.data[c_off + jj] += aik * B.data[b_off + jj];
                        }
                    }
                }
            }
        }
    }
}

double measureTime(void (*func)(Matrix&, const Matrix&, const Matrix&),
    Matrix& C, const Matrix& A, const Matrix& B, const string& name) {
    cout << "Bыполняется: " << name << "... " << flush;

    auto start = high_resolution_clock::now();
    func(C, A, B);
    auto end = high_resolution_clock::now();

    auto duration_ms = duration_cast<milliseconds>(end - start).count();
    double duration = duration_ms / 1000.0;

    double ops = 2.0 * N * N * N;
    double mflops = (ops / duration) * 1e-6;

    cout << "вpeмя = " << fixed << setprecision(3) << duration << " c, "
        << "производительность = " << setprecision(0) << mflops << " MFLOPS" << endl;

    return mflops;
}

bool checkResult(const Matrix& C1, const Matrix& C2) {
    float maxDiff = 0;
    int errors = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float diff = fabs(C1(i, j) - C2(i, j));
            if (diff > 0.001f) {
                errors++;
                if (diff > maxDiff) maxDiff = diff;
            }
        }
    }

    if (errors == 0) {
        cout << "Peзyльтaты идентичны [OK]" << endl;
        return true;
    }
    else {
        cout << "Haйдeнo " << errors << " различий. Максимальное отклонение: " << maxDiff << endl;
        return false;
    }
}

int main() {
    system("chcp 1251 > nul");
    cout << "========================================" << endl;
    cout << "        ПЕРЕМНОЖЕНИЕ МАТРИЦ" << endl;
    cout << "            " << N << " x " << N << endl;
    cout << "========================================" << endl << endl;

    cout << "[1] Инициализация матриц..." << endl;
    Matrix A(N), B(N), C_naive(N), C_opt(N);
    A.fillRandom();
    B.fillRandom();
    cout << "    Матрицы заполнены случайными числами" << endl << endl;

    cout << "[2] Запуск тестов производительности:" << endl;
    cout << "----------------------------------------" << endl;

    double mflops1 = measureTime(multiplyNaive, C_naive, A, B, "Классический алгоритм");
    double mflops2 = measureTime(multiplyOptimized, C_opt, A, B, "Оптимизированный алгоритм (блочный)");

    cout << "\n========================================" << endl;
    cout << "           РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ" << endl;
    cout << "========================================" << endl;
    cout << " Классический алгоритм:     " << fixed << setprecision(0) << mflops1 << " MFLOPS" << endl;
    cout << " Оптимизированный алгоритм: " << mflops2 << " MFLOPS" << endl;
    cout << "----------------------------------------" << endl;
    cout << " УСКОРЕНИЕ:                 x" << fixed << setprecision(2) << (mflops2 / mflops1) << endl;
    cout << "========================================" << endl;

    cout << "\n[3] ПРОВЕРКА КОРРЕКТНОСТИ:" << endl;
    cout << "----------------------------------------" << endl;
    checkResult(C_opt, C_naive);

    cout << "\n========================================" << endl;
    cout << "        ПРИМЕНЁННЫЕ ОПТИМИЗАЦИИ" << endl;
    cout << "========================================" << endl;
    cout << " * Блочная структура (tiling) с блоком 64x64" << endl;
    cout << " * Кэш-локальность (оптимизация обхода памяти)" << endl;
    cout << " * Минимизация обращений к памяти" << endl;
    cout << " * Использование регистровых переменных" << endl;
    cout << "========================================" << endl;

    cout << "\nПрограмма успешно завершена!" << endl;

    cout << "\nНажмите Enter для выхода...";
    cin.get();

    return 0;
}