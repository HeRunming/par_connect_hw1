#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <chrono>
#include <iomanip>

using namespace std;

namespace utils {
    int N; // number of rows (and columns) in the matrix
    vector<double> mat; // the input matrix
    vector<double> L; // lower triangular matrix
    vector<double> U; // upper triangular matrix

    // ... 其他函数不变 ...

    void lu_decomposition() {
        // 开始计时
        auto start = chrono::high_resolution_clock::now();

        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];
#pragma omp parallel for collapse(2)
                for (int k = j; k < N; ++k) {
                    U[i * N + k] -= L[i * N + j] * U[j * N + k];
                }
            }
        }

        // 结束计时
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

        // 输出并行计算所用的时间
        cout << "Parallel LU decomposition took " << duration.count() << " microseconds." << endl;
    }

    void write_outputs(const string& output_filename) {
        // ... 省略 ...
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3)
        utils::abort_with_error_message("ERROR: Invalid command. Usage: ./parallel_lu <input_file> <output_file>");

    string input_filename = argv[1];
    string output_filename = argv[2];

    utils::read_inputs_and_initialize(input_filename);

    // 设置OpenMP线程数
    int num_threads = 16; // 您可以根据需要修改这个值
    omp_set_num_threads(num_threads);

    utils::lu_decomposition();

    utils::write_outputs(output_filename);

    return 0;
}
