/*
 * Foundations of Parallel Computing II, Spring 2024.
 * Instructor: Chao Yang, Xiuhong Li @ Peking University.
 * This is a serial implementation of LU decomposition.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> 

using namespace std;

namespace utils {
    int N; // number of rows (and columns) in the matrix
    vector<double> mat; // the input matrix
    vector<double> L; // lower triangular matrix
    vector<double> U; // upper triangular matrix

    // 初始化 L0 和 U0
    vector<vector<double>> L0(N, vector<double>(N));
    vector<vector<double>> U0(N, vector<double>(N));



    void abort_with_error_message(const string& msg) {
        cerr << msg << endl;
        abort();
    }
    /*
     * The input file will be in the following format:
     * - The first line is an integer N, the number of rows and columns in the matrix.
     * - The following lines contain the matrix entries in row-major order.
     *
     * The output file will be in the following format:
     * - The lower triangular matrix L and the upper triangular matrix U obtained from the LU decomposition, in row-major order.
     *
     * Run:
     * $ ./serial <input_file> <output_file>
     */
    int read_inputs_and_initialize(const string& filename) {
        ifstream inputf(filename);
        if (!inputf.is_open())
            abort_with_error_message("ERROR: Unable to open input file.");

        inputf >> N;
        mat.resize(N * N);
        for (int i = 0; i < N * N; ++i) {
            inputf >> mat[i];
        }
        inputf.close();

        L.resize(N * N);
        U.resize(N * N);

        U = mat;

        // Initialize L to identity matrix
        for (int i = 0; i < N; ++i) {
            L[i * N + i] = 1.0;
        }
        return 0;
    }

    // 新尝试，用新的初始化方式，构建L0与U0，均为N维数组，每个维度是一个N维数组，代表L或U的一列
    int init_L0_U0() {
        L0.resize(N, vector<double>(N));
        U0.resize(N, vector<double>(N));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                L0[j][i] = L[i * N + j];
                U0[j][i] = U[i * N + j];
            }
        }
        return 0;
    }
        


    void lu_decomposition() {
        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];
                for (int k = j; k < N; ++k) {
                    U[i * N + k] -= L[i * N + j] * U[j * N + k];
                }
            }
        }
    }

        void lu_decomposition_0() {
        // 开始计时
        auto start = chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                L0[j][i] = U0[j][i] / U0[j][j];
                for (int k = j; k < N; ++k) {
                    U0[k][i] -= L0[j][i] * U0[k][j];
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
        ofstream outputf(output_filename);
        if (!outputf.is_open())
            abort_with_error_message("ERROR: Unable to open output file.");

        outputf << "L matrix:" << endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                outputf << L[i * N + j] << " ";
            }
            outputf << endl;
        }
        outputf << "U matrix:" << endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                outputf << U[i * N + j] << " ";
            }
            outputf << endl;
        }
        outputf.close();
    }

    // 将L0和U0写回L和U
    void writeback() {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                L[i * N + j] = L0[j][i];
                U[i * N + j] = U0[j][i];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3)
        utils::abort_with_error_message("ERROR: Invalid command. Usage: ./serial <input_file> <output_file>");

    string input_filename = argv[1];
    string output_filename = argv[2];

    utils::read_inputs_and_initialize(input_filename);
    utils::init_L0_U0();

    utils::lu_decomposition_0();
    utils::writeback();
    utils::write_outputs(output_filename);

    return 0;
}