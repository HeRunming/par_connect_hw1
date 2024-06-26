/*
 * Foundations of Parallel Computing II, Spring 2024.
 * Instructor: Chao Yang, Xiuhong Li @ Peking University.
 * This is a serial implementation of LU decomposition.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <omp.h>

using namespace std;

namespace utils {
    int N; // number of rows (and columns) in the matrix
    vector<double> mat; // the input matrix
    vector<double> L; // lower triangular matrix
    vector<double> U; // upper triangular matrix

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

        copy(mat.begin(), mat.end(), U.begin());

        // Initialize L to identity matrix
        for (int i = 0; i < N; ++i) {
            L[i * N + i] = 1.0;
        }
        return 0;
    }

    void lu_decomposition() {
        // 计时
        int num_threads = omp_get_max_threads();
        auto start = chrono::high_resolution_clock::now();
        #pragma omp parallel
        {
        int tid = omp_get_thread_num();
        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];}
        
            for (int i = j + 1 + tid; i < N; i += num_threads) {
                for (int k = j; k < N; ++k) {
                    U[i * N + k] -= L[i * N + j] * U[j * N + k];
                }
            }
        }
        }
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
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
}

int main(int argc, char* argv[]) {
    if (argc != 4)
        utils::abort_with_error_message("ERROR: Invalid command. Usage: ./serial <input_file> <output_file>");

    string input_filename = argv[1];
    string output_filename = argv[2];
    string num_threads_str = argv[3];
    int num_threads = stoi(num_threads_str);
    utils::read_inputs_and_initialize(input_filename);

    omp_set_num_threads(num_threads);
    utils::lu_decomposition();

    utils::write_outputs(output_filename);

    return 0;
}