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
    // 初始化 L0 和 U0
    vector<vector<double>> L0(N, vector<double>(N, 0));
    vector<vector<double>> U0(N, vector<double>(N, 0));
    int init_L0_U0() {
        L0.resize(N, vector<double>(N, 0));
        U0.resize(N, vector<double>(N, 0));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                L0[j][i] = L[i * N + j];
                U0[j][i] = U[i * N + j];
            }
        }
        return 0;
    }
        



    void lu_decomposition() {
        // 开始计时
        auto start = chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            #pragma omp parallel for schedule(dynamic)
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];
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

    void lu_decomposition_0() {
        // 开始计时
        auto start = chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            #pragma omp parallel for schedule(dynamic)
            for (int i = j + 1; i < N; ++i) {L0[j][i] = U0[j][i] / U0[j][j];}
            #pragma omp parallel for schedule(static, 64)
            for (int i = j + 1; i < N; ++i) {
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

void lu_decomposition_lab() {
        // 获取线程数
        int num_threads = omp_get_num_threads();
        // 开始计时
        auto start = chrono::high_resolution_clock::now();
        # pragma omp parallel
        {
        for (int j = 0; j < N; ++j) {
            //# pragma omp parallel for schedule(dynamic)
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];
            }
            # pragma omp for collapse(2) schedule(static, 8 * num_threads)
            for (int i = j + 1; i < N; ++i) {
                for (int k = j; k < N; ++k) {
                    U[i * N + k] -= L[i * N + j] * U[j * N + k];
                }
            }
        }
        }

        // 结束计时
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

        // 输出并行计算所用的时间
        cout << "Parallel LU decomposition took " << duration.count() << " microseconds." << endl;
}


void lu_decomposition_lab_plus() {
        // 开始计时
        auto start = chrono::high_resolution_clock::now();
        for (int j = 0; j < N; ++j) {
            //# pragma omp parallel for schedule(dynamic)
            for (int i = j + 1; i < N; ++i) {
                L[i * N + j] = U[i * N + j] / U[j * N + j];
            }
        int num_threads = 0; // 线程数
        // 细致的分配所有的任务，并行度更高
        # pragma omp parallel
        {
            int tid = omp_get_thread_num(); // 线程号
            if (tid == 0) num_threads = omp_get_num_threads(); // 获得线程数
            // 总任务量为(N - j - 1) * (N - j)
            for (int i = j + 1 + tid; i < N; i += num_threads) {
                for (int k = j; k < N; ++k) {
                    U[i * N + k] -= L[i * N + j] * U[j * N + k];
                }
            }
        }
        }

        // 结束计时
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

        // 输出并行计算所用的时间
        cout << "Parallel LU decomposition took " << duration.count() << " microseconds." << endl;
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
    if (argc != 3)
        utils::abort_with_error_message("ERROR: Invalid command. Usage: ./parallel_lu <input_file> <output_file>");

    string input_filename = argv[1];
    string output_filename = argv[2];


    int nlst[5] = {1,2,4,8,16};
    string slst[5] = {"1","2","4","8","16"};
    int num = 0;
    for (num = 0;num < 5;num++){
        utils::U.clear();
        utils::L.clear();
        utils::read_inputs_and_initialize(input_filename);
        //utils::init_L0_U0();
        int num_threads = nlst[num];
        omp_set_num_threads(num_threads);
        cout << "num_threads = " << num_threads << "  ";

        utils::lu_decomposition_lab();
        //utils::writeback();
        string out_name = output_filename + slst[num] + string(".txt");
        utils::write_outputs(out_name);
    }

    return 0;
}
