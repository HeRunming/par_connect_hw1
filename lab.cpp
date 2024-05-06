#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
int main() {
    int n = 100;
    // 二维vector，每一个元素代表一列，每一列是一个vector
    vector<vector<double>> matrix(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = i * j;
        }
    }
}