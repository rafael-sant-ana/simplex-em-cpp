#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

const double TOLERANCE = 1e-9;

/**
 * FAST IO
 */

static inline const char* skip_ws(const char* p, const char* e) {
    while (p < e && (unsigned char)*p <= ' ') ++p;
    return p;
}

static inline bool next_int(const char*& p, const char* e, int& out) {
    p = skip_ws(p, e);
    if (p >= e) return false;
    const char* s = p;
    long v = 0;
    while (p < e && *p >= '0' && *p <= '9') { v = v*10 + (*p - '0'); ++p; }
    if (p == s) return false;
    out = v;
    return true;
}

static inline bool next_double(const char*& p, const char* e, double& out) {
    p = skip_ws(p, e);
    if (p >= e) return false;
    char* endp = nullptr;
    out = std::strtod(p, &endp);
    if (endp == p) return false;
    p = endp;
    return true;
}

// otimizacao: tableau contiguo

double** allocate_tableau_contiguous(int rows, int cols) {
    double** tableau = new double*[rows];
    double* data = new double[rows * cols]();
    for (int i = 0; i < rows; ++i) {
        tableau[i] = data + i * cols;
    }
    return tableau;
}

void free_tableau_contiguous(double** tableau) {
    if (tableau) {
        delete[] tableau[0]; 
        delete[] tableau;   
    }
}

int find_entering_col(double** tableau, int m_vars) {
    int entering_col = -1;
    double max_positive_val = TOLERANCE;
    
    for (int j = 0; j < m_vars; ++j) {
        double val = tableau[0][j];
        if (val > max_positive_val) {
            max_positive_val = val;
            entering_col = j;
        }
    }
    return entering_col;
}

int find_leaving_row(double** tableau, int n_rows, int n_cols, int pivot_col) {
    int leaving_row = -1;
    double min_ratio = std::numeric_limits<double>::max();
    int rhs_col = n_cols - 1;
    
    for (int i = 1; i < n_rows; ++i) {
        double coeff = tableau[i][pivot_col];
        if (coeff > TOLERANCE) {
            double ratio = tableau[i][rhs_col] / coeff;
            if (ratio < min_ratio - TOLERANCE) {
                min_ratio = ratio;
                leaving_row = i;
            }
        }
    }
    return leaving_row;
}

void pivot(double** tableau, int n_rows, int n_cols, int pivot_row, int pivot_col) {
    double inv_pivot = 1.0 / tableau[pivot_row][pivot_col];
    
    double* piv_row = tableau[pivot_row];
    for (int j = 0; j < n_cols; ++j) {
        piv_row[j] *= inv_pivot;
    }
    
    for (int i = 0; i < n_rows; ++i) {
        if (i == pivot_row) continue;
        
        double factor = tableau[i][pivot_col];
        if (std::abs(factor) < TOLERANCE) continue;
        
        double* curr_row = tableau[i];
        curr_row[pivot_col] = 0.0;
        
        for (int j = 0; j < n_cols; ++j) {
            if (j != pivot_col) {
                curr_row[j] -= factor * piv_row[j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "uso: " << argv[0] << " <arquivo_entrada> <arquivo_saida>" << std::endl;
        return 1;
    }

    // maluquices de fast io
    FILE* file = std::fopen(argv[1], "rb");
    if (!file) {
        std::cerr << "erro: arquivo de entrada nao encontrado" << std::endl;
        return 1;
    }

    std::fseek(file, 0, SEEK_END);
    long buffer_size = std::ftell(file);
    std::fseek(file, 0, SEEK_SET);

    char* buffer = (char*)std::malloc(buffer_size + 1);
    std::fread(buffer, 1, buffer_size, file);
    buffer[buffer_size] = '\0';
    std::fclose(file);

    const char* p = buffer;
    const char* e = buffer + buffer_size;

    int n, m;
    if (!next_int(p, e, n) || !next_int(p, e, m)) {
        std::free(buffer);
        return 1;
    }

    double* c = new double[m];
    double* b = new double[n];
    double** A = allocate_tableau_contiguous(n, m);

    // Lê dados
    for (int j = 0; j < m; ++j) {
        next_double(p, e, c[j]);
    }
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            next_double(p, e, A[i][j]);
        }
        next_double(p, e, b[i]);
    }

    std::free(buffer);

    // faz b >= 0 sempre
    for (int i = 0; i < n; ++i) {
        if (b[i] < 0) {
            b[i] = -b[i];
            for (int j = 0; j < m; ++j) {
                A[i][j] = -A[i][j];
            }
        }
    }

    // ========== FASE I ==========
    int n_rows_aux = n + 1;
    int n_cols_aux = m + n + 1;
    int m_aux = m + n;

    double** aux_tableau = allocate_tableau_contiguous(n_rows_aux, n_cols_aux);

    for (int i = 1; i < n_rows_aux; ++i) {
        int i_data = i - 1;
        std::memcpy(aux_tableau[i], A[i_data], m * sizeof(double));
        aux_tableau[i][m + i_data] = 1.0; 
        aux_tableau[i][n_cols_aux - 1] = b[i_data];
    }

    for (int j = m; j < m_aux; ++j) {
        aux_tableau[0][j] = -1.0;
    }

    for (int i = 1; i < n_rows_aux; ++i) {
        for (int j = 0; j < n_cols_aux; ++j) {
            aux_tableau[0][j] += aux_tableau[i][j];
        }
    }

    // Simplex Fase I
    while (true) {
        int pivot_col = find_entering_col(aux_tableau, m_aux);
        if (pivot_col == -1) break;

        int pivot_row = find_leaving_row(aux_tableau, n_rows_aux, n_cols_aux, pivot_col);
        if (pivot_row == -1) break;

        pivot(aux_tableau, n_rows_aux, n_cols_aux, pivot_row, pivot_col);
    }

    double fo_aux = aux_tableau[0][n_cols_aux - 1];
    
    std::ofstream outputFile(argv[2]);
    
    if (fo_aux > TOLERANCE) {
        outputFile << "inviavel" << std::endl;
        
        delete[] c;
        delete[] b;
        free_tableau_contiguous(A);
        free_tableau_contiguous(aux_tableau);
        outputFile.close();
        return 0;
    }

    // ========== FASE II ==========
    double** tableau = allocate_tableau_contiguous(n + 1, m + 1);

    for (int i = 1; i <= n; ++i) {
        std::memcpy(tableau[i], aux_tableau[i], m * sizeof(double));
        tableau[i][m] = aux_tableau[i][n_cols_aux - 1];
    }

    int* basis = new int[n];
    for (int i = 1; i <= n; ++i) {
        basis[i-1] = -1;
        for (int j = 0; j < m; ++j) {
            if (std::abs(aux_tableau[i][j] - 1.0) < TOLERANCE) {
                bool is_basic = true;
                for (int k = 1; k <= n; ++k) {
                    if (i != k && std::abs(aux_tableau[k][j]) > TOLERANCE) {
                        is_basic = false;
                        break;
                    }
                }
                if (is_basic) {
                    basis[i-1] = j;
                    break;
                }
            }
        }
    }

    free_tableau_contiguous(aux_tableau);

    std::memcpy(tableau[0], c, m * sizeof(double));
    
    for (int i = 1; i <= n; ++i) {
        int basic_col = basis[i-1];
        if (basic_col != -1) {
            double factor = tableau[0][basic_col];
            if (std::abs(factor) > TOLERANCE) {
                for (int col = 0; col <= m; ++col) {
                    tableau[0][col] -= factor * tableau[i][col];
                }
            }
        }
    }

    delete[] basis;

    // Simplex Fase II
    bool is_unbounded = false;
    
    while (true) {
        int pivot_col = find_entering_col(tableau, m);
        if (pivot_col == -1) break;

        int pivot_row = find_leaving_row(tableau, n + 1, m + 1, pivot_col);
        if (pivot_row == -1) {
            is_unbounded = true;
            break;
        }

        pivot(tableau, n + 1, m + 1, pivot_row, pivot_col);
    }

    if (is_unbounded) {
        outputFile << "ilimitada" << std::endl;
    } else {
        double optimal_value = -tableau[0][m];
        outputFile << "otima" << std::endl;
        outputFile << std::fixed << std::setprecision(3) << optimal_value << std::endl;
    }

    delete[] c;
    delete[] b;
    free_tableau_contiguous(A);
    free_tableau_contiguous(tableau);
    outputFile.close();

    return 0;
}