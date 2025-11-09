#include <iostream> // lidar com io
#include <fstream>  // lidar com arquivos
#include <string>   // std::string
#include <iomanip>  // std::setw, std::setprecision
#include <limits>   // std::numeric_limits
#include <cmath>    // std::abs

const double TOLERANCE = 1e-9;

/**
 * Encontra a coluna para entrar na base (coluna pivô).
 * Procura pelo valor mais positivo na linha da função objetivo (linha 0).
 */
int find_entering_col(double** tableau, int m_vars) {
    int entering_col = -1;
    double max_positive_val = TOLERANCE; 

    for (int j = 0; j < m_vars; ++j) {
        if (tableau[0][j] > max_positive_val) {
            max_positive_val = tableau[0][j];
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
        double rhs = tableau[i][rhs_col];

        // Teste da razao so e valido para coeficientes > 0 !
        if (coeff > TOLERANCE) {
            double ratio = rhs / coeff;
            if (ratio < min_ratio - TOLERANCE) {
                min_ratio = ratio;
                leaving_row = i;
            }
        }
    }
    return leaving_row; // -1 se ilimitado
}

/**
 * Executa a operacao de pivotamento no tableau.
 * mto chato slk
 */
void pivot(double** tableau, int n_rows, int n_cols, int pivot_row, int pivot_col) {
    double pivot_val = tableau[pivot_row][pivot_col];
    // 1. Normalizar a linha pivo
    for (int j = 0; j < n_cols; ++j) {
        tableau[pivot_row][j] /= pivot_val;
    }

    for (int i = 0; i < n_rows; ++i) {
        if (i != pivot_row) {
            double factor = tableau[i][pivot_col];
            if (std::abs(factor) < TOLERANCE) {
                continue;
            }
            // L_i = L_i - fator * L_pivo
            // igual gauss jordan
            for (int j = 0; j < n_cols; ++j) {
                tableau[i][j] -= factor * tableau[pivot_row][j];
            }
        }
    }
}

/*
 * Funcao helper para imprimir o tableau no arquivo de saida.
 * comentei os usos na entrega final. mas era assim q eu tava fznd pra ver
 */
void debugPrintTableau(std::ofstream& outFile, const std::string& titulo,
                       double** tableau, int n_rows, int n_cols, int m) {
    
    outFile << "\n--- DEBUG: " << titulo << " ---\n";
    outFile << std::fixed << std::setprecision(3); 

    for (int i = 0; i < n_rows; ++i) {
        if (i == 1) {
            for(int k=0; k < n_cols + 2; ++k) outFile << "----------";
            outFile << "\n";
        }
        for (int j = 0; j < n_cols; ++j) {
            // Adiciona | para separar A, Artificiais, e b
            if ( (m > 0 && j == m) || j == n_cols - 1) { 
                 outFile << " | ";
            }
            outFile << std::setw(8) << tableau[i][j];
        }
        outFile << "\n";
    }
    outFile << "---------------------------------------\n";
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Uso: " << argv[0] << " <arquivo_entrada> <arquivo_saida>" << std::endl;
        return 1;
    }

    std::string nomeArquivoEntrada = argv[1];
    std::string nomeArquivoSaida = argv[2];

    std::ifstream inputFile(nomeArquivoEntrada);
    if (!inputFile.is_open()) {
        std::cerr << "erro: n deu pra abrir a entrada" << nomeArquivoSaida << "'" << std::endl;
        return 1;
    }

    std::ofstream outputFile(nomeArquivoSaida);
    if (!outputFile.is_open()) {
        std::cerr << "erro: n deu pra abrir a saida" << nomeArquivoSaida << "'" << std::endl;
        return 1;
    }

    int n, m;
    inputFile >> n >> m;

    if (inputFile.fail()) {
         std::cerr << "impossivel ler `n` e `m` na entrada" << std::endl;
         return 1;
    }

    double* c = nullptr;
    try { c = new double[m]; } 
    catch (const std::bad_alloc& e) { return 1; }
    
    double* b = nullptr;
     try { b = new double[n]; } 
     catch (const std::bad_alloc& e) { delete[] c; return 1; }

    double** A = nullptr;
    try {
        A = new double*[n]; 
        for (int i = 0; i < n; ++i) { A[i] = new double[m]; }
    } catch (const std::bad_alloc& e) {
        delete[] c; delete[] b;
        if (A) { 
             for (int i = 0; i < n; ++i) { if(A[i]) delete[] A[i]; }
             delete[] A;
        }
        return 1;
    }

    for (int j = 0; j < m; ++j) { inputFile >> c[j]; }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) { inputFile >> A[i][j]; }
        inputFile >> b[i];
    }
    if (inputFile.fail()) { std::cerr << "formato de arquivo invalido ou dados faltantes" << std::endl; } 
    // leu certin

    for (int i = 0; i < n; ++i) {
        if (b[i] < 0) {
            b[i] *= -1.0;
            for (int j = 0; j < m; ++j) { A[i][j] *= -1.0; }
        }
    }

    double** tableau = nullptr;
    try {
        tableau = new double*[n + 1];
        for (int i = 0; i < n + 1; ++i) { tableau[i] = new double[m + 1](); }
    } catch (const std::bad_alloc& e) {
        std::cerr << "falha ao alocar tableau Fase II." << std::endl;
        delete[] c; delete[] b;
        for (int i = 0; i < n; ++i) { delete[] A[i]; }
        delete[] A;
        return 1;
    }

    int n_rows_aux = n + 1;
    int n_cols_aux = m + n + 1;
    int m_aux = m + n; 

    double** aux_tableau = nullptr;
    try {
        aux_tableau = new double*[n_rows_aux];
        for (int i = 0; i < n_rows_aux; ++i) { aux_tableau[i] = new double[n_cols_aux](); }
    } catch (const std::bad_alloc& e) {
        std::cerr << "falha ao alocar aux_tableau Fase I." << std::endl;
        delete[] c; delete[] b;
        for (int i = 0; i < n; ++i) { delete[] A[i]; }
        delete[] A;
        for (int i = 0; i < n+1; ++i) { delete[] tableau[i]; }
        delete[] tableau;
        return 1;
    }

    // Popular aux_tableau [ A | I | b ] e F.O. [ 0 | -1 | 0 ] ---
    for (int i = 1; i < n_rows_aux; ++i) {
        int i_data = i - 1; 
        for (int j = 0; j < m; ++j) { aux_tableau[i][j] = A[i_data][j]; } // A
        aux_tableau[i][m + i_data] = 1.0; // I
        aux_tableau[i][n_cols_aux - 1] = b[i_data]; // b
    }
    // toda vez que eu falar f.o. eh funcao objetivo
    for (int j = m; j < m_aux; ++j) { aux_tableau[0][j] = -1.0; } // F.O.

    for (int i = 1; i < n_rows_aux; ++i) { 
        for (int j = 0; j < n_cols_aux; ++j) {
            aux_tableau[0][j] += aux_tableau[i][j];
        }
    }
    
    // debugPrintTableau(outputFile, "Tableau Fase I (Inicial)", aux_tableau, n_rows_aux, n_cols_aux, m);
    
    int iter = 0; 

    while (true) {
        int pivot_col = find_entering_col(aux_tableau, m_aux);
        if (pivot_col == -1) {
            break; 
        }

        int pivot_row = find_leaving_row(aux_tableau, n_rows_aux, n_cols_aux, pivot_col);
        if (pivot_row == -1) {
            break; 
        }

        pivot(aux_tableau, n_rows_aux, n_cols_aux, pivot_row, pivot_col);
    }
    // debugPrintTableau(outputFile, "Tableau Fase I (Final)", aux_tableau, n_rows_aux, n_cols_aux, m);

    double fo_aux = aux_tableau[0][n_cols_aux - 1];
    
    if (fo_aux > TOLERANCE) {
        // std::cout << "Resultado: Problema INVIAVEL (F.O. Fase I = " << fo_aux << ")" << std::endl;
        outputFile << "inviavel" << std::endl; 

        delete[] c; delete[] b;
        for (int i = 0; i < n; ++i) { delete[] A[i]; }
        delete[] A;
        if (tableau) {
            for (int i = 0; i < n + 1; ++i) { delete[] tableau[i]; }
            delete[] tableau;
        }
        if (aux_tableau) {
            for (int i = 0; i < n_rows_aux; ++i) { delete[] aux_tableau[i]; }
            delete[] aux_tableau;
        }
        inputFile.close();
        outputFile.close();
        // std::cout << "Processamento concluido." << std::endl;
        return 0;
    }

    // std::cout << "Fase I concluida: Solucao viavel encontrada." << std::endl;

    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j < m; ++j) {
            tableau[i][j] = aux_tableau[i][j]; // A'
        }
        tableau[i][m] = aux_tableau[i][n_cols_aux - 1]; // b'
    }
    for (int j = 0; j < m; ++j) {
        tableau[0][j] = c[j]; // -c
    }

    for (int i = 1; i <= n; ++i) { // Para cada linha
        int basic_col = -1;
        for (int j = 0; j < m; ++j) {
            if (std::abs(tableau[i][j] - 1.0) < TOLERANCE) {
                bool is_basic = true;
                for (int k = 1; k <= n; ++k) {
                    if (i != k && std::abs(tableau[k][j]) > TOLERANCE) {
                        is_basic = false;
                        break;
                    }
                }
                if (is_basic) {
                    basic_col = j;
                    break;
                }
            }
        }
        if (basic_col != -1) {
            // L0 = L0 - (custo_da_var_basica) * Li
            double factor = tableau[0][basic_col];
            if (std::abs(factor) < TOLERANCE) continue; 
            for (int col = 0; col <= m; ++col) {
                tableau[0][col] -= factor * tableau[i][col];
            }
        }
    }
    
    // debugPrintTableau(outputFile, "Tableau Fase II (Inicial, Apos Price-Out)", tableau, n + 1, m + 1, m);

    // std::cout << "Iniciando Fase II..." << std::endl;
    bool is_unbounded = false;
    iter = 0; 

    while (true) {
        int pivot_col = find_entering_col(tableau, m); 

        if (pivot_col == -1) {
            // std::cout << "Fase II: Solucao otima encontrada." << std::endl;
            break; 
        }

        int pivot_row = find_leaving_row(tableau, n + 1, m + 1, pivot_col);

        if (pivot_row == -1) {
            is_unbounded = true;
            break; 
        }
        
        pivot(tableau, n + 1, m + 1, pivot_row, pivot_col);
        
    }
    
    // debugPrintTableau(outputFile, "Tableau Fase II (Final)", tableau, n + 1, m + 1, m);

    if (is_unbounded) {
        // std::cout << "Resultado: Problema ILIMITADO." << std::endl;
        outputFile << "ilimitada" << std::endl; 
        
        delete[] c; delete[] b;
        for (int i = 0; i < n; ++i) { delete[] A[i]; }
        delete[] A;
        if (tableau) {
            for (int i = 0; i < n + 1; ++i) { delete[] tableau[i]; }
            delete[] tableau;
        }
        if (aux_tableau) {
            for (int i = 0; i < n_rows_aux; ++i) { delete[] aux_tableau[i]; }
            delete[] aux_tableau;
        }
        inputFile.close();
        outputFile.close();
        // std::cout << "Processamento concluido." << std::endl;
        return 0; 
    
    } else {
        // std::cout << "Resultado: Solucao OTIMA encontrada." << std::endl;
        double optimal_value = -1 * tableau[0][m];
        
        outputFile << "otima" << std::endl; 
        outputFile << std::fixed << std::setprecision(3) << optimal_value << std::endl; 
        
        // std::cout << "Valor Otimo: " << std::fixed << std::setprecision(3) << optimal_value << std::endl;

    }


    delete[] c;
    delete[] b;
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
    }
    delete[] A;

    if (tableau) {
        for (int i = 0; i < n + 1; ++i) {
            delete[] tableau[i];
        }
        delete[] tableau;
    }

    if (aux_tableau) {
        for (int i = 0; i < n_rows_aux; ++i) {
            delete[] aux_tableau[i];
        }
        delete[] aux_tableau;
    }

    inputFile.close();
    outputFile.close();

    return 0; 
}