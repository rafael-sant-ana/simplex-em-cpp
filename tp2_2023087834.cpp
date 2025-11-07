#include <iostream> 
#include <fstream>  
#include <string>   

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Erro: Numero incorreto de argumentos." << std::endl;
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
    try {
        c = new double[m];
    } catch (const std::bad_alloc& e) {
        return 1;
    }
    
    double* b = nullptr;
     try {
        b = new double[n];
    } catch (const std::bad_alloc& e) {
        delete[] c; 
        return 1;
    }

    double** A = nullptr;
    try {
        A = new double*[n]; 
        for (int i = 0; i < n; ++i) {
            A[i] = new double[m]; 
        }
    } catch (const std::bad_alloc& e) {
        delete[] c;
        delete[] b;
        if (A) { 
             for (int i = 0; i < n; ++i) {
                if(A[i]) delete[] A[i];
             }
             delete[] A;
        }
        return 1;
    }


    // lendo dados
    
    for (int j = 0; j < m; ++j) {
        inputFile >> c[j];
    }

    for (int i = 0; i < n; ++i) {

      for (int j = 0; j < m; ++j) {
            inputFile >> A[i][j];
        }

        inputFile >> b[i];
    }

    if (inputFile.fail()) {
        std::cerr << "formato de arquivo invalido ou dados faltantes" << std::endl;
    } else {
        std::cout << "leitura concluida com sucesso" << std::endl;
    }

    // Liberar memoria dinamica

    delete[] c;
    
    delete[] b;

    for (int i = 0; i < n; ++i) {
        delete[] A[i];
    }
    delete[] A;

    inputFile.close();
    outputFile.close();

    std::cout << "Processamento concluido." << std::endl;
    std::cout << "Dados lidos de '" << nomeArquivoEntrada << "' e resultado escrito em '" << nomeArquivoSaida << "'." << std::endl;

    return 0; 
}