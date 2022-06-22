#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <sstream>
#include <time.h>

#include "VCFParsingSorter.h"

int main(int argc, char **argv)
{
    // TODO: Larger version
    // main.cpp ../VCF_files/ ../VCF_files/reference.fa -n 1 ../VCF_files/test_4.vcf
    // For now
    if (argc < 6)
    {
        std::cout << "Usage: " << argv[0] << " destination_folder -n [NUMBER] [FILES] [OPTIONS]" << std::endl;
        std::cout << "       where [NUMBER] is the number of [FILES] to be given." << std::endl;
        return -1;
    }
    // INICIO PROCESO DE PARSING VIA PYTHON
    std::stringstream aux;
    aux << argv[4];
    int files_expected;
    aux >> files_expected;

    std::vector<std::string> py_params;

    py_params.push_back("python3 ./VCF_parsing/parsing_process.py");
    py_params.push_back(" ");
    py_params.push_back(argv[1]); // Destination folder
    py_params.push_back(" ");
    py_params.push_back(argv[2]); // reference file
    py_params.push_back(" "); 
    py_params.push_back(argv[3]); // -n
    py_params.push_back(" "); 
    py_params.push_back(argv[4]); // [NUMBER]

    for (int i = 5; i < 5 + files_expected; i++)
    {
        py_params.push_back(" ");
        py_params.push_back(argv[i]); // [FILES]
    }
    // TODO: Remain to include parsing options

    // Example: python3 parsing_process.py ../VCF_files/ -n 1 ../VCF_files/test_4.vcf
    std::string command = std::accumulate(py_params.begin(), py_params.end(), std::string(""));

    time_t begin, end;
    std::cout << "[RLZ] Start parsing process..." << std::endl;
    begin = clock();
    if (system(command.c_str()) == 0)
    {
        end = clock();
        std::cout << "[RLZ]\tParsing sucessful! Time elapsed: " << (float)(end - begin)/CLOCKS_PER_SEC << std::endl;
    } else
    {
        std::cout << "[RLZ]\tParsing failed, check log for more information :(" << std::endl;
        return -1;
    }

    VCFParsingSorter Sorter;
    // FIN PROCESO DE PARSING VIA PYTHON
    int result = Sorter.StartProcess(argv);
   

    return result;
}

// vim: et:ts=4:sw=4
