#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <functional>
#include <sstream>

// From RLZ
#include <relz/NanoTimer.h>
#include <relz/RelzIndexReference.h>

#include "VCFParsingSorter.h"
#include "VCFParsingInterpreter.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 6)
    {
        cout << "Usage: " << argv[0] << " destination_folder -n [NUMBER] [FILES] [OPTIONS]" << endl;
        cout << "       where [NUMBER] is the number of [FILES] to be given." << endl;
        return -1;
    }

    // ###########################################
    //    INICIO PROCESO DE PARSING VIA PYTHON
    // ###########################################

    // Preparar comando
    stringstream aux;
    aux << argv[4];
    int files_expected;
    aux >> files_expected;

    vector<string> py_params;
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
    string command = accumulate(py_params.begin(), py_params.end(), string(""));

    cout << "[RLZ] Start parsing process..." << endl;
    NanoTimer timer;
    if (system(command.c_str()) == 0)
    {
        cout << "[RLZ]\tParsing sucessful! Time elapsed: " << timer.getMilisec() << "ms" << endl;
    }
    else
    {
        cout << "[RLZ]\tParsing failed, check log for more information :(" << endl;
        return -1;
    }

    VCFParsingSorter* Sorter = new VCFParsingSorter();
    cout << "[RLZ] Start sorting process..." << endl;
    // FIN PROCESO DE PARSING VIA PYTHON
    Sorter->StartProcess(argv);

    cout << "[RLZ] Sorting process finished!" << endl;

    cout << "[RLZ] Start indexing process..." << endl;
    VCFParsingInterpreter* Interpreter = new VCFParsingInterpreter();

    Interpreter->InitializeFromParsing(argv[1]);

    cout << "[RLZ] Indexing ended!" << endl;
    
    string response;
    while(response != "2")
    {
        cout << "Puedes ejecutar las siguientes acciones:\n\t(1) Buscar un snippet\n\t(2) Salir" << endl;
        cout << "Ingresa el nro de tu opción:" << endl;
        cin >> response;
        if (response == "1")
        {
            cout << "Ingrese el snippet a consultar:" << endl;
            cin >> response;
            vector<pair<sampleID, ll>> result = Interpreter->FindSnippet(response, true);
        }
        
    }
    Interpreter->SaveInterpreter();

    return 0;
}

// vim: et:ts=4:sw=4
