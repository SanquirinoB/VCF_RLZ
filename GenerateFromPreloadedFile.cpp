#include <fstream>
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
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " source_folder" << endl;
        return -1;
    }
    VCFParsingInterpreter* Interpreter = new VCFParsingInterpreter();

    Interpreter->InitializeFromPreloadedFile(argv[1]);
    
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
            vector<pair<sampleID, unsigned int>> result = Interpreter->FindSnippet(response, true);
        }
        
    }

    return 0;
}

// vim: et:ts=4:sw=4
