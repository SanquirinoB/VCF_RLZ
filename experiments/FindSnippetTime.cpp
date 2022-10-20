#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <functional>
#include <sstream>

// From RLZ
#include <relz/NanoTimer.h>
// #include <relz/RelzIndexReference.h>

#include "VCFParsingSorter.h"
#include "VCFParsingInterpreter.h"

using namespace std;

int main(int argc, char **argv)
{
    // string patterns_path_list(argv[1]);

    // ifstream patterns_file(patterns_path_list);
    // string aux;
    // vector<string> pattern_paths;

    // while(getline(patterns_file, aux))
    // {
    //     pattern_paths.push_back(aux);
    // }

    // patterns_file.close();

    VCFParsingInterpreter* Interpreter = new VCFParsingInterpreter();

    Interpreter->InitializeFromPreloadedFile(argv[2]);
    NanoTimer timer;
    string aux;
    ll occ = 0;

    ifstream patterns_file(argv[1]);

    cout << "Init pattern search of: " << argv[1] << endl;

    timer.reset();

    while(getline(patterns_file, aux))
    {
        occ += Interpreter->FindSnippetExperimental(aux);
    }

    double timeElapsed = timer.getMilisec();

    if(occ != 0)
    {
        cout << timeElapsed << " (ms)" << endl;
        cout << occ << " (occurences)" << endl;
        cout << (timeElapsed / occ) << " (real elapsed (ms))" << endl;
    }

    cout << "   Finished pattern search of: " << argv[1] << endl;
    

    return 0;
}

// vim: et:ts=4:sw=4
