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

    cout << "Structure size: " << Interpreter->GetSize()/ (1024 * 1024) << " (MB)" << endl;

    Interpreter->Index->querytime_p1 = 0;
	Interpreter->Index->querytime_p2 = 0;
	Interpreter->Index->querytime_p3 = 0;
	Interpreter->Index->querytime_p4 = 0;
	Interpreter->Index->occs_a = 0;
	Interpreter->Index->occs_b = 0;
	Interpreter->Index->occs_c = 0;
    NanoTimer timer;
    double q_p1 = 0, q_p2 = 0, q_p3 = 0, q_p4 = 0;
    string aux;
    ll occ = 0, occ_a = 0, occ_b = 0, occ_c = 0;

    ifstream patterns_file(argv[1]);

    cout << "Init pattern search of: " << argv[1] << endl;


    while(getline(patterns_file, aux))
    {
        occ += Interpreter->FindSnippetExperimental(aux);
        // Recover values
        q_p1 += Interpreter->Index->querytime_p1;
        q_p2 += Interpreter->Index->querytime_p2;
        q_p3 += Interpreter->Index->querytime_p3;
        q_p4 += Interpreter->Index->querytime_p4;
        occ_a += Interpreter->Index->occs_a;
        occ_b += Interpreter->Index->occs_b;
        occ_c += Interpreter->Index->occs_c;

        Interpreter->Index->querytime_p1 = 0;
        Interpreter->Index->querytime_p2 = 0;
        Interpreter->Index->querytime_p3 = 0;
        Interpreter->Index->querytime_p4 = 0;
        Interpreter->Index->occs_a = 0;
        Interpreter->Index->occs_b = 0;
        Interpreter->Index->occs_c = 0;

    }


    if(occ != 0)
    {
        cout << "[TT] " << q_p4 << endl;
        cout << "[T1] " << q_p1 << endl;
        cout << "[T2] " << q_p2 << endl;
        cout << "[T3] " << q_p3 << endl;
        cout << "[OT] " << occ << endl;
        cout << "[O1] " << occ_a << endl;
        cout << "[O2] " << occ_b << endl;
        cout << "[O3] " << occ_c << endl;
    }

    cout << "   Finished pattern search of: " << argv[1] << endl;
    

    return 0;
}

// vim: et:ts=4:sw=4
