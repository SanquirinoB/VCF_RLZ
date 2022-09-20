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
    string commands_file_path(argv[1]);
    ifstream commands_file(commands_file_path);
    string command_aux;
    vector<string> commands;
    while(getline(commands_file, command_aux))
    {
        commands.push_back(command_aux);
    }

    NanoTimer timer;
    double timeElapsed, size;
    string results_file_path(argv[2]);
    ofstream results_file(results_file_path);
    results_file << "ID\ttime\tfinalSize\n";
    results_file.flush();

    for (string command:commands)
    {
        string destPath = command.substr(0, command.find(" "));
        cout << "Executing: " << endl;
        cout << "\t" << ("python3 ./VCF_parsing/parsing_process.py " + command) << endl;
        
        timer.reset();
        // Python parsing
        if (system(("python3 ./VCF_parsing/parsing_process.py " + command).c_str()) != 0)
        {
            cout << "[RLZ]\tParsing failed, check log for more information :(" << endl;
            return -1;
        }

        VCFParsingSorter* Sorter = new VCFParsingSorter();
        Sorter->StartProcess(destPath);

        VCFParsingInterpreter* Interpreter = new VCFParsingInterpreter();
        Interpreter->InitializeFromParsing(destPath);

        timeElapsed = timer.getMilisec();
        size = Interpreter->GetSize();
        results_file << destPath << "\t" << timeElapsed << "\t" << size << "\n";
        results_file.flush();
        cout << "   Ended: " << endl;
        cout << destPath << "\t" << timeElapsed << "\t" << size << endl;
        Interpreter->SaveInterpreter();
    }

    

    return 0;
}

// vim: et:ts=4:sw=4
