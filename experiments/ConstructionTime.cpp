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
    string skipParse(argv[3]);
    string skipSort(argv[4]);
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
    ofstream results_file(results_file_path, ios_base::app);
    // results_file << "ID\tprocess\ttime\tfinalSize\n";
    results_file.flush();

    for (string command:commands)
    {
        string destPath = command.substr(0, command.find(" "));
        if(strcmp(skipParse.c_str(),"TRUE") != 0)
        {
            cout << "Executing: " << endl;
            cout << "\t" << ("python3 ./VCF_parsing/parsing_process.py " + command) << endl;
            cout << " PARSE starting..." << endl;
            timer.reset();
            // Python parsing
            if (system(("python3 ./VCF_parsing/parsing_process.py " + command).c_str()) != 0)
            {
                cout << "[RLZ]\tParsing failed, check log for more information :(" << endl;
                return -1;
            }
            timeElapsed = timer.getMilisec();

            results_file << destPath << "\t" << "PARSE" << "\t" << timeElapsed << "\t" << "?" << "\n";
            results_file.flush();
            cout << "   PARSE Ended: " << endl;
            cout << destPath << "\t" << "PARSE" << "\t" << timeElapsed << "\t" << "?" << "\n";
        }
        else
        {
            cout << "   PARSE skipped..." << endl;
        }

        if(strcmp(skipSort.c_str(),"TRUE") != 0)
        {
            cout << " SORT starting..." << endl;
            timer.reset();
            VCFParsingSorter* Sorter = new VCFParsingSorter();
            Sorter->StartProcess(destPath);
            timeElapsed = timer.getMilisec();
            results_file << destPath << "\t" << "SORT" << "\t" << timeElapsed << "\t" << "?" << "\n";
            results_file.flush();
            cout << "   SORT Ended: " << endl;
            cout << destPath << "\t" << "SORT" << "\t" << timeElapsed << "\t" << "?" << "\n";
            delete Sorter;
        }
        else
        {
            cout << "   SORT skipped..." << endl;
        }

        cout << " BUILD starting..." << endl;
        timer.reset();
        VCFParsingInterpreter* Interpreter = new VCFParsingInterpreter();
        Interpreter->InitializeFromParsing(destPath);

        timeElapsed = timer.getMilisec();
        size = Interpreter->GetSize();
        results_file << destPath << "\t" << "BUILD" << "\t" << timeElapsed << "\t" << size << "\n";
        results_file.flush();
        cout << "   BUILD Ended: " << endl;
        cout << destPath << "\t" << "BUILD" << "\t" << timeElapsed << "\t" << size << "\n";
        Interpreter->SaveInterpreter();
    }

    return 0;
}

// vim: et:ts=4:sw=4
