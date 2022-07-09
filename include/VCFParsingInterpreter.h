#ifndef _VCF_PARSING_INTERPRETER_H
#define _VCF_PARSING_INTERPRETER_H

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include <relz/NanoTimer.h>
#include <VCFCommon.h>

using namespace std;

class VCFParsingInterpreter
{
private:
    string Reference_file_subpath = "Tmp/Parsing/Reference.tmprlz";
    string Phrases_file_subpath = "Tmp/Parsing/phrases.tmprlz";
    string MetaReference_file_subpath = "Tmp/Meta_data/Reference.metarlz";
    string MetaParsing_file_subpath = "Tmp/Meta_data/Meta_info.metarlz";

    string Reference_file_path;
    string Phrases_file_path;
    string MetaReference_file_path;
    string MetaParsing_file_path;

    ifstream Reference_file;
    ifstream Phrases_file;
    ifstream MetaReference_file;
    ifstream MetaParsing_file;

    metareference metareference_buffer;
    metainfo metainfo_buffer;
    phrase phrase_buffer;

    int Reference_len;
    int n_References;
    int n_Phrases;
    int n_Samples;

public:
    VCFParsingInterpreter(char *destination_path);

    void StartProcess();

private:
    void buildFactorFromVCFParserPhrase(char *phrases_path, vector<pair<unsigned int, unsigned int>> &factors);

    // File related
    void OpenBinaryInput(ifstream *stream, string path);

    // Consume files
    void ProcessReference();
    void ProcessMetaReference();

    void ProcessPhrases();
    void ProcessMetaParsing();
};

#endif