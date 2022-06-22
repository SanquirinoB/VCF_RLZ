#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>

using namespace std;

class VCFParsingInterpreter{
private:

    string Reference_file_subpath = "Tmp/Parsing/Reference.tmprlz";
    string Phrases_file_subpath = "Tmp/Parsing/phrases.tmprlz";
    string MetaReference_file_subpath = "Tmp/Meta_data/Reference.metarlz";
    string MetaParsing_file_subpath = "Tmp/Meta_data/Meta_info.metarlz";

    int n_Phrases;
    int n_Samples;

public:

    VCFParsingInterpreter(char* destination_path);

    void buildFactorFromVCFParserPhrase(char * phrases_path, vector<pair<unsigned int, unsigned int> > &factors);
};

