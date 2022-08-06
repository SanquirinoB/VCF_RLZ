#ifndef _VCF_PARSING_INTERPRETER_H
#define _VCF_PARSING_INTERPRETER_H

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
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
    string IDInfo_file_subpath = "Tmp/Meta_data/ID_info.metarlz";

    string Reference_file_path;
    string Phrases_file_path;
    string MetaReference_file_path;
    string MetaParsing_file_path;
    string IDInfo_file_path;

    ifstream Reference_file;
    ifstream Phrases_file;
    ifstream MetaReference_file;
    ifstream MetaParsing_file;
    ifstream IDInfo_file;

    metareference metareference_buffer;
    metainfo metainfo_buffer;

    // Reference related
    int Reference_len;
    int n_References;
    map<int, metareference> dict_metareference{};
    map<int, sampleID> dict_samples{};

    // Parsing related
    ll n_Phrases;
    int n_Samples;
    vector_type Phrases;

public:
    VCFParsingInterpreter(char *destination_path);

    void Initialize();
    ll buildFactorFromVCFParserPhrase(vector<pair<unsigned int, unsigned int>> &factors);
    pair<const char *, ll> GetReference();

private:
    // Consume files
    void ProcessMetaReference();
    void ProcessMetaParsing();

    void ProcessReference();
    void ProcessPhrases();

    void UpdateSampleData(ll index, phrase data);
    void AddFactor(vector<pair<unsigned int, unsigned int>> &factors, ll pos, ll len);
};

#endif