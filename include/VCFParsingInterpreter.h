#ifndef _VCF_PARSING_INTERPRETER_H
#define _VCF_PARSING_INTERPRETER_H

#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <relz/NanoTimer.h>
#include <relz/RelzIndexReference.h>
#include <VCFCommon.h>
#include <sdsl/bit_vectors.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

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
    int n_chromosomes;
    unsigned int ploidy = 2;
    vector_type Phrases;

    // Actual S related
    ll S_size;
    ll n_S_i;
    // Vector with checkpoints of sample position start inside S
    vector<ll> S_i_pos;
    vector<pair<unsigned int, unsigned int>> factors;
    ll ref_offset;
    ll rel_pos_chrom;

    // Reconstruction structures
    RelzIndexReference* Index;
    rrr_vector<127> bit_vector_S_i;
    rrr_vector<127>::rank_1_type rank_S_i;
    rrr_vector<127>::select_1_type select_S_i;

public:
    VCFParsingInterpreter(char *destination_path);
    vector<pair<sampleID, unsigned int>> FindSnippet(string snippet);
    void Initialize();

private:
    // Consume files
    void ProcessMetaReference();
    void ProcessMetaParsing();

    void ProcessReference();
    void ProcessPhrases();

    void BuildFactors();
    char* GetReference();
    string GetSampleRealName(unsigned int internal_id);

    void HigienicFactorPushBack(pair<unsigned int, unsigned int> factor);

    bool PhraseIsValidInit(phrase curr);
    bool NeedsToFullFill(phrase ref, phrase curr);
    bool ChangeInSameIndvChromAlele(phrase ref, phrase curr);
    bool ChangeInSameIndvChromDiffAlele(phrase ref, phrase curr);
    bool ChangeInSameIndvDiffChrom(phrase ref, phrase curr);
    bool ChangeInDiffIndv(phrase ref, phrase curr);

    pair<unsigned int, unsigned int> CalculateAleleInitFactor(phrase Phrase);
    pair<unsigned int, unsigned int> CalculateAleleInterFactor(phrase Phrase1, phrase Phrase2);
    pair<unsigned int, unsigned int> CalculateAleleEndFactor(phrase Phrase);
    pair<unsigned int, unsigned int> CalculateFullAleleFactor(ll chromosome);

    void BuildReconstructionStructures();

    void InduceFillFactors(phrase last_phrase, phrase curr_phrase, bool last_is_dummy, bool curr_is_dummy);

    void PrintState(phrase curr_phrase, ll last_pos, ll last_l, ll last_l_e);
};

#endif