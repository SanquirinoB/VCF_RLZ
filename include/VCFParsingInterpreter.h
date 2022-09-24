#ifndef _VCF_PARSING_INTERPRETER_H
#define _VCF_PARSING_INTERPRETER_H

#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <NanoTimer.h>
#include <RelzIndexReference.h>
#include <VCFCommon.h>
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <ctime>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace sdsl;

class VCFParsingInterpreter
{
private:
    // Required to be reminded
    string Destination_path;
    string IDInfo_file_path;

    // For internal construction (clear not required)
    string Reference_file_subpath = "Tmp/Parsing/Reference.tmprlz";
    string Phrases_file_subpath = "Tmp/Parsing/phrases.tmprlz";
    string MetaReference_file_subpath = "Tmp/Meta_data/Reference.metarlz";
    string MetaParsing_file_subpath = "Tmp/Meta_data/Meta_info.metarlz";
    string IDInfo_file_subpath = "Tmp/Meta_data/ID_info.metarlz";
    ifstream Reader;
    ll ploidy = (ll) 2;

    // For internal construction (clear required)
    ll Reference_len;
    map<int, metareference> dict_metareference{};

    // Parsing related
    ll n_Phrases;
    ll n_chromosomes;
    vector_type Phrases;

    // Actual S related
    ll S_size;
    ll n_S_i;
    // Vector with checkpoints of sample position start inside S
    vector<ll> S_i_pos;
    vector<pair<ull, ull>> factors;
    ll ref_offset;
    ll rel_pos_chrom;

    // Reconstruction structures
    RelzIndexReference *Index;
    sd_vector<> bit_vector_S_i;
    sd_vector<>::rank_1_type rank_S_i;
    sd_vector<>::select_1_type select_S_i;

public:
    VCFParsingInterpreter();
    vector<pair<sampleID, ll>> FindSnippet(string snippet, bool show = false);
    void InitializeFromParsing(char *destination_path);
    void InitializeFromParsing(string destination_path);
    void InitializeFromPreloadedFile(char *folder_path);

    void SaveInterpreter();

    double GetSize();

private:
    // Consume files
    void ProcessMetaReference();
    void ProcessMetaParsing();

    void ProcessReference();
    void ProcessPhrases();

    void BuildFactors();
    string GetSampleRealName(ll internal_id);

    void HigienicFactorPushBack(pair<ll, ll> factor);

    bool PhraseIsValidInit(phrase curr);
    bool NeedsToFullFill(phrase ref, phrase curr);
    bool ChangeInSameIndvChromAlele(phrase ref, phrase curr);
    bool ChangeInSameIndvChromDiffAlele(phrase ref, phrase curr);
    bool ChangeInSameIndvDiffChrom(phrase ref, phrase curr);
    bool ChangeInDiffIndv(phrase ref, phrase curr);

    pair<ll, ll> CalculateAleleInitFactor(phrase Phrase);
    pair<ll, ll> CalculateAleleInterFactor(phrase Phrase1, phrase Phrase2);
    pair<ll, ll> CalculateAleleEndFactor(phrase Phrase);
    pair<ll, ll> CalculateFullAleleFactor(ll chromosome);

    void BuildReconstructionStructures();

    void InduceFillFactors(phrase last_phrase, phrase curr_phrase, bool last_is_dummy, bool curr_is_dummy);

    void PrintState(phrase curr_phrase, ll last_pos, ll last_l, ll last_l_e);

};

#endif