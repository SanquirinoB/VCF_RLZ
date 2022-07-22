#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

#define BINARY_FILE ifstream::in | ifstream::binary

VCFParsingInterpreter::VCFParsingInterpreter(char *destination_path, vector_type &sorted_phrases)
{
    // Something like .../VCF_files/
    string Destination_path(destination_path);

    // TODO: Assert if ./Tmp is a valid folder
    cout << "1" << endl;
    Reference_file_path = Destination_path + Reference_file_subpath;
    Phrases_file_path = Destination_path + Phrases_file_subpath;
    MetaReference_file_path = Destination_path + MetaReference_file_subpath;
    MetaParsing_file_path = Destination_path + MetaParsing_file_subpath;
    IDInfo_file_path = Destination_path + IDInfo_file_subpath;

    Phrases = sorted_phrases;
    // Open all required files
    // OpenBinaryInput(Reference_file, Reference_file_path);
    // cout << "2" << endl;
    // OpenBinaryInput(Phrases_file, Phrases_file_path);
    // OpenBinaryInput(MetaReference_file, MetaReference_file_path);
    // cout << "3" << endl;
    // OpenBinaryInput(MetaParsing_file, MetaParsing_file_path);
}

void VCFParsingInterpreter::Initialize()
{
    ProcessReference();
    ProcessMetaParsing();
}

void VCFParsingInterpreter::ProcessReference()
{
    // Open required files
    MetaReference_file.open(MetaReference_file_path, BINARY_FILE);

    // The first structure contains summary info
    // cout << "Lee la estructura" << endl;
    MetaReference_file.read((char *)&metareference_buffer, sizeof(metareference));
    Reference_len = metareference_buffer.n_bases();
    n_References = metareference_buffer.rel_pos();

    // cout << "Largo: " << Reference_len << "|Referencias: " << n_References << endl;

    // Collect all metadata
    for (int i = 0; i < n_References; i++)
    {
        MetaReference_file.read((char *)&metareference_buffer, sizeof(metareference));
        // cout << i << endl;
        // cout << "\tReference: " << metareference_buffer.ID() << "|" << metareference_buffer.n_bases() << "|" << metareference_buffer.rel_pos() << endl;
        // char ref_content[metareference_buffer.n_bases() + 1];
        // ref_content[metareference_buffer.n_bases()] = '\0';
        // Reference_file.seekg(metareference_buffer.rel_pos());
        // Reference_file.readsome((char *)ref_content, metareference_buffer.n_bases());
        // cout << "\tContains: " << ref_content << endl;
        dict_metareference[metareference_buffer.ID()] = metareference_buffer;
    }
    MetaReference_file.close();
}

void VCFParsingInterpreter::ProcessMetaParsing()
{
    MetaParsing_file.open(MetaParsing_file_path, BINARY_FILE);
    MetaParsing_file.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_Phrases = metainfo_buffer.n_phrases();
    MetaParsing_file.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_Samples = metainfo_buffer.n_phrases();
    MetaParsing_file.close();
}

void VCFParsingInterpreter::buildFactorFromVCFParserPhrase(vector<pair<ll, ll>> &factors)
{
    ll new_factors = 0;
    // ID variables
    ll curr_indv = -1, curr_chrom = -1, curr_alele = -1;
    // Factor variables
    ll last_pos = 0, last_l = 0, last_l_e = 0;
    // Final text variable
    ll full_size = 0, full_size_snapshot = 0;
    ll n_append = -1;
    // Vector which indicates the end of each sample/chrom/alele in the final S
    vector<ll> end_pos;

    for (ll i = 0; i < n_Phrases; i++)
    {
        phrase_buffer = Phrases[i];

        ll delta_pos = phrase_buffer.pos() - last_pos;
        cout << "DeltaPos: " << delta_pos << endl;
        // Check Indv|Chrom|Alele as ID
        if (curr_alele != phrase_buffer.alele() ||
            curr_chrom != phrase_buffer.chrom() ||
            curr_indv != phrase_buffer.indv())
        {
            curr_alele = phrase_buffer.alele();
            curr_chrom = phrase_buffer.chrom();
            curr_indv = phrase_buffer.indv();
            end_pos.push_back(full_size);
            n_append += 1;
            full_size_snapshot = full_size;
        }

        if (delta_pos > 1) // delta = 0 or delta = 1, assert delta > 0?
        {

            // Entre cambios de ID no se debe hacer delta_pos
            ll inter_factor_len = phrase_buffer.pos() - last_pos + last_l;
            pair<ll, ll> inter_factor = make_pair(last_pos + last_l + full_size,
                                                  inter_factor_len);
            factors.push_back(inter_factor);
            full_size += inter_factor_len;
            new_factors += 1;
            cout << " relleno (" << last_pos + last_l << "," << inter_factor_len << ")" << endl;
        }
        else if (delta_pos == 1)
        {
            full_size += 1;
        }

        // Create current factor
        pair<ll, ll> curr_factor = make_pair(phrase_buffer.pos_e(), phrase_buffer.len_e());
        factors.push_back(curr_factor);
        new_factors += 1;

        cout << "(" << phrase_buffer.pos_e() << "," << phrase_buffer.len_e() << ")" << endl;

        full_size += phrase_buffer.len_e();

        // Update values for next iteration
        last_pos = phrase_buffer.pos();
        last_l = phrase_buffer.len();
        last_l_e = phrase_buffer.len_e();
    }

    Phrases_file.close();
    cout << new_factors << endl;
    return;
}