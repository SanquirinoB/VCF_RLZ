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
        cout << "Ref: " << metareference_buffer.ID() << ", n_bases: " << metareference_buffer.n_bases() << endl;
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

// void VCFParsingInterpreter::buildFactorFromVCFParserPhrase(vector<pair<ll, ll>> &factors)
// {
//     ll new_factors = 0;
//     // ID variables
//     ll curr_indv = -1, curr_chrom = -1, curr_alele = -1;
//     // Factor variables
//     ll last_pos = 0, last_l = 0, last_l_e = 0, internal_pos = 0;
//     // Final text variable
//     ll full_size = 0, full_size_snapshot = 0;
//     ll n_append = -1;
//     // Vector which indicates the end of each sample/chrom/alele in the final S
//     vector<ll> end_pos;

//     for (ll i = 0; i < n_Phrases; i++)
//     {
//         phrase_buffer = Phrases[i];
//         cout << "Leo la phrase " << i << endl;
//         cout << "\t" << phrase_buffer.indv() << "|" << phrase_buffer.chrom() << "|" << phrase_buffer.alele() << "|" << phrase_buffer.pos() << "|"
//              << phrase_buffer.len() << "|" << phrase_buffer.edit() << "|" << phrase_buffer.pos_e() << "|" << phrase_buffer.len_e() << endl;
//         cout << "\t Estado: last_pos|last_l|last_l_e|full_size|internal_pos" << endl;
//         cout << "\t Estado:" << last_pos << "|" << last_l << "|" << last_l_e << "|" << full_size << "|" << internal_pos << endl;

//         ll delta_pos = (phrase_buffer.pos() + internal_pos) - last_pos;
//         cout << "DeltaPos: " << delta_pos << endl;

//         // Check Indv|Chrom|Alele as ID
//         if (curr_alele != phrase_buffer.alele() ||
//             curr_chrom != phrase_buffer.chrom() ||
//             curr_indv != phrase_buffer.indv())
//         {

//             // Reset position variables
//             if (last_pos != 0)
//             {
//                 // Tenemos que hacer un relleno entre el ultimo edit y el final del cromosoma asociado
//                 ll end_factor_len = dict_metareference[curr_chrom].n_bases() - (last_pos - internal_pos);
//                 pair<ll, ll> inter_factor = make_pair(internal_pos + 1,
//                                                       end_factor_len);
//                 factors.push_back(inter_factor);
//                 full_size += end_factor_len;
//                 internal_pos += end_factor_len;
//                 new_factors += 1;
//                 cout << " End: (" << inter_factor.first << "," << inter_factor.second << ")" << endl;
//             }

//             last_pos = 0;
//             last_l = 0;
//             last_l_e = 0;
//             // Update ID variables
//             curr_alele = phrase_buffer.alele();
//             curr_chrom = phrase_buffer.chrom();
//             curr_indv = phrase_buffer.indv();
//             // Report position of change
//             end_pos.push_back(full_size);
//             // Update length variables
//             n_append += 1;
//             full_size_snapshot = full_size;
//         }

//         if (delta_pos > 1) // delta = 0 or delta = 1, assert delta > 0?
//         {

//             // Entre cambios de ID no se debe hacer delta_pos
//             ll inter_factor_len = (phrase_buffer.pos() + internal_pos) - last_pos + last_l;
//             pair<ll, ll> inter_factor = make_pair(last_pos + last_l + internal_pos,
//                                                   inter_factor_len);
//             factors.push_back(inter_factor);
//             full_size += inter_factor_len;
//             internal_pos += inter_factor_len;
//             new_factors += 1;
//             cout << " Inter: (" << inter_factor.first << "," << inter_factor.second << ")" << endl;
//         }
//         else if (delta_pos == 1)
//         {
//             full_size += 1;
//         }

//         // Create current factor
//         pair<ll, ll> curr_factor = make_pair(phrase_buffer.pos_e(), phrase_buffer.len_e());
//         factors.push_back(curr_factor);
//         new_factors += 1;

//         cout << "(" << phrase_buffer.pos_e() << "," << phrase_buffer.len_e() << ")" << endl;

//         full_size += phrase_buffer.len_e();

//         // Update values for next iteration
//         last_pos = phrase_buffer.pos() + internal_pos;
//         last_l = phrase_buffer.len();
//         last_l_e = phrase_buffer.len_e();
//     }

//     Phrases_file.close();
//     cout << new_factors << endl;
//     return;
// }

void VCFParsingInterpreter::buildFactorFromVCFParserPhrase(vector<pair<ll, ll>> &factors)
{
    // Position variables: Offset checkpoints
    ll last_pos = 0, last_l = 0, last_l_e = 0, pos_ref_consumed = 0;
    // Counter variables: S info
    ll n_S_i = 0;
    vector<ll> S_i_pos;

    // Helper variables: S construction
    ll S_ref_actual_pos = 0, S_size = 0, factor_pos = 0, factor_len = 0;
    // Helper variables: Phrase ID
    ll last_indv = -1, last_chrom = -1, last_alele = -1;

    // Initialization
    phrase_buffer = Phrases[0];
    cout << phrase_buffer.indv() << "|" << phrase_buffer.chrom() << "|" << phrase_buffer.alele() << "|" << phrase_buffer.pos() << "|"
         << phrase_buffer.len() << "|" << phrase_buffer.edit() << "|" << phrase_buffer.pos_e() << "|" << phrase_buffer.len_e() << endl;
    cout << "State: " << last_pos << "," << last_l << "," << last_l_e << endl;

    last_indv = phrase_buffer.indv();
    last_chrom = phrase_buffer.chrom();
    last_alele = phrase_buffer.alele();

    if (phrase_buffer.pos() != 0)
    {
        // If the first phrase edit doesnt occurr in the first base
        // we need to induce a completition factor
        factor_len = phrase_buffer.pos() + 1;
        factors.push_back(make_pair(factor_pos, factor_len));
        cout << " Init: (" << factor_pos << "," << factor_len << ")->" << endl;
        S_ref_actual_pos += factor_len;
        pos_ref_consumed += factor_len;
        S_size += factor_len;
    }
    // And now create the factor associated to the current phrase
    factors.push_back(make_pair(phrase_buffer.pos_e(), phrase_buffer.len_e()));
    cout << " Actu: (" << phrase_buffer.pos_e() << "," << phrase_buffer.len_e() << ")" << endl;
    S_size += phrase_buffer.len_e();

    // Update with current values
    last_pos = phrase_buffer.pos();
    last_l = phrase_buffer.len();
    last_l_e = phrase_buffer.len_e();

    for (ll i = 1; i < n_Phrases; i++)
    {
        // For each phrase
        phrase_buffer = Phrases[i];
        cout << phrase_buffer.indv() << "|" << phrase_buffer.chrom() << "|" << phrase_buffer.alele() << "|" << phrase_buffer.pos() << "|"
             << phrase_buffer.len() << "|" << phrase_buffer.edit() << "|" << phrase_buffer.pos_e() << "|" << phrase_buffer.len_e() << endl;
        cout << "State: " << last_pos << "," << last_l << "," << last_l_e << endl;

        // If the current sample changes
        if (last_indv != phrase_buffer.indv() ||
            last_chrom != phrase_buffer.chrom() ||
            last_alele != phrase_buffer.alele())
        {
            // Create an end completition factor for the last sample
            factor_pos = S_ref_actual_pos;
            factor_len = dict_metareference[last_chrom].n_bases() - (pos_ref_consumed + 1);
            factors.push_back(make_pair(factor_pos, factor_len));
            cout << " ECmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
            cout << "\tState: " << S_ref_actual_pos << endl;

            S_ref_actual_pos += factor_len;
            S_size += factor_len;

            // Update positon checkpoint variables
            last_pos = 0;
            last_l = 0;
            last_l_e = 0;
            pos_ref_consumed = 0;

            if (last_indv != phrase_buffer.indv())
            {
                S_ref_actual_pos = 0;
            }

            // Update ID variables
            last_indv = phrase_buffer.indv();
            last_chrom = phrase_buffer.chrom();
            last_alele = phrase_buffer.alele();

            // Then create anothe completition factor for the start
            factor_pos = S_ref_actual_pos;
            factor_len = phrase_buffer.pos() + 1;
            factors.push_back(make_pair(factor_pos, factor_len));
            cout << " ICmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
            cout << "\tState: " << S_ref_actual_pos << endl;

            S_ref_actual_pos += factor_len;
            pos_ref_consumed += factor_len;
            S_size += factor_len;
        }
        else
        {
            // Just create a completition factor if it required
            if (phrase_buffer.pos() - last_pos > 1)
            {
                factor_pos = S_ref_actual_pos;
                factor_len = phrase_buffer.pos() + 1;
                factors.push_back(make_pair(factor_pos, factor_len));
                cout << " MCmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                cout << "\tState: " << S_ref_actual_pos << endl;

                S_ref_actual_pos += factor_len;
                pos_ref_consumed += factor_len;
                S_size += factor_len;
            }
        }
        // Then create the actual factor
        factors.push_back(make_pair(phrase_buffer.pos_e(), phrase_buffer.len_e()));
        cout << "\t Actu: (" << phrase_buffer.pos_e() << "," << phrase_buffer.len_e() << ")" << endl;
        cout << "\t\tState: " << S_ref_actual_pos << endl;

        S_size += phrase_buffer.len_e();

        // Update variables for the next iteration
        last_pos = phrase_buffer.pos();
        last_l = phrase_buffer.len();
        last_l_e = phrase_buffer.len_e();
        // We need to discard the slice of ref not used in this edit
        S_ref_actual_pos += phrase_buffer.len();
        pos_ref_consumed += phrase_buffer.len();
    }
    // For the last edit, do the same as we change to a new sample
    factor_pos = S_ref_actual_pos;
    factor_len = dict_metareference[last_chrom].n_bases() - (pos_ref_consumed + 1);
    factors.push_back(make_pair(factor_pos, factor_len));
    cout << " Last: ->(" << factor_pos << "," << factor_len << ")" << endl;
    cout << "\tState: " << S_ref_actual_pos << endl;

    S_ref_actual_pos += factor_len;
}