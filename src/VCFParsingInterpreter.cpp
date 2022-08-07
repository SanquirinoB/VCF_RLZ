#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

#define BINARY_FILE ifstream::in | ifstream::binary

VCFParsingInterpreter::VCFParsingInterpreter(char *destination_path)
{
    // Something like .../VCF_files/
    string Destination_path(destination_path);

    // TODO: Assert if ./Tmp is a valid folder
    Reference_file_path = Destination_path + Reference_file_subpath;
    Phrases_file_path = Destination_path + Phrases_file_subpath;
    MetaReference_file_path = Destination_path + MetaReference_file_subpath;
    MetaParsing_file_path = Destination_path + MetaParsing_file_subpath;
    IDInfo_file_path = Destination_path + IDInfo_file_subpath;

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
        if (i <= 4)
        {
            cout << "Ref: " << metareference_buffer.ID() << ", n_bases: " << metareference_buffer.n_bases()
                 << ", rel_pos: " << metareference_buffer.rel_pos() << endl;
        }
    }
    MetaReference_file.clear();
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
    MetaParsing_file.close();
}

void VCFParsingInterpreter::UpdateSampleData(ll index, phrase data)
{
    dict_samples[index] = sampleID(data.indv(), data.chrom(), data.alele());
}

void VCFParsingInterpreter::AddFactor(vector<pair<unsigned int, unsigned int>> &factors, ll pos, ll len)
{
    factors.push_back(make_pair((unsigned int)pos, (unsigned int)len));
}

bool VCFParsingInterpreter::NeedsToFullFill(phrase ref, phrase curr)
{
    bool same_sample = ref.indv() == curr.indv();
    bool same_chrom = ref.chrom() == curr.chrom();
    bool same_alele = ref.alele() == curr.alele();

    bool is_change_in_same_alele = (same_sample && same_chrom && same_alele);

    bool is_contigous_alele = (same_sample && same_chrom && ((ref.alele() == 1) && (curr.alele() == 0)));

    bool is_contigous_end_chrom = (same_sample &&
                                   (ref.chrom() == curr.chrom() - 1) &&
                                   ((ref.alele() == 1) && (curr.alele() == 0)));

    bool is_contigous_end_sample = ((ref.indv() == curr.indv() - 1) &&
                                    ((ref.chrom() == 23) && (curr.chrom() == 1)) &&
                                    ((ref.alele() == 1) && (curr.chrom() == 0)));

    return !(is_change_in_same_alele || is_contigous_alele || is_contigous_end_chrom || is_contigous_end_sample);
}

ll VCFParsingInterpreter::buildFactorFromVCFParserPhrase(vector<pair<unsigned int, unsigned int>> &factors)
{
    // Position variables: Offset checkpoints
    ll last_pos = 0, last_l = 0, last_l_e = 0;
    // Counter variables: S info
    vector<ll> S_i_pos;
    ll n_S_i = 0;

    // Helper variables: S construction
    ll ref_offset = 0, S_size = 0, factor_pos = 0, factor_len = 0, rel_pos_chrom = 0;
    // Helper variables: Phrase ID
    ll last_indv = -1, last_chrom = -1, last_alele = -1;
    bool is_same_sample = false;

    // PRINT
    bool is_debug = false;

    // Helpers
    phrase buffer;
    phrase curr_phrase;

    // Initialization
    Phrases_file.open(Phrases_file_path, BINARY_FILE);
    Phrases_file.read((char *)&curr_phrase, sizeof(phrase));
    // curr_phrase = buffer;

    if (is_debug)
    {
        cout << curr_phrase.indv() << "|" << curr_phrase.chrom() << "|" << curr_phrase.alele() << "|" << curr_phrase.pos() << "|"
             << curr_phrase.len() << "|" << curr_phrase.edit() << "|" << curr_phrase.pos_e() << "|" << curr_phrase.len_e() << endl;
        cout << "State: " << last_pos << "," << last_l << "," << last_l_e << endl;
    }

    last_indv = curr_phrase.indv();
    last_chrom = curr_phrase.chrom();
    last_alele = curr_phrase.alele();
    rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();

    if (curr_phrase.pos() != 0)
    {
        // If the first phrase edit doesnt occurr in the first base
        // we need to induce a completition factor
        factor_pos = rel_pos_chrom;
        factor_len = curr_phrase.pos() - 1;
        AddFactor(factors, factor_pos, factor_len);
        ref_offset += factor_len;
        S_size += factor_len;

        if (is_debug)
        {
            cout << " Init: (" << factor_pos << "," << factor_len << ")->" << endl;
        }
    }
    // And now create the factor associated to the current phrase
    AddFactor(factors, curr_phrase.pos_e(), curr_phrase.len_e());
    S_size += curr_phrase.len_e();

    if (is_debug)
    {
        cout << " Actu: (" << curr_phrase.pos_e() << "," << curr_phrase.len_e() << ")" << endl;
    }

    // Update with current values
    last_pos = curr_phrase.pos();
    last_l = curr_phrase.len();
    last_l_e = curr_phrase.len_e();
    rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();

    n_S_i += 1;
    UpdateSampleData(n_S_i, curr_phrase);
    S_i_pos.push_back(S_size);

    for (ll i = 1; i < n_Phrases; i++)
    {
        // For each phrase
        Phrases_file.read((char *)&curr_phrase, sizeof(phrase));
        // curr_phrase = buffer;
        is_same_sample = (last_indv == curr_phrase.indv() &&
                          last_chrom == curr_phrase.chrom() &&
                          last_alele == curr_phrase.alele());
        if (is_debug)
        {
            cout << curr_phrase.indv() << "|" << curr_phrase.chrom() << "|" << curr_phrase.alele() << "|" << curr_phrase.pos() << "|"
                 << curr_phrase.len() << "|" << curr_phrase.edit() << "|" << curr_phrase.pos_e() << "|" << curr_phrase.len_e() << endl;
            cout << "State: " << last_pos << "," << last_l << "," << last_l_e << endl;
        }

        // Discard slice not required from the previous phrase and go to the next position
        //  We only wont do it if there are multiple edits over the same position (sequence)
        if (!(is_same_sample && (last_pos == curr_phrase.pos())))
        {
            ref_offset += last_l + 1;
        }

        // If the current sample changes
        if (!is_same_sample)
        {
            // Create an end completition factor for the last sample
            factor_pos = rel_pos_chrom + ref_offset;
            factor_len = dict_metareference[last_chrom].n_bases() - (ref_offset + 1);
            AddFactor(factors, factor_pos, factor_len);

            S_size += factor_len;

            // Fill data info
            n_S_i += 1;
            UpdateSampleData(n_S_i, curr_phrase);
            S_i_pos.push_back(S_size);

            if (is_debug)
            {
                cout << " ECmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                cout << "\tState: " << ref_offset << endl;
            }

            // If we skip one alele, we need to induce it
            // We skip first alele of the new chromosome
            if (last_alele == 1 && curr_phrase.alele() == 1)
            {
                factor_pos = dict_metareference[curr_phrase.chrom()].rel_pos();
                factor_len = dict_metareference[curr_phrase.chrom()].n_bases() - 1;
                AddFactor(factors, factor_pos, factor_len);
                S_size += factor_len;

                if (is_debug)
                {
                    cout << " Indc: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                }
            }
            // Or we skip the second alele of the last chromosome
            else if (last_alele == 0 && curr_phrase.alele() == 0)
            {
                factor_pos = dict_metareference[last_chrom].rel_pos();
                factor_len = dict_metareference[last_chrom].n_bases() - 1;
                AddFactor(factors, factor_pos, factor_len);
                S_size += factor_len;

                if (is_debug)
                {
                    cout << " Indc: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                }
            }

            // Update positon checkpoint variables
            last_pos = 0;
            last_l = 0;
            last_l_e = 0;
            ref_offset = 0;

            // Update ID variables
            last_indv = curr_phrase.indv();
            last_chrom = curr_phrase.chrom();
            last_alele = curr_phrase.alele();
            rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();

            // Then create anothe completition factor for the start
            factor_pos = rel_pos_chrom + ref_offset;
            factor_len = curr_phrase.pos() - 1;
            AddFactor(factors, factor_pos, factor_len);

            if (is_debug)
            {
                cout << " ICmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                cout << "\tState: " << ref_offset << endl;
            }

            ref_offset += factor_len;
            S_size += factor_len;
        }
        else
        {
            // Just create a completition factor if it required
            if (curr_phrase.pos() - last_pos > 1)
            {
                factor_pos = rel_pos_chrom + ref_offset;
                factor_len = (curr_phrase.pos() - 1) - ref_offset;
                AddFactor(factors, factor_pos, factor_len);
                if (is_debug)
                {
                    cout << " MCmp: ->(" << factor_pos << "," << factor_len << ")->" << endl;
                    cout << "\tState: " << ref_offset << endl;
                }
                ref_offset += factor_len;
                S_size += factor_len;
            }
        }
        // Then create the actual factor
        AddFactor(factors, curr_phrase.pos_e(), curr_phrase.len_e());
        if (is_debug)
        {
            cout << "\t Actu: (" << curr_phrase.pos_e() << "," << curr_phrase.len_e() << ")" << endl;
            cout << "\t\tState: " << ref_offset << endl;
        }

        S_size += curr_phrase.len_e();

        // Update variables for the next iteration
        last_pos = curr_phrase.pos();
        last_l = curr_phrase.len();
        last_l_e = curr_phrase.len_e();
    }
    // For the last edit, do the same as we change to a new sample
    factor_pos = rel_pos_chrom + ref_offset;
    factor_len = dict_metareference[last_chrom].n_bases() - (ref_offset + 1);
    AddFactor(factors, factor_pos, factor_len);
    if (is_debug)
    {
        cout << " Last: ->(" << factor_pos << "," << factor_len << ")" << endl;
        cout << "\tState: " << ref_offset << endl;
    }
    S_size += factor_len;

    n_S_i += 1;
    UpdateSampleData(n_S_i, curr_phrase);
    S_i_pos.push_back(S_size);

    Phrases_file.clear();
    Phrases_file.close();

    return S_size;
}

pair<const char *, ll> VCFParsingInterpreter::GetReference()
{
    cout << 1 << endl;
    Reference_file.open(Reference_file_path, ios::in);
    cout << Reference_len << endl;

    string reference_aux;
    if (Reference_file.is_open())
    {
        getline(Reference_file, reference_aux);
    }
    else
    {
        cout << 69 << endl;
    }
    cout << 3 << endl;

    Reference_file.clear();
    Reference_file.close();
    cout << 5 << endl;

    const char *reference = reference_aux.c_str();

    return make_pair(reference, Reference_len);
}