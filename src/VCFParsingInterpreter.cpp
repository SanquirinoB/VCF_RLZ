#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

#define BINARY_FILE (ifstream::in | ifstream::binary);

VCFParsingInterpreter::VCFParsingInterpreter(char *destination_path)
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

    // Open all required files
    // OpenBinaryInput(Reference_file, Reference_file_path);
    // cout << "2" << endl;
    // OpenBinaryInput(Phrases_file, Phrases_file_path);
    // OpenBinaryInput(MetaReference_file, MetaReference_file_path);
    // cout << "3" << endl;
    // OpenBinaryInput(MetaParsing_file, MetaParsing_file_path);
}

void VCFParsingInterpreter::StartProcess()
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
    MetaParsing_file.read((char *)metainfo_buffer, sizeof(metainfo));
    n_Phrases = metainfo_buffer.n_phrases();
    MetaParsing_file.read((char *)metainfo_buffer, sizeof(metainfo));
    n_Samples = metainfo_buffer.n_phrases();
    MetaParsing_file.close();
}

void VCFParsingInterpreter::buildFactorFromVCFParserPhrase(vector<pair<unsigned int, unsigned int>> &factors)
{
    ProcessReference();
    ProcessMetaParsing();
    // Reference_file.open(Reference_file_path, ifstream::in | ifstream::binary);
    Phrases_file.open(Phrases_file_path, BINARY_FILE);

    int full_size = 0, last_pos = 0, curr_pos = 0, last_l = 0, last_l_e = 0;
    vector<int> end_pos;

    for (ten_d i = 0; i < n_Phrases; i++)
    {
        // TODO: Check Indv|Chrom|Alele as ID
        Phrases_file.read((char *)&phrase_buffer, sizeof(phrase));
        if (i == 0)
        {
            last_pos = phrase_buffer.pos();
            last_l = phrase_buffer.len();
            last_l_e = phrase_buffer.len_e();
            continue;
        }
        // Problema: este calculo asume que los edits son continuos entre si, falta contemplar que hay que mantener la parte de la referencia
        full_size += (phrase_buffer.pos() - last_l) - last_l + last_l_e;
        // Update values for next iteration
        last_pos = phrase_buffer.pos();
        last_l = phrase_buffer.len();
        last_l_e = phrase_buffer.len_e();
    }
}