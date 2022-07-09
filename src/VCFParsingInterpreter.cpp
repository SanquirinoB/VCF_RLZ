#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

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

    // Open all required files
    // OpenBinaryInput(Reference_file, Reference_file_path);
    // cout << "2" << endl;
    // OpenBinaryInput(Phrases_file, Phrases_file_path);
    // OpenBinaryInput(MetaReference_file, MetaReference_file_path);
    // cout << "3" << endl;
    // OpenBinaryInput(MetaParsing_file, MetaParsing_file_path);
}

void VCFParsingInterpreter::OpenBinaryInput(ifstream *stream, string path)
{
    cout << path << endl;
    stream->open(path, ifstream::in | ifstream::binary);
}

void VCFParsingInterpreter::StartProcess()
{
    ProcessReference();
}

void VCFParsingInterpreter::ProcessReference()
{
    // Open required files
    Reference_file.open(Reference_file_path, ifstream::in | ifstream::binary);
    MetaReference_file.open(MetaReference_file_path, ifstream::in | ifstream::binary);

    // The first structure contains summary info
    cout << "Lee la estructura" << endl;
    MetaReference_file.read((char *)&metareference_buffer, sizeof(metareference));
    Reference_len = metareference_buffer.m_n_bases;
    n_References = metareference_buffer.m_rel_pos;

    // TEST
    cout << "Largo: " << Reference_len << "|Referencias: " << n_References << endl;
}

void VCFParsingInterpreter::buildFactorFromVCFParserPhrase(char *phrases_path, vector<pair<unsigned int, unsigned int>> &factors)
{
}