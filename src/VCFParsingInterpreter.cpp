#include "VCFParsingInterpreter.h"

VCFParsingInterpreter::VCFParsingInterpreter(char* destination_path)
{
    string Destination_path(destination_path);

    string Reference_file_path = Destination_path + Reference_file_subpath;
    string Phrases_file_path = Destination_path + Phrases_file_subpath;
    string MetaReference_file_path = Destination_path + MetaReference_file_subpath;
    string MetaParsing_file_path = Destination_path + MetaParsing_file_subpath;
    


} 