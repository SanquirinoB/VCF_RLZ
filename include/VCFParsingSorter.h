#ifndef _VCF_PARSING_SORTER_H
#define _VCF_PARSING_SORTER_H
#include <VCFCommon.h>

class VCFParsingSorter
{
public:
    VCFParsingSorter();

    ~VCFParsingSorter();

    void StartProcess(char **argv);
};

#endif