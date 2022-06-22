#include "VCFParsingSorter.h"

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/vector>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <VCFCommon.h>

VCFParsingSorter::VCFParsingSorter(){}

VCFParsingSorter::~VCFParsingSorter(){}

int VCFParsingSorter::StartProcess(char **argv)
{
     // INICIO PROCESO DE SORTING

    std::ifstream meta_info_file;
    metainfo info;
    phrase p1;

    std::string path_raw(argv[1]);
    std::string name_meta("Tmp/Meta_data/Meta_info.metarlz");
    std::string name_parsing("Tmp/Parsing/phrases.tmprlz");

    meta_info_file.open(path_raw + name_meta, std::ifstream::in | std::ifstream::binary);
    meta_info_file.read((char *)&info, sizeof(metainfo));

    std::cout << "[RLZ] Number of identified edits: " << info.n_phrases() << std::endl;

#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif

    stxxl::syscall_file f(path_raw + name_parsing, stxxl::file::DIRECT | stxxl::file::RDWR);

    const stxxl::unsigned_type block_size = sizeof(phrase) * 4096;
    // TODO: Cual es sentido del nro extra?
    unsigned memory_to_use = 16 * sizeof(phrase) * info.n_phrases();

    typedef stxxl::vector<phrase, 1, stxxl::lru_pager<8>, block_size> vector_type;
    vector_type v(&f);

    STXXL_MSG("Checking order...");
    STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "\tOK" : "\tWRONG"));

    std::cout << "[RLZ] Unsorted state sample:" << std::endl;
    std::cout << "\tIndv|Chrom|Alele|Pos|Len|Edit|Pos_E|Len_E" << std::endl;
    for (int i = 1; i < 10; i++)
    {
        p1 = v[i];
        std::cout << "\t" << p1.indv() << "|" << p1.chrom() << "|" << p1.alele() << "|" << p1.pos() << "|"
                  << p1.len() << "|" << p1.edit() << "|" << p1.pos_e() << "|" << p1.len_e() << std::endl;
    }

    STXXL_MSG("Sorting...");
    begin = clock();
    stxxl::sort(v.begin(), v.end(), Cmp(), memory_to_use);
    end = clock();
    std::cout << "[RLZ]\tSort sucessfull! Time elpased: " << (float)(end - begin)/CLOCKS_PER_SEC << std::endl;

    STXXL_MSG("Checking order...");
    STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "\tOK" : "\tWRONG"));

    std::cout << "[RLZ] Sorted state sample:" << std::endl;
    std::cout << "\tIndv|Chrom|Alele|Pos|Len|Edit|Pos_E|Len_E" << std::endl;
    for (int i = 1; i < 10; i++)
    {
        p1 = v[i];
        std::cout << "\t" << p1.indv() << "|" << p1.chrom() << "|" << p1.alele() << "|" << p1.pos() << "|"
                  << p1.len() << "|" << p1.edit() << "|" << p1.pos_e() << "|" << p1.len_e() << std::endl;
    }

    return 0;
}