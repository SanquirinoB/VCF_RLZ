#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/vector>
#include <iostream>
#include <fstream>
#include <string>
#include "VCFParsingSorter.h"
#include "NanoTimer.h"

using namespace std;

VCFParsingSorter::VCFParsingSorter() {}

VCFParsingSorter::~VCFParsingSorter() {}

void VCFParsingSorter::StartProcess(char **argv)
{
    string path_raw(argv[1]);
    StartProcess(path_raw);
}

void VCFParsingSorter::StartProcess(string path_raw)
{
    // INICIO PROCESO DE SORTING
    NanoTimer timer;
    ifstream meta_info_file;
    metainfo info;
    phrase p1;

    string name_meta("Tmp/Meta_data/Meta_info.metarlz");
    string name_parsing("Tmp/Parsing/phrases.tmprlz");
    string name_parsing_sorted("Tmp/Parsing/phrases_sorted.tmprlz");

    meta_info_file.open(path_raw + name_meta, ifstream::in | ifstream::binary);
    meta_info_file.read((char *)&info, sizeof(metainfo));
    meta_info_file.close();

    cout << "[RLZ] Number of identified edits: " << info.n_phrases() << endl;

#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif

    stxxl::syscall_file f(path_raw + name_parsing, stxxl::file::DIRECT | stxxl::file::RDWR);

    // const stxxl::unsigned_type block_size = sizeof(phrase) * 4096;
    // TODO: Cual es sentido del nro extra?
    ull memory_to_use = 16ULL * sizeof(phrase) * ((ull) info.n_phrases());
    // typedef stxxl::vector<phrase, 1, stxxl::lru_pager<8>, sizeof(phrase) * 4096> vector_type;
    vector_type v(&f);

    STXXL_MSG("Checking order...");
    STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "\tOK" : "\tWRONG"));

    // cout << "[RLZ] Unsorted state sample:" << endl;
    // cout << "\tIndv|Chrom|Alele|Pos|Len|Edit|Pos_E|Len_E" << endl;
    // for (int i = 0; i < 10; i++)
    // {
    //     p1 = v[i];
    //     cout << "\t" << p1.indv() << "|" << p1.chrom() << "|" << p1.alele() << "|" << p1.pos() << "|"
    //          << p1.len() << "|" << p1.edit() << "|" << p1.pos_e() << "|" << p1.len_e() << endl;
    // }

    STXXL_MSG("Sorting...");
    timer.reset();
    stxxl::sort(v.begin(), v.end(), Cmp(), memory_to_use);
    cout << "[RLZ]\tSort sucessfull! Time elpased: " << timer.getMilisec() << "ms" << endl;
    // cout << "[RLZ]\tSort sucessfull! Time elpased "<< timer.getMilisec << endl;

    STXXL_MSG("Checking order...");
    STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "\tOK" : "\tWRONG"));

    // cout << "[RLZ] Sorted state sample:" << endl;
    // cout << "\tIndv|Chrom|Alele|Pos|Len|Edit|Pos_E|Len_E" << endl;
    // for (int i = 0; i < 10; i++)
    // {
    //     p1 = v[i];
    //     cout << "\t" << p1.indv() << "|" << p1.chrom() << "|" << p1.alele() << "|" << p1.pos() << "|"
    //          << p1.len() << "|" << p1.edit() << "|" << p1.pos_e() << "|" << p1.len_e() << endl;
    // }

    ofstream sorted_phrases;
    sorted_phrases.open(name_parsing_sorted, ofstream::out | ofstream::binary);
    for (ll i = 0; i < info.n_phrases(); i++)
    {
        sorted_phrases.write((char *)&v[i], sizeof(phrase));
    }
    sorted_phrases.close();
}