/***************************************************************************
 *  examples/algo/sort_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example algo/sort_file.cpp
//! This example imports a file into an \c stxxl::vector without copying its
//! content and then sorts it using stxxl::sort / stxxl::ksort / ...

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/vector>
#include <iostream>
#include <fstream>
#include <string>

typedef stxxl::uint16 four_d; // 2 byte
typedef unsigned ten_d;       // 4 bytes

struct metainfo
{
    ten_d m_n_phrases;
    ten_d n_phrases() const { return m_n_phrases; }
    metainfo() {}
    metainfo(ten_d n_phrases) : m_n_phrases(n_phrases) {}
};

struct phrase
{

    // Optimal aligment
    four_d m_indv, m_chrom, m_alele;
    ten_d m_pos, m_pos_e;
    four_d m_edit;
    ten_d m_len, m_len_e;

    four_d indv() const { return m_indv; }
    four_d chrom() const { return m_chrom; }
    four_d alele() const { return m_alele; }
    ten_d pos() const { return m_pos; }
    ten_d pos_e() const { return m_pos_e; }
    four_d edit() const { return m_edit; }
    ten_d len() const { return m_len; }
    ten_d len_e() const { return m_len_e; }

    phrase() {}
    phrase(four_d v_indv, four_d v_chrom, four_d v_alele, ten_d v_pos,
           ten_d v_pos_e, four_d v_edit, ten_d v_len, ten_d v_len_e) : m_indv(v_indv),
                                                                       m_chrom(v_chrom),
                                                                       m_alele(v_alele),
                                                                       m_pos(v_pos),
                                                                       m_pos_e(v_pos_e),
                                                                       m_edit(v_edit),
                                                                       m_len(v_len),
                                                                       m_len_e(v_len_e)
    {
    }

    static phrase min_value()
    {
        return phrase((four_d)0, (four_d)0, (four_d)0, (ten_d)0, (ten_d)0, (four_d)0, (ten_d)0, (ten_d)0);
    }

    static phrase max_value()
    {
        return phrase((four_d)9999, (four_d)999, (four_d)99, (ten_d)9999999999, (ten_d)9999999999, (four_d)0, (ten_d)0, (ten_d)0);
    }
};

inline bool operator<(const phrase &a, const phrase &b)
{
    if (a.indv() == b.indv())
    {
        if (a.chrom() == b.chrom())
        {
            if (a.alele() == b.alele())
            {
                if (a.pos() == b.pos())
                {
                    return a.pos_e() < b.pos_e();
                }
                else
                {
                    return a.pos() < b.pos();
                }
            }
            else
            {
                return a.alele() < b.alele();
            }
        }
        else
        {
            return a.chrom() < b.chrom();
        }
    }
    else
    {
        return a.indv() < b.indv();
    }
}

inline bool operator==(const phrase &a, const phrase &b)
{
    return (a.indv() == b.indv() &&
            a.chrom() == b.chrom() &&
            a.alele() == b.alele() &&
            a.pos() == b.pos() &&
            a.pos_e() == b.pos_e()); // Tecnically is not possible to have different edit/len/len_e if the rest are the same
}

struct Cmp
{
    typedef phrase first_argument_type;
    typedef phrase second_argument_type;
    typedef bool result_type;
    bool operator()(const phrase &a, const phrase &b) const
    {
        return a < b;
    }
    static phrase min_value()
    {
        return phrase::min_value();
    }
    static phrase max_value()
    {
        return phrase::max_value();
    }
};

std::ostream &operator<<(std::ostream &o, const phrase &obj)
{
    phrase suitable_obj = obj;
    o.write((char *)&suitable_obj, sizeof(phrase));
    return o;
}

template <typename T>
T digit_cast(char *c_number, four_d n_digits)
{
    T result = 0;

    for (four_d i = 0; i < n_digits; i++)
    {
        result += (c_number[i] - '0') * pow(10, i);
    }

    return result;
}

int main(int argc, char **argv)
{
    // TODO: Larger version
    // main.cpp ../VCF_files/ -n 1 ../VCF_files/test_4.vcf
    // For now
    if (argc < 5)
    {
        std::cout << "Usage: " << argv[0] << " destination_folder -n [NUMBER] [FILES] [OPTIONS]" << std::endl;
        std::cout << "       where [NUMBER] is the number of [FILES] to be given." << std::endl;
        return -1;
    }
    // INICIO PROCESO DE PARSING VIA PYTHON
    std::stringstream aux;
    aux << argv[3];
    int files_expected;
    aux >> files_expected;

    std::vector<std::string> py_params;

    py_params.push_back("python3 ./VCF_parsing/parsing_process.py");
    py_params.push_back(" ");
    py_params.push_back(argv[1]); // Destination folder
    py_params.push_back(" ");
    py_params.push_back(argv[2]); // -n
    py_params.push_back(" "); 
    py_params.push_back(argv[3]); // [NUMBER]

    for (int i = 4; i < 4 + files_expected; i++)
    {
        py_params.push_back(" ");
        py_params.push_back(argv[i]); // [FILES]
    }
    // TODO: Remain to include parsing options

    // Example: python3 parsing_process.py ../VCF_files/ -n 1 ../VCF_files/test_4.vcf
    std::string command = std::accumulate(py_params.begin(), py_params.end(), std::string(""));
    std::cout << "[RLZ] Start parsing process..." << std::endl;
    if (system(command.c_str()) == 0)
    {
        std::cout << "[RLZ]\tParsing sucessfull!" << std::endl;
    } else
    {
        std::cout << "[RLZ]\tParsing failed, check log for more information :(" << std::endl;
        return -1;
    }

    // FIN PROCESO DE PARSING VIA PYTHON

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
    stxxl::sort(v.begin(), v.end(), Cmp(), memory_to_use);

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

// vim: et:ts=4:sw=4
