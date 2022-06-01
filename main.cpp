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
#include <stxxl/ksort>
#include <stxxl/sort>
#include <stxxl/stable_ksort>
#include <stxxl/vector>
#include <iostream> 

typedef unsigned char three_d; // 1 byte
typedef stxxl::uint16 four_d;  // 2 byte
typedef unsigned ten_d;        // 4 bytes

struct phrase
{

    // Optimal aligment
    four_d m_indv;
    three_d m_chrom, m_alele;
    ten_d m_pos;
    ten_d m_pos_e;
    char m_edit[4], m_len[6], m_len_e[6];

    four_d indv() const { return m_indv; }
    three_d chrom() const { return m_chrom; }
    three_d alele() const { return m_alele; }
    ten_d pos() const { return m_pos; }
    ten_d pos_e() const { return m_pos_e; }
    char[] edit() const { return m_edit; }
    char[] len() const { return m_len; }
    char[] len_e() const { return m_len_e; }

    phrase() {}
    phrase(four_d v_indv, three_d v_chrom, three_d v_alele, ten_d v_pos,
           ten_d v_pos_e, char *v_edit, char *v_len, char *v_len_e) : indv(v_indv),
                                                                      chrom(v_chrom),
                                                                      alele(v_alele),
                                                                      pos(v_pos),
                                                                      pos_e(v_pos_e),
                                                                      edit(v_edit),
                                                                      len(v_len),
                                                                      len_e(v_len_e) {}

    static phrase min_value()
    {
        return phrase((four_d)0, (three_d)0, (three_d)0, (ten_d)0, (ten_d)0, '0', '0', '0');
    }

    static phrase max_value()
    {
        phrase((four_d)9999, (three_d)999, (three_d)99, (ten_d)9999999999, (ten_d)9999999999, '0', '0', '0');
    }
};

inline bool operator<(const phrase &a, const phrase &b)
{
    return (a.indv() < b.indv() &&
            a.chrom() < b.chrom() &&
            a.alele() < b.alele() &&
            a.pos() < b.pos() &&
            a.pos_e() < b.pos_e());
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
    o.write((char *) obj, sizeof(phrase));
    return o;
}

template <typename T>
T digit_cast(char *c_number, three_d n_digits)
{
    T result = 0;

    for (three_d i = 0; i < n_digits; i++)
    {
        result += (c_number[i] - '0') * pow(10, i);
    }

    return result;
}

int fphrase_to_sphrase(std::iostream &o, const stxxl::int64 n_phrases)
{
    char *line;
    four_d indv;
    three_d chrom, alele;
    ten_d pos;
    ten_d pos_e;
    char edit[4], len[6], len_e[6];

    o.read(line, 32);
    std::cout << "Size of line " << sizeof(*line) << std::endl;
    // Parse line
    // Start again
    o.seekg(ios_base::beg);
    // Read indv
    o.read(line, 4);
    indv = digit_cast<four_d>(line, 4);
    // Read chrom
    o.read(line, 3);
    chrom = digit_cast<three_d>(line, 3);
    // Read alele
    o.read(line, 2);
    alele = digit_cast<three_d>(line, 2);
    // Read pos
    o.read(line, 10);
    pos = digit_cast<ten_d>(line, 10);
    // Read len
    o.read(line, 6);
    len = *line;
    // Read edit
    o.read(line, 4);
    edit = *line;

    // Need to test type of line
    o.read(line, 3);
    if (strcmp( &line[2], '\n'))
    {
        // If is a short phrase
        pos_e = 0;
        len_e = '0';
    }
    else
    {
        // Roll back and continue
        o.seekg(o.tellg() - 3);
        // Read pos_e
        o.read(line, 10);
        pos_e = digit_cast<ten_d>(line, 10);
        // Read len_e
        o.read(line, 6);
        len_e = *line;
    }

    phrase test = {indv, chrom, alele, pos, pos_e, edit, len, len_e};

    std::cout << "Size of line in struct " << sizeof(test) << std::endl;

    o.close();

    return 0;
}

int mainu(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " action file" << std::endl;
        std::cout << "       where action is one of generate, sort, ksort, stable_sort, stable_ksort" << std::endl;
        return -1;
    }
    // Definimos el tamano del bloque de memoria
    const stxxl::unsigned_type block_size = sizeof(phrase) * 4096;

    // Si buscamos generar los
    if (strcmp(argv[1], "generate") == 0)
    {
        const phrase::key_type num_elements = 1 * 1024 * 1024;
        const stxxl::unsigned_type records_in_block = block_size / sizeof(phrase);
        stxxl::syscall_file f(argv[2], stxxl::file::CREAT | stxxl::file::RDWR);
        phrase *array = (phrase *)stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(block_size);
        memset(array, 0, block_size);

        phrase::key_type cur_key = num_elements;
        for (unsigned i = 0; i < num_elements / records_in_block; i++)
        {
            for (unsigned j = 0; j < records_in_block; j++)
                array[j].m_key = cur_key--;

            stxxl::request_ptr req = f.awrite((void *)array, stxxl::int64(i) * block_size, block_size);
            req->wait();
        }
        stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(array);
    }
    else
    {
#if STXXL_PARALLEL_MULTIWAY_MERGE
        STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif
        stxxl::syscall_file f(argv[2], stxxl::file::DIRECT | stxxl::file::RDWR);
        unsigned memory_to_use = 50 * 1024 * 1024;
        typedef stxxl::vector<phrase, 1, stxxl::lru_pager<8>, block_size> vector_type;
        vector_type v(&f);

        /*
        STXXL_MSG("Printing...");
        for(stxxl::int64 i=0; i < v.size(); i++)
            STXXL_MSG(v[i].key());
         */

        STXXL_MSG("Checking order...");
        STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "OK" : "WRONG"));

        STXXL_MSG("Sorting...");
        if (strcmp(argv[1], "sort") == 0)
        {
            stxxl::sort(v.begin(), v.end(), Cmp(), memory_to_use);
#if 0 // stable_sort is not yet implemented
        }
        else if (strcmp(argv[1], "stable_sort") == 0) {
            stxxl::stable_sort(v.begin(), v.end(), memory_to_use);
#endif
        }
        else if (strcmp(argv[1], "ksort") == 0)
        {
            stxxl::ksort(v.begin(), v.end(), memory_to_use);
        }
        else if (strcmp(argv[1], "stable_ksort") == 0)
        {
            stxxl::stable_ksort(v.begin(), v.end(), memory_to_use);
        }
        else
        {
            STXXL_MSG("Not implemented: " << argv[1]);
        }

        STXXL_MSG("Checking order...");
        STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "OK" : "WRONG"));
    }

    return 0;
}

int main()
{
    iostream *file;
    file.open("../VCF_files/Tmp/Parsing/test_4.tmprlz", ios::in);
    // std::cout << "Tamano de estructura phrase " << sizeof(phrase) << std::endl;
    // return 0;

    return fphrase_to_sphrase(file, 7);
}

// vim: et:ts=4:sw=4
