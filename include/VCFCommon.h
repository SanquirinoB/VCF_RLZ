#include <stxxl/bits/common/types.h>

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
        return phrase((four_d)9999, (four_d)9999, (four_d)9999, (ten_d)9999999999, (ten_d)9999999999, (four_d)0, (ten_d)0, (ten_d)0);
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