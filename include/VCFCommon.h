#ifndef _VCF_COMMON_H
#define _VCF_COMMON_H

#include <stxxl/bits/common/types.h>
#include <stxxl/vector>
#include <iostream>
#include <cstdint>

typedef long long ll; // 4 bytes

struct metareference
{
    ll m_ID;
    ll m_n_bases;
    ll m_rel_pos;

    ll ID() const { return m_ID; }
    ll n_bases() const { return m_n_bases; }
    ll rel_pos() const { return m_rel_pos; }

    metareference() {}
    metareference(ll ID, ll n_bases, ll rel_pos) : m_ID(ID), m_n_bases(n_bases), m_rel_pos(rel_pos) {}
};

struct metainfo
{
    ll m_n_phrases;
    ll n_phrases() const { return m_n_phrases; }
    metainfo() {}
    metainfo(ll n_phrases) : m_n_phrases(n_phrases) {}
};

struct phrase
{
    ll m_indv, m_chrom, m_alele;
    ll m_pos, m_pos_e;
    ll m_edit;
    ll m_len, m_len_e;

    ll indv() const { return m_indv; }
    ll chrom() const { return m_chrom; }
    ll alele() const { return m_alele; }
    ll pos() const { return m_pos; }
    ll pos_e() const { return m_pos_e; }
    ll edit() const { return m_edit; }
    ll len() const { return m_len; }
    ll len_e() const { return m_len_e; }

    phrase() {}
    phrase(ll v_indv, ll v_chrom, ll v_alele, ll v_pos,
           ll v_pos_e, ll v_edit, ll v_len, ll v_len_e) : m_indv(v_indv),
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
        return phrase((ll)0, (ll)0, (ll)0, (ll)0, (ll)0, (ll)0, (ll)0, (ll)0);
    }

    static phrase max_value()
    {
        return phrase((ll)9999, (ll)9999, (ll)9999, (ll)9999999999, (ll)9999999999, (ll)0, (ll)0, (ll)0);
    }
};

typedef stxxl::vector<phrase, 1, stxxl::lru_pager<8>, sizeof(phrase) * 4096> vector_type;

inline bool operator<(const phrase &a, const phrase &b)
{
    if (a.indv() == b.indv())
    {
        if (a.chrom() == b.chrom())
        {
            if (a.alele() == b.alele())
            {
                return a.pos() < b.pos();
                // if (a.pos() == b.pos())
                // {
                //     return a.pos_e() < b.pos_e();
                // }
                // else
                // {
                //     return a.pos() < b.pos();
                // }
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

inline std::ostream &operator<<(std::ostream &o, const phrase &obj)
{
    phrase suitable_obj = obj;
    o.write((char *)&suitable_obj, sizeof(phrase));
    return o;
}

#endif // _VCF_COMMON_H