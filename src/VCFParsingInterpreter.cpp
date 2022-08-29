#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

#define BINARY_FILE ifstream::in | ifstream::binary

VCFParsingInterpreter::VCFParsingInterpreter(char *destination_path)
{
    // Something like .../VCF_files/
    string Destination_path(destination_path);

    // TODO: Assert if ./Tmp is a valid folder
    Reference_file_path = Destination_path + Reference_file_subpath;
    Phrases_file_path = Destination_path + Phrases_file_subpath;
    MetaReference_file_path = Destination_path + MetaReference_file_subpath;
    MetaParsing_file_path = Destination_path + MetaParsing_file_subpath;
    IDInfo_file_path = Destination_path + IDInfo_file_subpath;
}

void VCFParsingInterpreter::Initialize()
{
    cout << "[RLZ] Reading metadata" << endl;
    NanoTimer timer;
    ProcessReference();
    ProcessMetaParsing();
    cout << "[RLZ]      Metadata readed in " << timer.getMilisec() << " ms" << endl;
    cout << "[RLZ] Building factors from phrases" << endl;
    timer.reset();
    BuildFactors();
    cout << "[RLZ] Building index" << endl;
    char* reference = GetReference();
    timer.reset();
    for (pair<unsigned int, unsigned int> factor:factors)
    {
        cout << "(" << factor.first << "," << factor.second << ")" << endl;
    }
    cout << reference << endl;
    Index = new RelzIndexReference(factors, reference, S_size, reference, Reference_len);
    cout << "[RLZ]      Index finished in " << timer.getMilisec() << " ms" << endl;
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
        dict_metareference[metareference_buffer.ID()] = metareference_buffer;
        // if (i <= 4)
        // {
        //     cout << "Ref: " << metareference_buffer.ID() << ", n_bases: " << metareference_buffer.n_bases()
        //          << ", rel_pos: " << metareference_buffer.rel_pos() << endl;
        // }
    }
    MetaReference_file.clear();
    MetaReference_file.close();
}

void VCFParsingInterpreter::ProcessMetaParsing()
{
    MetaParsing_file.open(MetaParsing_file_path, BINARY_FILE);
    MetaParsing_file.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_Phrases = metainfo_buffer.n_phrases();
    MetaParsing_file.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_Samples = metainfo_buffer.n_phrases();
    MetaParsing_file.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_chromosomes = metainfo_buffer.n_phrases();
    MetaParsing_file.close();
    MetaParsing_file.close();
}

void VCFParsingInterpreter::HigienicFactorPushBack(pair<unsigned int, unsigned int> factor)
{
    if (factor.second > 0)
    {
        factors.push_back(make_pair((unsigned int)factor.first, (unsigned int)factor.second));
        S_size += factor.second;
    }
}

bool VCFParsingInterpreter::ChangeInSameIndvChromAlele(phrase ref, phrase curr)
{
    bool same_sample = ref.indv() == curr.indv();
    bool same_chrom = ref.chrom() == curr.chrom();
    bool same_alele = ref.alele() == curr.alele();

    return (same_sample && same_chrom && same_alele);
}

bool VCFParsingInterpreter::ChangeInSameIndvChromDiffAlele(phrase ref, phrase curr)
{
    bool same_sample = ref.indv() == curr.indv();
    bool same_chrom = ref.chrom() == curr.chrom();
    return (same_sample && same_chrom && ((ref.alele() == 0) && (curr.alele() == 1)));
}

bool VCFParsingInterpreter::ChangeInSameIndvDiffChrom(phrase ref, phrase curr)
{
    bool same_sample = ref.indv() == curr.indv();
    bool diff_chrom = ref.chrom() != curr.chrom();
    return (same_sample && diff_chrom);
}

bool VCFParsingInterpreter::ChangeInDiffIndv(phrase ref, phrase curr)
{
    return ref.indv() != curr.indv();
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateAleleInitFactor(phrase Phrase)
{
    rel_pos_chrom = dict_metareference[Phrase.chrom()].rel_pos();
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom,     // Init position of current reference chromosome
                                                    Phrase.pos() - 1); // Copy until reach the actual position
    // cout << "(" << aux.first << "," << aux.second << ") : init" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateAleleInterFactor(phrase Phrase1, phrase Phrase2)
{
    rel_pos_chrom = dict_metareference[Phrase1.chrom()].rel_pos();
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom + Phrase1.pos() + Phrase1.len(),          // Start after the discarded slice of Phrase1
                     (Phrase2.pos() - 1) - (Phrase1.pos() + Phrase1.len())); // Fill until reach pos of Phrase2
    // cout << "(" << aux.first << "," << aux.second << ") : inter" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateAleleEndFactor(phrase Phrase)
{
    rel_pos_chrom = dict_metareference[Phrase.chrom()].rel_pos();
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom + Phrase.pos() + Phrase.len(),                                   // Start after the discarded slice of Phrase
                     dict_metareference[Phrase.chrom()].n_bases() - (Phrase.pos() + Phrase.len())); // Fill until reach the end of the current chromosome-alele
    // cout << "(" << aux.first << "," << aux.second << ") : end" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateFullAleleFactor(ll chromosome)
{
    rel_pos_chrom = dict_metareference[chromosome].rel_pos();
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom,                                 // Init position of current reference chromosome
                     dict_metareference[chromosome].n_bases()); // Fill until reach the end of the current chromosome-alele
    // cout << "(" << aux.first << "," << aux.second << ") : full" << endl;
    return aux;
}

void VCFParsingInterpreter::InduceFillFactors(phrase last_phrase, phrase curr_phrase, bool last_is_dummy, bool curr_is_dummy)
{
    // Consideration, we asume for now all chromosomes are diploid
    // Tecnically we dont properly process X and Y chromosomes, we will treat them as any other chromosome
    //  See https://www.biostars.org/p/59419/ , i dont have time for implement this handling or even understand it
    // TODO: Induce checkopoints over samples
    
    pair<unsigned int, unsigned int> aux_factor;
    if (ChangeInSameIndvChromAlele(last_phrase, curr_phrase))
    {
        if (last_phrase.pos() == curr_phrase.pos()) return;

        if (last_is_dummy)
        {
            // We need to induce the init factor (0, distanceToCurrPos)
            aux_factor = CalculateAleleInitFactor(curr_phrase);
            S_i_pos.push_back(S_size);
        }
        else if (curr_is_dummy)
        {
            // We need to induce an end factor, closing the current alele
            aux_factor = CalculateAleleEndFactor(last_phrase);
        }
        else
        {
            // Induce a factor which ignore the discarded section of lp and continue to curr pos
            aux_factor = CalculateAleleInterFactor(last_phrase, curr_phrase);
        }
        HigienicFactorPushBack(aux_factor);
        
    }
    else if (ChangeInSameIndvChromDiffAlele(last_phrase, curr_phrase))
    {
        if (last_is_dummy)
        {
            // We need to induce the full factor for the first alele
            aux_factor = CalculateFullAleleFactor(last_phrase.chrom());
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
            // Then induce the init factor for the current alele
            aux_factor = CalculateAleleInitFactor(curr_phrase);
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
        }
        else if (curr_is_dummy)
        {
            // We need to induce an end factor, closing the last alele
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);
            // We need to induce the full factor for the curr and last alele
            aux_factor = CalculateFullAleleFactor(last_phrase.chrom());
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
        }
        else
        {
            // We need to induce an end factor, closing the last alele
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);
            // We need to induce the init factor for the current alele
            aux_factor = CalculateAleleInitFactor(curr_phrase);
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
        }
    }
    else if (ChangeInSameIndvDiffChrom(last_phrase, curr_phrase))
    {
        int actual_chrom, actual_alele;
        actual_chrom = last_phrase.chrom();
        actual_alele = last_phrase.alele();

        if (!last_is_dummy)
        {
            // First induce an end fill
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);
            if (actual_alele == 1)
            {
                actual_chrom++;
                actual_alele--;
            }
            else
            {
                actual_alele++;
            }
        }

        // Full fill until reach the actual chrom and alele
        while (!(actual_chrom == curr_phrase.chrom() && actual_alele == curr_phrase.alele()))
        {
            aux_factor = CalculateFullAleleFactor(actual_chrom);
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
            if (actual_alele == 1)
            {
                actual_chrom++;
                actual_alele--;
            }
            else
            {
                actual_alele++;
            }
        }

        if (curr_is_dummy)
        {
            // Just complete the last alele
            aux_factor = CalculateFullAleleFactor(actual_chrom);
        }
        else
        {
            // Then induce an init factor for current phrase
            aux_factor = CalculateAleleInitFactor(curr_phrase);
        }
        S_i_pos.push_back(S_size);
        HigienicFactorPushBack(aux_factor);
        
    }
    else if (ChangeInDiffIndv(last_phrase, curr_phrase))
    {
        int actual_indv, actual_chrom, actual_alele;
        actual_indv = last_phrase.indv();
        actual_chrom = last_phrase.chrom();
        actual_alele = last_phrase.alele();

        if (!last_is_dummy)
        {
            // First induce an end fill
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);

            if (actual_alele == 1)
            {
                // If reach chrom limit
                if (actual_chrom == n_chromosomes - 1)
                {
                    // Move to next indv
                    actual_indv++;
                    // Reset chromosome
                    actual_chrom = 0;
                }
                else
                {
                    // Move to next chromosome
                    actual_chrom++;
                }
                // Reset alele
                actual_alele--;
            }
            else
            {
                actual_alele++;
            }
        }

        // Full fill until reach the actual indv, chrom and alele
        while (!(actual_indv == curr_phrase.indv() && actual_chrom == curr_phrase.chrom() && actual_alele == curr_phrase.alele()))
        {
            aux_factor = CalculateFullAleleFactor(actual_chrom);
            S_i_pos.push_back(S_size);
            HigienicFactorPushBack(aux_factor);
            // If reach alele limit
            if (actual_alele == 1)
            {
                // If reach chrom limit
                if (actual_chrom == n_chromosomes - 1)
                {
                    // Move to next indv
                    actual_indv++;
                    // Reset chromosome
                    actual_chrom = 0;
                }
                else
                {
                    // Move to next chromosome
                    actual_chrom++;
                }
                // Reset alele
                actual_alele--;
            }
            else
            {
                // Move to next alele
                actual_alele++;
            }
        }

        if (curr_is_dummy)
        {
            // Just complete the last alele
            aux_factor = CalculateFullAleleFactor(actual_chrom);
        }
        else
        {
            // Then induce an init factor for current phrase
            aux_factor = CalculateAleleInitFactor(curr_phrase);
        }
        
        S_i_pos.push_back(S_size);
        HigienicFactorPushBack(aux_factor);

    }
    else
    {
        cout << "Oh no, this is not supposed to happen" << endl;
    }
}

void VCFParsingInterpreter::BuildReconstructionStructures()
{
    bit_vector b(S_size, 0);
    // Turn on for all initial pos of each S_i
    for (ll i = 0; i < S_i_pos.size(); i++)
    {
        b[S_i_pos[i]] = 1;
    }

    bit_vector_S_i = new rrr_vector<127>(b);
    rank_S_i = new rrr_vector<127>::rank_1_type(bit_vector_S_i);
    select_S_i = new rrr_vector<127>::select_1_type(bit_vector_S_i);
}

void VCFParsingInterpreter::BuildFactors()
{
    // Phrase buffer
    phrase curr_phrase;
    // Actual size of S
    S_size = 0;
    n_S_i = 0;
    // Helpers
    ref_offset = 0;
    ll factor_pos = 0, factor_len = 0;
    bool last_is_dummy = true, curr_is_dummy = false;
    bool is_debug = true;

    // Open phrases file
    Phrases_file.open(Phrases_file_path, BINARY_FILE);
    // Read first dummy phrase
    Phrases_file.read((char *)&curr_phrase, sizeof(phrase));

    // Initialize indentifier variants
    phrase last_phrase = curr_phrase;
    rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();

    for (ll i = 1; i < n_Phrases; i++)
    {
        // Read the next phrase (Read the first actual phrase)
        Phrases_file.read((char *)&curr_phrase, sizeof(phrase));
        if (i == n_Phrases - 1)
        {
            curr_is_dummy = true;
        }

        // First check if with respect to the last phrase we need to induce full filling phrases
        InduceFillFactors(last_phrase, curr_phrase, last_is_dummy, curr_is_dummy);

        if (!curr_is_dummy)
        {
            // Add actual edit inside current phrase
            pair<unsigned int, unsigned int> aux = make_pair(curr_phrase.pos_e(), curr_phrase.len_e());
            // cout << "(" << aux.first << "," << aux.second << ") : actual" << endl;
            HigienicFactorPushBack(aux);
        }
        last_is_dummy = false;

        // Update variables for the next iteration
        last_phrase = curr_phrase;
        rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();
    }

    Phrases_file.clear();
    Phrases_file.close();

    int actual_indv = 0, actual_chrom = 0, actual_alele = 0;
    // for (int i = 0; i < S_i_pos.size(); i++)
    // {
    //     cout << "(" << actual_indv << "," << actual_chrom << "," << actual_alele << "): " << S_i_pos[i] << endl;
    //     // If reach alele limit
    //     if (actual_alele == 1)
    //     {
    //         // If reach chrom limit
    //         if (actual_chrom == n_chromosomes - 1)
    //         {
    //             // Move to next indv
    //             actual_indv++;
    //             // Reset chromosome
    //             actual_chrom = 0;
    //         }
    //         else
    //         {
    //             // Move to next chromosome
    //             actual_chrom++;
    //         }
    //         // Reset alele
    //         actual_alele--;
    //     }
    //     else
    //     {
    //         // Move to next alele
    //         actual_alele++;
    //     }
    // }

    BuildReconstructionStructures();
}

char* VCFParsingInterpreter::GetReference()
{
    Reference_file.open(Reference_file_path, ios::in);

    string reference_aux;
    if (Reference_file.is_open())
    {
        getline(Reference_file, reference_aux);
    }

    Reference_file.clear();
    Reference_file.close();

    const char *reference_const = reference_aux.c_str();
    char *reference = new char[Reference_len + 1];
    strcpy(reference, reference_const);

    return reference;
}

vector<pair<sampleID, unsigned int>> VCFParsingInterpreter::FindSnippet(string snippet)
{
    vector<unsigned int> raw_positions;
    vector<pair<sampleID, unsigned int>> result;
    Index->findTimes(snippet, raw_positions);
    unsigned int n_chrom = (unsigned int) n_chromosomes;
    unsigned int near_init, sample, chrom, alele, last_i, ini_s_i, aux_1, aux_2;
    string sample_name;

    IDInfo_file.open(IDInfo_file_path);

    sort(raw_positions.begin(), raw_positions.end());
    // Initialize ID
    last_i = 1;
    getline(IDInfo_file, sample_name);

    for(unsigned int raw_pos : raw_positions)
    {
        raw_pos += 1;
        near_init = (*rank_S_i)(raw_pos);
        // Cuantos Indv lleno hasta este punto?
        aux_1 = near_init / (n_chrom * ploidy);
        // Si lleno justo, entonces me quedo en aux, sino paso al sgte
        cout << (aux_1 * (n_chrom * ploidy) < near_init) << endl;
        sample = aux_1 * (n_chrom * ploidy) < near_init ? aux_1 + 1 : aux_1;

        if(sample != last_i)
        {
            getline(IDInfo_file, sample_name);
            last_i = sample;
        }

        // De lo que me sobro dentro de este sample
        aux_1 = near_init - aux_1 * (n_chrom * ploidy);
        // Arreglo en caso de llenar justo, equivale a ultima seccion
        aux_1 = aux_1 == 0? (n_chrom * ploidy) : aux_1;
        // De lo que sobra, que cromosoma llena
        aux_2 = aux_1 / ploidy;
        chrom = aux_2 * ploidy < aux_1 ? aux_2 + 1 : aux_2;
        // Recuperamos lo que sobra dentro
        aux_2 = aux_1 - aux_2 * ploidy;
        // Arreglamos
        alele = aux_2 == 0? ploidy : aux_2;

        // Recuperamos posicion de inicio real
        ini_s_i = (*select_S_i)(near_init);

        cout << raw_pos << "," << near_init << "," << sample << "," << chrom << "," << alele << "," << ini_s_i << endl;
        cout << "Sample: " << sample_name << endl;
        cout << "\t - Chrom: " << chrom << endl;
        cout << "\t - Alele: " << alele << endl;
        cout << "\t - Pos: " << (raw_pos - ini_s_i) << endl;

        result.push_back(make_pair(sampleID(sample_name, chrom, alele), raw_pos - ini_s_i));

        // Si no sobro, entonces estoy en en el resto (n_chrom * ploidy) - 1 de mi sample
        // sino 
        // // primero buscamos cuantas muestras hay hasta esta posicion
        // n_S_i_before = (unsigned int)(*rank_S_i)(raw_pos);
        // aux = n_S_i_before % (n_chrom * ploidy);

        // sample = n_S_i_before / (n_chrom * ploidy);
        // sample = (aux == 0  && sample != 0) ? sample - 1 : sample;
        
        // if(sample != last_i)
        // {
        //     getline(IDInfo_file, sample_name);
        //     last_i = sample;
        // }

        // chrom = aux / ploidy;

        // // aux = n_S_i_before % ploidy;

        // alele = n_S_i_before % ploidy;
        
        // ini_s_i = n_S_i_before == 0 ? 0 : (unsigned int) (*select_S_i)(n_S_i_before);
        
        // cout << raw_pos << "," << n_S_i_before << "," << sample << "," << chrom << "," << alele << "," << ini_s_i << endl;
        
        // cout << "Sample: " << sample_name << endl;
        // cout << "\t - Chrom: " << chrom << endl;
        // cout << "\t - Alele: " << alele << endl;
        // cout << "\t - Pos: " << (raw_pos - ini_s_i) << endl;
        // result.push_back(make_pair(sampleID(sample_name, chrom, alele), raw_pos - ini_s_i));
    }

    IDInfo_file.close();
    return result;
}
