#include "VCFParsingInterpreter.h"

using namespace std;
using namespace sdsl;

#define INPUT_BINARY_FILE ios::in | ios::binary
#define OUTPUT_BINARY_FILE ios::out | ios::binary

VCFParsingInterpreter::VCFParsingInterpreter()
{
}

void VCFParsingInterpreter::InitializeFromPreloadedFile(char *folder_path)
{
    // Recover index
    string Destination_folder_name(folder_path);
    string Destination_aux = Destination_folder_name + "vcf-rlz-index";
    Index = new RelzIndexReference();
    Index->load(Destination_aux);

    // Recover ID name data
    IDInfo_file_path = Destination_folder_name + ".idinfo";

    // Recover numeric data
    fstream src(Destination_aux + ".data", INPUT_BINARY_FILE);
    src.read((char*)&n_chromosomes, sizeof(int));
    src.read((char*)&n_S_i, sizeof(int));
    src.read((char*)&S_size, sizeof(int));
    src.close();

    // Recover sdsl structurescd ..
    bit_vector_S_i = new sd_vector<>();
    load_from_file((*bit_vector_S_i), Destination_aux + ".sdv");

    rank_S_i = new sd_vector<>::rank_1_type(bit_vector_S_i);
    select_S_i = new sd_vector<>::select_1_type(bit_vector_S_i);
}

void VCFParsingInterpreter::InitializeFromParsing(char *destination_path)
{
    // Something like .../VCF_files/
    string Destination_path_aux(destination_path);
    InitializeFromParsing(Destination_path);
}

void VCFParsingInterpreter::InitializeFromParsing(string destination_path)
{
    // Something like .../VCF_files/
    string Destination_path_aux = destination_path;
    Destination_path = Destination_path_aux;

    // TODO: Assert if ./Tmp is a valid folder
    IDInfo_file_path = Destination_path + IDInfo_file_subpath;

    cout << "[RLZ] Reading metadata" << endl;
    NanoTimer timer;
    ProcessReference();
    ProcessMetaParsing();
    cout << "[RLZ]      Metadata readed in " << timer.getMilisec() << " ms" << endl;
    cout << "[RLZ] Building factors from phrases" << endl;
    timer.reset();
    BuildFactors();
    // for (pair<unsigned int, unsigned int> factor : factors)
    // {
    //     cout << "(" << factor.first << "," << factor.second << ")" << endl;
    // }
    cout << "[RLZ] Building index" << endl;
    char* reference = GetReference();
    timer.reset();
    Index = new RelzIndexReference(factors, reference, S_size, reference, Reference_len);
    cout << "[RLZ]      Index finished in " << timer.getMilisec() << " ms" << endl;
    FreeCacheVariables();
}

void VCFParsingInterpreter::FreeCacheVariables()
{

}

void VCFParsingInterpreter::ProcessReference()
{
    // Open required files
    Reader.open(Destination_path + MetaReference_file_subpath, INPUT_BINARY_FILE);

    // The first structure contains summary info
    // cout << "Lee la estructura" << endl;
    metareference metareference_buffer;
    Reader.read((char *)&metareference_buffer, sizeof(metareference));
    Reference_len = metareference_buffer.n_bases();
    int n_References = metareference_buffer.rel_pos();

    // cout << "Largo: " << Reference_len << "|Referencias: " << n_References << endl;

    // Collect all metadata
    for (int i = 0; i < n_References; i++)
    {
        Reader.read((char *)&metareference_buffer, sizeof(metareference));
        dict_metareference[metareference_buffer.ID()] = metareference_buffer;
        // if (i <= 4)
        // {
        //     cout << "Ref: " << metareference_buffer.ID() << ", n_bases: " << metareference_buffer.n_bases()
        //          << ", rel_pos: " << metareference_buffer.rel_pos() << endl;
        // }
    }
    Reader.clear();
    Reader.close();
}

void VCFParsingInterpreter::ProcessMetaParsing()
{
    Reader.open(Destination_path + MetaParsing_file_subpath, INPUT_BINARY_FILE);

    metainfo metainfo_buffer;

    Reader.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_Phrases = metainfo_buffer.n_phrases();
    Reader.read((char *)&metainfo_buffer, sizeof(metainfo));
    n_chromosomes = metainfo_buffer.n_phrases();
    Reader.close();
    Reader.close();
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
    // 0-index to 1-index
    rel_pos_chrom = dict_metareference[Phrase.chrom()].rel_pos() + 1;
    // cout << "1-rel_pos_chrom " << rel_pos_chrom << endl;
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom,     // Init position of current reference chromosome
                                                    Phrase.pos());      // Copy until before the actual position
    // cout << "(" << aux.first << "," << aux.second << ") : init" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateAleleInterFactor(phrase Phrase1, phrase Phrase2)
{
    // 0-index to 1-index
    rel_pos_chrom = dict_metareference[Phrase1.chrom()].rel_pos() + 1;
    // cout << "2-rel_pos_chrom " << rel_pos_chrom << endl;
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom + Phrase1.pos() + (Phrase1.len() - 1),          // Start after the discarded slice of Phrase1
                                        Phrase2.pos() - ((Phrase1.pos() + 1) + (Phrase1.len() - 1))); // Fill until reach before pos of Phrase2
    // cout << "(" << aux.first << "," << aux.second << ") : inter" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateAleleEndFactor(phrase Phrase)
{
    // 0-index to 1-index
    rel_pos_chrom = dict_metareference[Phrase.chrom()].rel_pos() + 1;
    // cout << "3-rel_pos_chrom " << rel_pos_chrom << endl;
    pair<unsigned int, unsigned int> aux = make_pair(rel_pos_chrom + Phrase.pos() + (Phrase.len() - 1),  // Start after the discarded slice of Phrase
                     dict_metareference[Phrase.chrom()].n_bases() - ((Phrase.pos() + 1) + (Phrase.len() - 1) + 1)); // Fill until reach the end of the current chromosome-alele
    // cout << "(" << aux.first << "," << aux.second << ") : end" << endl;
    return aux;
}

pair<unsigned int, unsigned int> VCFParsingInterpreter::CalculateFullAleleFactor(ll chromosome)
{
    // 0-index to 1-index
    rel_pos_chrom = dict_metareference[chromosome].rel_pos() + 1;
    // cout << "4-rel_pos_chrom " << rel_pos_chrom << endl;
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
            S_i_pos.push_back(S_size );
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
            S_i_pos.push_back(S_size );
            HigienicFactorPushBack(aux_factor);
            // Then induce the init factor for the current alele
            aux_factor = CalculateAleleInitFactor(curr_phrase);
            S_i_pos.push_back(S_size );
            HigienicFactorPushBack(aux_factor);
        }
        else if (curr_is_dummy)
        {
            // We need to induce an end factor, closing the last alele
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);
            // We need to induce the full factor for the curr and last alele
            aux_factor = CalculateFullAleleFactor(last_phrase.chrom());
            S_i_pos.push_back(S_size );
            HigienicFactorPushBack(aux_factor);
        }
        else
        {
            // We need to induce an end factor, closing the last alele
            aux_factor = CalculateAleleEndFactor(last_phrase);
            HigienicFactorPushBack(aux_factor);
            // We need to induce the init factor for the current alele
            aux_factor = CalculateAleleInitFactor(curr_phrase);
            S_i_pos.push_back(S_size );
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
            S_i_pos.push_back(S_size );
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
        S_i_pos.push_back(S_size );
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
            S_i_pos.push_back(S_size );
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

    n_S_i = S_i_pos.size();

    bit_vector_S_i = new sd_vector<>(b);
    rank_S_i = new sd_vector<>::rank_1_type(bit_vector_S_i);
    select_S_i = new sd_vector<>::select_1_type(bit_vector_S_i);
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
    Reader.open(Destination_path + Phrases_file_subpath, INPUT_BINARY_FILE);
    // Read first dummy phrase
    Reader.read((char *)&curr_phrase, sizeof(phrase));

    // Initialize indentifier variants
    phrase last_phrase = curr_phrase;
    rel_pos_chrom = dict_metareference[curr_phrase.chrom()].rel_pos();

    for (ll i = 1; i < n_Phrases; i++)
    {
        // Read the next phrase (Read the first actual phrase)
        Reader.read((char *)&curr_phrase, sizeof(phrase));
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

    Reader.clear();
    Reader.close();

    int actual_indv = 0, actual_chrom = 0, actual_alele = 0;

    BuildReconstructionStructures();
}

char* VCFParsingInterpreter::GetReference()
{
    Reader.open(Destination_path + Reference_file_subpath, ios::in);

    string reference_aux;
    if (Reader.is_open())
    {
        getline(Reader, reference_aux);
    }

    Reader.clear();
    Reader.close();

    const char *reference_const = reference_aux.c_str();
    char *reference = new char[Reference_len + 1];
    strcpy(reference, reference_const);

    return reference;
}

vector<pair<sampleID, unsigned int>> VCFParsingInterpreter::FindSnippet(string snippet, bool show)
{
    vector<pair<sampleID, unsigned int>> result;
    unsigned int n_chrom = (unsigned int) n_chromosomes;
    unsigned int near_init, sample, chrom, alele, last_i, ini_s_i, next_s_i, aux_1, aux_2;
    string sample_name;

    vector<unsigned int> raw_positions;
    // cout << 1 << endl;
    Index->findTimes(snippet, raw_positions);
    // cout << 2 << endl;

    if (raw_positions.empty())
    {
        return result;
    } else {
        cout << "Occ: "<< raw_positions.size() << endl;
    }

    Reader.open(IDInfo_file_path);

    sort(raw_positions.begin(), raw_positions.end());

    // Initialize ID
    last_i = 1;
    getline(Reader, sample_name);
    cout << "Sample: " << sample_name << endl;

    // TODO
    int max_query = n_S_i;

    for(unsigned int raw_pos : raw_positions)
    {
        // cout << " Raw_pos: " << raw_pos << endl;
        near_init = (*rank_S_i)(raw_pos);
        // Recuperamos posicion de inicio real
        ini_s_i = (*select_S_i)(near_init);
        if(near_init < max_query)
        {
            next_s_i = (*select_S_i)(near_init + 1);
        }
        else
        {
            next_s_i = S_size;
        }

        // cout << near_init << "," << ini_s_i << "," << raw_pos << "," << snippet.size() << "," << next_s_i << endl;

        // Si mi snippet no esta completamente dentro de mi seccion, entonces se descarta
        if (raw_pos + snippet.size() - 1 > next_s_i)
        {
            // It means a pattern between S_i
            // cout << "   Discarded! " << next_s_i << endl;
            continue;
        }

        // Cuantos Indv lleno hasta este punto?
        aux_1 = near_init / (n_chrom * ploidy);
        // Si lleno justo, entonces me quedo en aux, sino paso al sgte
        sample = aux_1 * (n_chrom * ploidy) < near_init ? aux_1 + 1 : aux_1;

        // Update sample name if required
        if(sample != last_i)
        {
            getline(Reader, sample_name);
            last_i = sample;
            if (show) cout << "Sample: " << sample_name << endl;
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

        // cout << raw_pos << "," << near_init << "," << sample << "," << chrom << "," << alele << "," << ini_s_i << endl;
        if (show)
        {
            cout << " - Chrom: " << chrom << endl;
            cout << "   - Alele: " << alele << endl;
            cout << "     - Pos: " << (raw_pos - ini_s_i) << endl;
        }

        result.push_back(make_pair(sampleID(sample_name, chrom, alele), raw_pos - ini_s_i));
    }

    Reader.close();

    return result;
}

void VCFParsingInterpreter::SaveInterpreter()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y_%H:%M:%S/",timeinfo);
    string current_time(buffer);
    string Destination_folder_name = Destination_path + "VCF_RLZ_" + current_time;

    if(mkdir(Destination_folder_name.c_str(), 0777) == -1)
    {
        cerr << "[RLZ] Save Error : " << strerror(errno) << endl;
        return;
    }
    // Save index
    string Destination_aux = Destination_folder_name + "vcf-rlz-index";
    Index->save(Destination_aux);

    // Save ID name data
    fstream  src(IDInfo_file_path, INPUT_BINARY_FILE);
    fstream  dst(Destination_aux + ".idinfo", OUTPUT_BINARY_FILE);
    dst << src.rdbuf();
    src.close();
    dst.close();

    // Save numeric data
    dst.open(Destination_aux + ".data", OUTPUT_BINARY_FILE);
    dst.write((char*)&n_chromosomes, sizeof(int));
    dst.write((char*)&n_S_i, sizeof(ll));
    dst.write((char*)&S_size, sizeof(ll));
    dst.close();

    // Save sdsl structures
    store_to_file((*bit_vector_S_i), Destination_aux + ".sdv");

    // Eliminar la carpeta Tmp
    string command = "rm -r " + Destination_path + "Tmp/";
    if(system(command.c_str()) == -1)
    {
        cerr << "[RLZ] Delete Error : Please try to delete manually the following folder to free space: " << Destination_path + "Tmp/" << endl;
        cerr << "\t Error caused by: " << strerror(errno) << endl;
    }
}

double VCFParsingInterpreter::GetSize()
{
    double total_size = 0;

    total_size += Index->getSize();

    total_size += size_in_bytes(*bit_vector_S_i);

    return total_size;
}
