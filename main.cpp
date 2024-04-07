#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <stdlib.h>

#include "src/util.hpp"
#include "src/hierarchical_overlap_graph.hpp"

using namespace std;

void random_replace_ACGT(string &in, const string &search_value);

int main(int argc, char **argv)
{
    // Print memory occupation at the very start of the program
    memory_usage_MB("Memory occupation at the very start: ");

    // Set random seed for all random operations to make the execution replicable
    srand(1);

    // Where the reads are stored
    vector<string> reads;

    // Parse all files in input as argument
    for (int i = 0; i < argc; i++)
    {
        string filename{argv[i]};
        ifstream input_reads(filename);

        /**
         * fastq file format:
         * file contains records of reads with their metadata. Every read
         * occupies 4 rows. In the first there is an ID. The SECOND contains
         * the sequence of nucleotides. The thirds contains a '+' and some
         * optional desctioption. The fourth contains the quality of the sequence.
         *
         *
         * Example:
         *
         * @SEQ_ID
         * GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
         * +
         * !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
         * ...
         *
         */
        if (input_reads.is_open() && ends_with(filename, string(".fastq")))
        {
            string read;
            int row_number = 0;
            while (getline(input_reads, read))
            {
                if ((row_number++) % 4 == 1)
                {
                    read = trim(read, "N");
                    random_replace_ACGT(read, "N");
                    reads.push_back(read);
                }
            }
        }
        /**
         * text file format:
         * file contains records of reads. Every line contains a
         * read sequence. Nucleotide representation can be lower
         * case or upper case.
         *
         * Example:
         *
         * GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
         * cagtatcgatcaaatagtaacgaagtaacgataacgatcaaat
         * GgtTcAAAGCaAtcGaTCAAatAGtaAatcCaTTTG
         * AGTATCGATCAAATAGTAAAGCAGTATCGATCAAATCCATTTGTCAACTCAC
         * ...
         *
         */
        else if (input_reads.is_open())
        {
            string read;
            while (getline(input_reads, read))
            {
                to_upper(read);

                if (read.find_first_of("ACGTN") == 0)
                {
                    read = trim(read, "N");
                    random_replace_ACGT(read, "N");
                    reads.push_back(read);
                }
            }
        }
        input_reads.close();
    }

    cout << "Number of reads: " << reads.size() << std::endl;  // Print number of imported reads
    memory_usage_MB("Memory occupation after reads import: "); // Print memory footprint after reads import

    auto start_elab = std::chrono::high_resolution_clock::now();
    HOG hog(reads);
    auto end_elab = std::chrono::high_resolution_clock::now();

    // t.print();

    cout << "Overall elaboration time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_elab - start_elab).count() << "ms" << endl;
    cout << "Trie size: " << hog.size() << endl;
    memory_usage_MB("Memory occupation at the end: ");

    return 0;
}

/* Useful to replace all search_value occurrences of in_str with a random from {'A', 'C', 'G', 'T'} */
void random_replace_ACGT(string &in_str, const string &search_value)
{
    size_t start_pos;
    while ((start_pos = in_str.find_first_of(search_value)) != string::npos)
    {
        int random_number = rand() % 4;

        switch (random_number)
        {
        case 0:
            in_str.replace(start_pos, 1, "A");
            break;
        case 1:
            in_str.replace(start_pos, 1, "C");
            break;
        case 2:
            in_str.replace(start_pos, 1, "G");
            break;
        case 3:
            in_str.replace(start_pos, 1, "T");
        }
    }
}
