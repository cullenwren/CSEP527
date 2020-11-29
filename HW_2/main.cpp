/**
 *   @file: main.cpp
 *
 *   @author: Cullen Billhartz
 *            Student ID: 1877747
 *
 *   @brief: Assignment #2 
 *           CSEP 527 - Computational Biology
 *           Computes Smith Waterman local alignments 
 *
 *
 *
 */


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <map> 
#include <queue>
#include <stack>
using namespace std;
/**
 * Prints full contents of a queue as a single string
 * 
 * @param q desired queue to print
 */
queue<char> printQueue(queue<char> q)
{
    while (!q.empty())
    {
        cout << " " << q.front();
        q.pop();
    }
    cout << endl;
    return q;
}

/**
* Reads a standard .fasta file and returns a string of the
 * protein
 * 
 * @param filepath desired .fasta file to read
 */
string read_fasta(string filepath)
{
    ifstream inFile(filepath);
    string line;
    string acid = "";

    while (getline(inFile, line)){
        size_t found = line.find(">");
        if (found == string::npos)
        {
            acid.append(line);
        }
    }
    return acid;
}

/**
 * Reads in a BLOSUM62 table and returns it as a map
 * where a string of two acids maps to the corresponding 
 * pair score. Could be any scoring table
 *
 * @param filepath filepath to the BLOSUM62 table 
 */
map<string, int> loadBLOSUM62(string filepath){
    map<string, int> BLOSUM62;
    map<int, char> char_idx;
    ifstream inFile(filepath);
    string line;
    while (getline(inFile, line))
    {
        size_t found = line.find("#");
        if (found == string::npos)
        {
            size_t found = line.find(" ");
            if (found == 0)
            {
                int idx = 0;
                for (int i = 0; i < line.length(); i++)
                {
                    if (line.at(i) != ' ')
                    {
                        char_idx.insert(pair<int, char>(idx, line.at(i)));
                        idx = idx + 1;
                    }
                }
            }
            else
            {
                char letter = line.at(0);
                if (letter != 'B' && letter != 'Z' && letter != 'X' && letter != '*')
                {
                    vector<string> substr;
                    string num = "";
                    for (int i = 0; i < line.length(); i++)
                    {
                        if (line.at(i) != ' ')
                            num = num + line.at(i);
                        else
                        {
                            if ((int)num.size() != 0)
                                substr.push_back(num);
                            num = "";
                        }
                    }
                    for (int i = 1; i < substr.size(); i++)
                    {
                        if (char_idx[i - 1] != 'B' && char_idx[i - 1] != 'Z' && 
                            char_idx[i - 1] != 'X' && char_idx[i - 1] != '*')
                        {
                            string acid = "";
                            acid.push_back(letter);
                            acid.push_back(char_idx[i - 1]);
                            BLOSUM62[acid] = stoi(substr[i]);
                        }
                    }
                }
            }
        }
    }
    return BLOSUM62;
}

/**
 * 
 * Implementation of the smith waterman local alignment alogorithm
 * 
 * @param series_1 string presentation of first protein
 * @param series_2 string presentation of second protein
 * @param scoring map object that accepts acid pair as key and returns pair alignment score
 * @param print boolean on whether to print score and alignment
 * @param print_table boolean on whether to print the alignment table
 */
int smith_water(string series_1, string series_2, map<string, int> scoring, bool print, bool print_table)
{

    int m = series_1.length();
    int n = series_2.length();
    int table[m][n];

    //initialize zeroth column and row
    for (int i = 0; i < m; i++){
        table[i][0] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        table[0][i] = 0;
    }

    int max_val = 0;
    int max_i = 0;
    int max_j = 0;
    //Implementation will fill row by row
    for (int i = 1; i < m; i++){
        for (int j = 1; j < n; j++){

            string acid = "";
            acid.push_back(series_1.at(i));
            acid.push_back(series_2.at(j));
            
            //score for aligning acid_1 with acid_2
            int match = scoring[acid];      
           
            int V_11 = table[i - 1][j - 1]; //diagonal adjacent score
            int V_10 = table[i - 1][j];     //row adjacent score
            int V_01 = table[i][j - 1];     //column adjacent score
            
            //Assume gap penalty of -4 
            int scores[3] = {V_11 + match, V_10 - 4, V_01 - 4};
            int max = 0;

            //Find max of values and store in table
            for (int a = 0; a < 3; a++)
            {
                if (scores[a] > max)
                    max = scores[a];
            }
            table[i][j] = max;
            
            //store if max is larger than previous max for traceback
            if (max > max_val){
                max_val = max;
                max_i = i;
                max_j = j;
            }
        }
    }
    if (print)
        cout << "Max Score of Alignment is: " << max_val << endl;
    int i = max_i;
    int j = max_j;

    if (print_table){
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::cout << table[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    }

    //Stacks of alignment traceback
    stack<char> alignment_1;    
    stack<char> alignment_2;

    //Traceback until an edge of table is hit
    while (i > 1 && j > 1){

        //Asses term of three adjacent blocks
        string acid = "";
        acid.push_back(series_1.at(i));
        acid.push_back(series_2.at(j));
        int V_11 = table[i - 1][j - 1] + scoring[acid];
        int V_10 = table[i - 1][j] - 4;
        int V_01 = table[i][j - 1] - 4;

        //Tie breaker for equal terms will go to the adjacent diagonal term
        //which is the aligned protein
        if (V_11 >= V_10 && V_11 >= V_01){
            alignment_1.push(series_1.at(i));
            alignment_2.push(series_2.at(j));
            i = i - 1;
            j = j - 1;
        }
        //Next tiebreaker goes to adjacent row
        else if (V_10 >= V_11 & V_10 >= V_01){
            alignment_1.push(series_1.at(i));
            alignment_2.push('-');
            i = i -1;
        }
        else if (V_01 >= V_11 && V_01 >= V_10){
            alignment_1.push('-');
            alignment_2.push(series_2.at(j));
            j = j - 1;
        }
    }
    if (print){

        i = 1;
        queue<char> tmp_q_1;
        queue<char> tmp_q_2;
        queue<char> tmp_q_mid;
        string acid = "";
        char a1;
        char a2;
        while (!alignment_1.empty())
        {
            if (i % 60 == 0)    //Breaking alignments into 60 char chunks
            {
                tmp_q_1 = printQueue(tmp_q_1);
                tmp_q_mid = printQueue(tmp_q_mid);
                tmp_q_2 = printQueue(tmp_q_2);
                cout << endl;
                i = i+1;
            } 
            else
            {
                a1 = alignment_1.top();
                a2 = alignment_2.top();
                tmp_q_1.push(a1);
                tmp_q_2.push(a2);
                acid = "";
                acid.push_back(a1);
                acid.push_back(a2);

                // If the acids are the same or the score is positive, 
                // Highlight the alignment in the middle row
                if (a1 == a2)
                    tmp_q_mid.push(alignment_1.top());
                else if (scoring[acid] > 0){
                    tmp_q_mid.push('+');
                }
                else
                {
                    tmp_q_mid.push(' ');
                }
                i = i + 1;
                alignment_1.pop();
                alignment_2.pop();
            }
        }
        tmp_q_1 = printQueue(tmp_q_1);
        tmp_q_mid = printQueue(tmp_q_mid);
        tmp_q_2 = printQueue(tmp_q_2);
    }

    return max_val;
}

/**
 * Random permutes a string by swapping letters 
 * 
 * @param s string to be permuted
 */ 
string permute(string s){
    int n = s.length();
    for (int i = n-1; i > 0; i--){
        int j = rand() % n;
        char tmp = s.at(j);
        s[j] = s.at(i);
        s[i] = tmp;
    }
    return s;
}

/**
 * Calcuates smith waterman and the emprical p value after 999 permutations 
 * Does not print alignment or score of each permutation alignment
 * 
 * @param sequence_1 first protein to align, all permutations will be compared to this protein 
 * @param sequence_2 second protein to align, will be permuted to p-value
 * @param scoring map to determine alignment scores
 */
void smith_water_with_permute(string sequence_1, string sequence_2, map<string, int> scoring, bool print, bool print_table, float n){
    string permutation;
    int score = smith_water(sequence_1, sequence_2, scoring, print, print_table);
    int tmp_score;
    float k = 0;
    for (int i = 0; i < n; i++){
        permutation = permute(sequence_2);
        tmp_score = smith_water(sequence_1, permutation, scoring, false, false);
        if (tmp_score > score) k = k + 1;
    } 
    float p_value = (k+1) / (n+1);
    cout.precision(5);
    cout << "P_Value for alignment is: " << scientific << p_value << endl;
}

int main()
{
    // Load scoring table
    map<string, int> BLOSUM62 = loadBLOSUM62("BLOSUM62.txt");

    // 3a) Compare deadly and ddgearlyk. Print Score matrix and empirical p-value
    //      with 999 permutations

    cout << "Aligning Sequences deadly and ddgearlyk" << endl;
    smith_water_with_permute("DEADLY", "DDGEARLYK", BLOSUM62, true, true, 999);

    // 3b) align and score each protein sequence from I(b) against the others (10 choose 2 non-identical pairs)
    string proteins[] = {"P15172.fasta", "P17542.fasta", "P10085.fasta", "P16075.fasta", "P13904.fasta", 
                        "Q90477.fasta", "Q8IU24.fasta", "P22816.fasta", "Q10574.fasta", "O95363.fasta"};
    for (int i = 0; i < 10; i++){
        for (int j = i + 1; j < 10; j++){
            cout << endl << "Aligning Sequences " << proteins[i] << " and " << proteins[j] << endl;
            string protein_1 = read_fasta(proteins[i]);
            string protein_2 = read_fasta(proteins[j]);
            //Compare P-Value on requested pairs
            if (proteins[i].compare("P15172.fasta") == 0 &&
                (proteins[j].compare("Q10574.fasta") == 0 || proteins[j].compare("O95363.fasta") == 0))
                smith_water_with_permute(protein_1, protein_2, BLOSUM62, true, false, 99999);
            else smith_water(protein_1, protein_2, BLOSUM62, true, false);
        }
    }
}