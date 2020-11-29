#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <queue>
#include <stack>
#include <math.h>
#include <algorithm>
#include <chrono>
using namespace std;
using namespace std::chrono;
struct Traceback
{
    stack<int> states;
    stack<array<int, 3>> hits;
};

struct StateParameters
{
    vector<map<char, double>> emissions;
    vector<vector<double>> transitions;
};

/**
* Reads a standard .fasta file and returns a string of the
 * protein
 * 
 * @param filepath desired .fasta file to read
 */
vector<string> read_fasta(string filepath)
{
    ifstream inFile(filepath);
    string line;
    string acid = "";
    vector<string> acid_vector;

    while (getline(inFile, line))
    {
        size_t found = line.find(">");
        if (found == string::npos)
        {
            acid.append(line);
        }
        else
        {
            acid_vector.push_back(acid);
            acid = "";
        }
    }
    acid_vector.push_back(acid);
    return acid_vector;
}

vector<array<double, 2>> calcScores(StateParameters params, array<double, 2> score, string genome){
    vector<map<char, double>> emissions = params.emissions;
    vector<vector<double>> transitions = params.transitions;
    vector<array<double, 2>> scores;
    double score_1;
    double score_2;
    char c;
    for (int i = 0; i < genome.length(); i++)
    {
        c = genome.at(i);
        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T'))
        {
            c = 'T';
        }
        if (i == 0)
        {
            score[0] = score[0] + log(emissions[0][c]);
            score[1] = score[1] + log(emissions[1][c]);
        }
        else
        {
            score_1 = scores.back()[0] + log(transitions[0][0]) + log(emissions[0][c]);
            score_2 = scores.back()[1] + log(transitions[1][0]) + log(emissions[0][c]);

            if (score_1 > score_2)
                score[0] = score_1;
            else
                score[0] = score_2;

            score_1 = scores.back()[0] + log(transitions[0][1]) + log(emissions[1][c]);
            score_2 = scores.back()[1] + log(transitions[1][1]) + log(emissions[1][c]);
            if (score_1 > score_2)
                score[1] = score_1;
            else
                score[1] = score_2;
        }
        scores.push_back(score);

    }
    return scores;
}

Traceback getTraceBack(vector<array<double, 2>> scores, string genome, vector<vector<double>> transitions)
{
    Traceback trace;
    array<double, 2> score;
    int state;
    char c;
    stack<int> path;
    score = scores.back(); //Determining ending state
    scores.pop_back();
    if (score[0] > score[1])
        state = 0;
    else
        state = 1;

    path.push(state);

    cout << "Overall Viterbi path log probability: " << score[state] << endl;

    int prev_state = state;
    int count = 0;
    int length = 0;
    bool x = false;
    bool y = false;
    array<int, 3> hit = {0, 0, 0};
    stack<array<int, 3>> hits;
    // Traceback on rest of states
    for (int i = scores.size(); i > 0; i--)
    {
        score = scores.back();
        scores.pop_back();

        if (prev_state == 1) // if previous state was state 1
        {
            //transition from 1 -> 1 more likely than 0 --> 1
            x = (score[1] + log(transitions[1][1]) >= score[0] + log(transitions[0][1]));
        }
        else //previous state was 0
        {
            //transition from 1 -> 0 more likely than 0 --> 0
            x = (score[1] + log(transitions[1][0]) >= score[0] + log(transitions[0][0]));
        }

        if (!(x)) // Transition was more likely from 0 -> 1/0
        {
            //Next State: 0
            state = 0;
            path.push(state);
            if (prev_state == 1)
            {
                // add to count and reset length var
                count += 1;
                prev_state = 0;
                hits.push({i + 1, i + length, length});
                length = 0;
            }
        }
        else
        {
            //Next State: 1
            state = 1;
            prev_state = 1;
            length += 1;
            path.push(state);
        }
    }
    if (prev_state == 1)
    {
        // add to count and reset length var
        count += 1;
        prev_state = 0;
        length = 0;
    }
    trace.hits = hits;
    trace.states = path;
    return trace;
}

void printHits(Traceback trace, int k){
    array<int, 3> hit;
    if (k == 0) k = trace.hits.size();
    for (int i = 0; i < k; i++)
    {
        if (!(trace.hits.empty())){
            hit = trace.hits.top();
            cout << hit[0] << "  " << hit[1] << "  Length: " << hit[2] << endl;
            trace.hits.pop();
        }
    }
    cout << endl;
}

void printParams(StateParameters params){
    vector<map<char,double>> emissions = params.emissions;
    vector<vector<double>> transitions = params.transitions;
    cout.precision(10);
    cout << "Emissions  " << "A " << "C " << "G " << "T " << endl;
    cout << "State 1    " << emissions[0]['A'] << " " << emissions[0]['C'] << " " << emissions[0]['G'] << " " << emissions[0]['T'] << endl;
    cout << "State 2    " << emissions[1]['A'] << " " << emissions[1]['C'] << " " << emissions[1]['G'] << " " << emissions[1]['T'] << endl;
    cout << endl;
    cout << "Transitions  " << "State 1 " << "State 2 " << endl;
    cout << "State 1    " << transitions[0][0] << " " << transitions[0][1] << endl;
    cout << "State 2    " << transitions[1][0] << " " << transitions[1][1] << endl;
    cout << endl;
}

StateParameters getNewParams(stack<int> states, string genome)
{
    StateParameters newParams;
    vector<vector<double>> transitions{{0, 0},{0, 0}};
    map<char, double> emission_1;
    emission_1['A'] = 0;
    emission_1['C'] = 0;
    emission_1['G'] = 0;
    emission_1['T'] = 0;

    map<char, double> emission_2;
    emission_2['A'] = 0;
    emission_2['C'] = 0;
    emission_2['G'] = 0;
    emission_2['T'] = 0;

    vector<map<char, double>> emissions = { emission_1,
                                            emission_2 };
    int prev_state = 0;
    int state = 0;
    char c;
    float i_t = 0; //State 0 transition count
    float j_t = 0; //State 1 transition count
    float i_e = 0; //State 0 emission count
    float j_e = 0; //State 1 emission count
    while (!(states.empty())){

        c = genome.at(i_e+j_e);

        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T'))
        {
            c = 'T';
        }

        state = states.top();
        if (state == 0 && prev_state == 0){
            emissions[0][c] += 1;
            transitions[0][0] += 1;
            i_t += 1;
            i_e += 1;
        }
        else if (state == 0 && prev_state == 1){
            emissions[0][c] += 1;
            transitions[1][0] += 1;
            j_t += 1;
            i_e += 1;
        }
        else if (state == 1 && prev_state == 0)
        {
            emissions[1][c] += 1;
            transitions[0][1] += 1;
            i_t += 1;
            j_e += 1;
        }
        else if (state == 1 && prev_state == 1)
        {
            emissions[1][c] += 1;
            transitions[1][1] += 1;
            j_t += 1;
            j_e += 1;
        }
        prev_state = state;
        states.pop();
    }
    transitions[0][0] = transitions[0][0] / i_t;
    transitions[0][1] = transitions[0][1] / i_t;
    transitions[1][0] = transitions[1][0] / j_t;
    transitions[1][1] = transitions[1][1] / j_t;
    emissions[0]['A'] = emissions[0]['A'] / i_e;
    emissions[0]['C'] = emissions[0]['C'] / i_e;
    emissions[0]['G'] = emissions[0]['G'] / i_e;
    emissions[0]['T'] = emissions[0]['T'] / i_e;
    emissions[1]['A'] = emissions[1]['A'] / j_e;
    emissions[1]['C'] = emissions[1]['C'] / j_e;
    emissions[1]['G'] = emissions[1]['G'] / j_e;
    emissions[1]['T'] = emissions[1]['T'] / j_e;
    newParams.transitions = transitions;
    newParams.emissions = emissions;
    return newParams;
}

int main()
{
    cout.precision(10);
    string genome = read_fasta("GCF_000091665.1_ASM9166v1_genomic.fna")[1];
    cout << genome.size() << endl;
    //string genome = read_fasta("test.fasta")[1];
    vector<map<char,double>> emissions;
    vector<vector<double>> transitions;
    array<double, 2> score; 
    map<char, double> emission_1;
    emission_1['A'] = .25;
    emission_1['C'] = .25;
    emission_1['G'] = .25;
    emission_1['T'] = .25;

    map<char, double> emission_2;
    emission_2['A'] = .20;
    emission_2['C'] = .30;
    emission_2['G'] = .30;
    emission_2['T'] = .20;
    
    emissions = {emission_1, emission_2};
    transitions = {{0.9999, 0.0001}, {0.01, 0.99}};
    score = {{log(0.9999), log(0.0001)}};

    StateParameters params;
    params.emissions = emissions;
    params.transitions = transitions;
    printParams(params);
    Traceback trace;
    vector<array<double, 2>> scores;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 9; i++)
    {
        cout << "Iteration: " << i+1 << endl;
        scores = calcScores(params, score, genome);
        trace = getTraceBack(scores, genome, transitions);
        printHits(trace, 5);

        params = getNewParams(trace.states, genome);
        printParams(params);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<microseconds>(stop - start);

    cout << "Iteration: " << 10 << endl;
    scores = calcScores(params, score, genome);
    trace = getTraceBack(scores, genome, transitions);
    printHits(trace, 0);
    params = getNewParams(trace.states, genome);
    printParams(params);
    
}