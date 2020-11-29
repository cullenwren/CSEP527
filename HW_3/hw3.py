"""
@file: main.cpp

@author: Cullen Billhartz
         Student ID: 1877747

@brief: Assignment #3 
           CSEP 527 - Computational Biology
           Executes MEME Algo 
"""
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)

"""
Reads in a source fasta file and returns list of strings
"""
def read_fasta(path):
    acids = []
    file1 = open(path, 'r') 
    Lines = file1.readlines() 
    acid = ""
    for line in Lines:
        if '>' in line:
            acids.append(acid)
            acid = ""
        else:
            acid = acid + line.replace('\n', '')

    acids.append(acid)

    return acids

"""
Constructs the count matrix from a list of sequences 
"""
def makeCountMatrix(sequences):
    k = len(sequences[0])
    countMatrix = np.zeros(shape=(4,k))

    for sequence in sequences:
        for i, letter in enumerate(sequence):
            if letter == 'A':
                countMatrix[0][i] += 1
            elif letter == 'C':
                countMatrix[1][i] += 1
            elif letter == 'G':
                countMatrix[2][i] += 1
            elif letter == 'T':
                countMatrix[3][i] += 1
    return countMatrix

"""
Adds a pseudoVector to a countMatrix
"""
def addPseudo(countMatrix, pseudoVec):
    for i, row in enumerate(countMatrix):
        for j, col in enumerate(row):
            countMatrix[i, j] = col + pseudoVec[i] 

    return countMatrix

"""
Constructs the frequency matrix from a count matrix
"""
def makeFrequencyMatrix(countMatrix):

    frequencyMatrix = np.zeros(shape=countMatrix.shape)
    
    for col in range(countMatrix.shape[1]):
        frequencyMatrix[:, col] = countMatrix[:, col] / sum(countMatrix[:, col])

    return frequencyMatrix

"""
Calculates the total entropy of a frequency matrix relative to a 
background vector
"""
def entropy(frequencyMatrix, background):
    entr = 0
    for i, row in enumerate(frequencyMatrix):
        for freq in row:
            entr += freq * math.log2(freq / background[i])
    
    return entr

"""
Constructs the weight matrix model from a frequency matrix and 
background vector.

Also returns the entropy of the matrix
"""
def makeWMM(frequencyMatrix, background):
    WMM = np.zeros(shape=frequencyMatrix.shape)

    for i, row in enumerate(frequencyMatrix):
        for j, freq in enumerate(row):
            WMM[i,j] = math.log2(freq / background[i]) 

    entr = entropy(frequencyMatrix, background)
    return entr, WMM

"""
Scans and scores a list of sequences with a given weight
matrix model
"""
def scanWMM(WMM, sequences):
    scores = []
    seq_len = WMM.shape[1]
    for i, sequence in enumerate(sequences):
        scores.append(np.zeros(shape=(len(sequence) - seq_len + 1)))

        for j in range(len(sequence) - seq_len+1):
            window = sequence[j:j + seq_len]
            score = 0
            for k, letter in enumerate(window):
                if letter.upper() == 'A':
                    score += WMM[0][k]
                elif letter.upper() == 'C':
                    score += WMM[1][k]
                elif letter.upper() == 'G':
                    score += WMM[2][k]
                elif letter.upper() == 'T':
                    score += WMM[3][k]
            scores[i][j] = score
    return scores

"""
Estimates the likilhood of a motif starting at a position 
from a list of sequences based on the weight matrix model
"""
def Estep(WMM, sequences):    
    scores = scanWMM(WMM, sequences)
    EMatrix = []
    for i, row in enumerate(scores):
        E_row = np.zeros(shape=len(row))
        total_score = 0
        for score in row:
            total_score = total_score + 2**score
        for j, score in enumerate(row):
            E_row[j] = 2**score / total_score
        EMatrix.append(E_row)
    return EMatrix

"""
Constructs a new count matrix (and then weight matrix) based 
on the probability associated with the new motif starting at a position
"""
def Mstep(EMatrix, sequences, WMM):
    seq_len = WMM.shape[1]
    countMatrix = np.zeros(shape = WMM.shape)
    for i, sequence in enumerate(sequences):
        for j in range(len(sequence) - seq_len+1):
            window = [sequence[j:j + seq_len]]
            prob = EMatrix[i][j]
            countMatrix_tmp = (prob * makeCountMatrix(window))
            countMatrix = countMatrix + countMatrix_tmp

    #print(countMatrix)
    countMatrix_new = addPseudo(countMatrix, [0.25, 0.25, 0.25, 0.25])
    frequencyMatrix = makeFrequencyMatrix(countMatrix_new)
    entr, WMM = makeWMM(frequencyMatrix, [0.25, 0.25, 0.25, 0.25])
    return entr, WMM, frequencyMatrix

"""
Initalizes a list of weight matrix models of length k from the
sequence provided
"""
def initialize(sequence, k):
    background = [0.25, 0.25, 0.25, 0.25]
    WMM_list = []
    i = 0
    while i <= len(sequence) - k:
        window = sequence[i:i+k]
        countMatrix = makeCountMatrix([window])
        pseudoVec = np.full((4), 0.15 / 0.85 / 4)
        countMatrix = addPseudo(countMatrix, pseudoVec)
        frequencyMatrix = makeFrequencyMatrix(countMatrix)
        entr, WMM_tmp = makeWMM(frequencyMatrix, background)
        WMM_list.append(WMM_tmp)
        i = int(i + k/2)
        
    return WMM_list
        
if __name__ == "__main__":
    #Read in train and test sequences
    acids_train = read_fasta("hw3-train.fasta")[1:]
    acids_test = read_fasta("hw3-test.fasta")[1:]
    #Initialize on the training sequence
    WMMs = initialize(acids_train[0], 10)

    #Begin by training with 3 iterations
    entr_list = []
    freqMatrix_list = []
    WMM_list = []
    for WMM in WMMs: 
        for i in range(3):
            EMatrix = Estep(WMM, acids_train)
            entr, WMM, frequencyMatrix = Mstep(EMatrix, acids_train, WMM)
        
        WMM_list.append(WMM)
        entr_list.append(entr)
        freqMatrix_list.append(frequencyMatrix)

    #Store indexs of Best, median, and worst Entropies
    A = np.argmax(entr_list)
    B = np.argsort(entr_list)[len(entr_list)//2]
    C = np.argmin(entr_list)

    WMM_ABCD = [WMM_list[A]]
    WMM_ABCD.append(WMM_list[B])
    WMM_ABCD.append(WMM_list[C])
    freq_ABCD = [freqMatrix_list[A]]
    freq_ABCD.append(freqMatrix_list[B])
    freq_ABCD.append(freqMatrix_list[C])

    #Data structures to store calculated data
    entropies = pd.DataFrame(columns = range(10))
    WMM_list = []
    entr_list = []
    freqMatrix_list = []

    #Execute EM 10 times on training data
    for s, WMM in enumerate(WMMs):
        for i in range(10):
            EMatrix = Estep(WMM, acids_train)
            entr, WMM, frequencyMatrix = Mstep(EMatrix, acids_train, WMM)
            entropies.loc[s, i] = entr

        entr_list.append(entr)
        freqMatrix_list.append(frequencyMatrix)
        WMM_list.append(WMM)

    #Print entropies of each iteration 
    print(entropies)
    
    #Find highest entropy model
    D = np.argmax(entr_list)
    WMM_ABCD.append(WMM_list[D])
    freq_ABCD.append(freqMatrix_list[D])

    #Print Frequency matricies
    freq_A = pd.DataFrame(freq_ABCD[0], index = ['A', 'C', 'G', 'T'], columns=range(10))
    freq_B = pd.DataFrame(freq_ABCD[1], index = ['A', 'C', 'G', 'T'], columns=range(10))
    freq_C = pd.DataFrame(freq_ABCD[2], index = ['A', 'C', 'G', 'T'], columns=range(10))
    freq_D = pd.DataFrame(freq_ABCD[3], index = ['A', 'C', 'G', 'T'], columns=range(10))
    print('Frequency Matrix for Matrix A, index {} \n'.format(A))
    print(freq_A.round(3))
    print('\n')
    print('Frequency Matrix for Matrix B, index {} \n'.format(B))
    print(freq_B.round(3))
    print('\n')
    print('Frequency Matrix for Matrix C, index {} \n'.format(C))
    print(freq_C.round(3))
    print('\n')
    print('Frequency Matrix for Matrix D, index {} \n'.format(D))
    print(freq_D.round(3))
    print('\n')

    #Evaluate the four models against the training acids
    scores_list = []
    starting_pos = []
    legends = ['WMM_A', 'WMM_B', 'WMM_C', 'WMM_D']
    #Determine the most common predicted starting location
    for model, wmm in enumerate(WMM_ABCD): 
        scores = scanWMM(wmm, acids_test)
        scores_list.append(scores)
        hist = np.zeros(shape = (len(scores[0])))
        for i, row in enumerate(scores):
            idx = np.argmax(row)
            hist[idx] += 1
        starting_pos.append(np.argmax(hist))
        plt.bar(range(len(scores[0])), hist)
        plt.title('Most common start position: {}'.format(np.argmax(hist)))
        plt.xlabel('Index of Estimated Starting Position')
        plt.ylabel('Frequency')
        plt.legend([legends[model]])
        plt.savefig('WMM_{}.png'.format(model))
        plt.close()

    #Construct the ROC and calcuate AUC 
    point = [0, 1, 0]
    for model, score_mat in enumerate(scores_list):
        pairs = []
        counts = pd.DataFrame()
        for i, seq in enumerate(score_mat):
            for j, score in enumerate(seq):
                counts.loc[j, i] = score
                
        sort_scores = sorted(counts.values.flatten())
        
        AUC = 0
        for score in sort_scores:
            tmp = counts > score
            motif = tmp.loc[starting_pos[model], :]
            not_motif = tmp.drop(starting_pos[model] ,axis = 0)
            
            tp = np.sum(motif.values.flatten())
            tn = np.sum(~not_motif.values.flatten())
            fp = np.sum(not_motif.values.flatten())
            fn = np.sum(~motif.values.flatten())
            tpr = tp / (tp + fn)
            fpr = fp / (fp + tn)
            if (model == 2) and (tpr == 1) and (fpr < point[1]):
                #print(point)
                point = [tpr, fpr, score]
            pairs.append([tpr, fpr])
        
        pairs_np = np.array(pairs)
        AUC = np.trapz(pairs_np[:, 0], pairs_np[:, 1]) * -1
        plt.plot(pairs_np[:,1], pairs_np[:,0])
        print('AUC for Model {}: {}'.format(model, AUC))

    plt.scatter(point[1], point[0])
    plt.legend(['WMM_A', 'WMM_B','WMM_C', 'WMM_D'])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig('ROCs.png')
    plt.close()

    print('Point of interest at ({}, {}) with score of {}'.format(point[1], point[0], point[2]))
