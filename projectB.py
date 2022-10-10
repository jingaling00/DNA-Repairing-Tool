# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:10:49 2022

@author: jingy
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:05:08 2022

@author: jingy
"""

codons_dict = {"UUU":"Phe", "UUC":"Phe", "UUA":"Leu", "UUG":"Leu",
          "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
          "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
          "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
          "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
          "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
          "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
          "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
          "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
          "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
          "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
          "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
          "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
          "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
          "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
          "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}
    
def readFile(fileName):
    """
    Reads a text file.
    
    Parameters
    ----------
    fileName : str
        File path to read.

    Returns
    -------
    str
        Text from file.
    """
    with open(fileName,'r') as dnaFile:
        dna = "".join(dnaFile.readlines()).strip()
    return dna
    
    
def writeFile(fileName,text):
    """
    Writes a text file.

    Parameters
    ----------
    fileName : str
        File path to write.
    text : str
        Text to write.

    Returns
    -------
    None.

    """
    with open(fileName,'w') as textFile:
        textFile.write(text)

def transcribe(dna): # perform 5' to 3' transcription: convert Thymine to Uracil
    mrna = []
    for base in dna:
        if base == 'T':
            base = 'U'
        mrna.append(base)
    mrna = ''.join(mrna)
    return mrna

def translate(mrna): # perform mRNA translation according to codons
    n = len(mrna)
    parsed_mrna = [mrna[i:i+3] for i in range(0,n-1,3)]
    amino_acids = []
    for codon in parsed_mrna:
        if len(codon) == 3:
            codon = codons_dict[codon]
            amino_acids.append(codon)
    aa_sequence = ' '.join(amino_acids)
    return aa_sequence

def synonymous(sub,ref): # evaluate if two DNA sequences encode the same amino acids
    transcribe_sub = transcribe(sub)
    protein_sub = translate(transcribe_sub)
    
    transcribe_ref = transcribe(ref)
    protein_ref = translate(transcribe_ref)
    
    if protein_sub == protein_ref:
        return True
    else:
        return False

def delete(dna,i): # perform nucleotide deletion at index i
    dna = dna[:i] + dna[i+1:]
    return dna # delete dna at index i

def insert(dna,i,base): # perform nucleotide insertion at index i
    dna = dna[:i] + base + dna[i:]
    return dna

def substitute(dna,i,base): # perform nucleotide substitution at index i
    dna = dna[:i] + base + dna[i+1:]
    return dna

def diff(sub,ref): # evaluate first instance of nucleotide difference
    if len(sub) == len(ref): 
        for i in range(len(sub)):
            j = ref[i]
            k = sub[i]
            if j != k:
                return i
                break
            else:
                continue
        return -1
    elif len(sub) != len(ref): # use smaller length to account for length differences
        default = min(len(ref),len(sub))
        for i in range(default):
            j = ref[i]
            k = sub[i]
            if j != k:
                return i
                break
            else:
                continue
        return default

def repair(sub,ref): # repair sub sequence using optimal adjustment function
    index_diff = diff(sub,ref)
  
    if index_diff == -1:
        return sub
    
    new_sub = substitute(sub,index_diff,ref[index_diff])
    substituted = diff(new_sub,ref)
        
    deletion = delete(sub,index_diff)
    deletion_diff = diff(deletion,ref)
    
    insertion = insert(sub,index_diff,ref[index_diff])
    insertion_diff = diff(insertion,ref)
        
    repairs = {index_diff:sub, substituted:new_sub, 
               deletion_diff:deletion, insertion_diff:insertion}
                
    if -1 in repairs.keys():
        bestSub = repairs[-1]
        return bestSub
    else:
        bestSub = max(index_diff, substituted, deletion_diff, insertion_diff)
        return repairs[bestSub]

def count(sub,ref): # count how many times repair function is called, until ref sequence is achieved
    if sub == ref:
        return 0
    else:
        counter = 0
        while sub != ref:
            sub = repair(sub,ref)
            counter += 1
            if sub == ref:
                break
        return counter
    
def main():
    ref = readFile('ref.txt')
    dna1 = readFile('dna1.txt')
    dna2 = readFile('dna2.txt')   
    dna3 = readFile('dna3.txt') 
    dna_samples = [dna1,dna2,dna3]
    
    for dna in dna_samples: # iterate through sample files
        n = count(dna,ref)
        if synonymous(dna,ref):
            synon = 'synonymous'
        else:
            synon = 'not synonymous'
        
        print(f'Subject {dna_samples.index(dna) + 1} DNA has {n} mutations and is {synon}')
    
if __name__ == "__main__":
    main()
