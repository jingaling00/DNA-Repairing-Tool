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
        
def pointMutation(file):
    return None

def transcribe(dna):
    mrna = []
    for base in dna:
        if base == 'T':
            base = 'U'
        mrna.append(base)
    mrna = ''.join(mrna)
    return mrna

# perform transcription Thymine to Uracil

def translate(mrna):
    n = len(mrna)
    parsed_mrna = [mrna[i:i+3] for i in range(0,n-1,3)]
    amino_acids = []
    for codon in parsed_mrna:
        if len(codon) == 3:
            codon = codons_dict[codon]
            amino_acids.append(codon)
    aa_sequence = ' '.join(amino_acids)
    return aa_sequence
# perform mRNA translation

def synonymous(sub,ref):
    transcribe_sub = transcribe(sub)
    protein_sub = translate(transcribe_sub)
    
    transcribe_ref = transcribe(ref)
    protein_ref = translate(transcribe_ref)
    
    if protein_sub == protein_ref:
        return True
    else:
        return False
# encode same AA

def delete(dna,i):
    dna = dna[:i] + dna[i+1:]
    return dna # delete dna at index i

def insert(dna,i,base):
    dna = dna[:i] + base + dna[i:]
    return dna

def substitute(dna,i,base):
    dna = dna[:i] + base + dna[i+1:]
    return dna

def diff(sub,ref):
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
    elif len(sub) != len(ref):
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
            # if all bases are same up to extraneous
    #index of first difference 

def repair(sub,ref):
    
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

def count(sub,ref):
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
    
 # total number of point mutations
    # count number of times repair runs
    

if __name__ == '__main__':
    ref = readFile('ref.txt')
    dna1 = readFile('dna1.txt')
    dna2 = readFile('dna2.txt')   
    dna3 = readFile('dna3.txt') 
    dna_samples = [dna1,dna2,dna3]
    
    for dna in dna_samples:
        n = count(dna,ref)
        
        if synonymous(dna,ref):
            synon = 'synonymous'
        else:
            synon = 'not synonymous'
        
        print(f'Subject {dna_samples.index(dna) + 1} DNA has {n} and is {synon}')
