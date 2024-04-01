import random
# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

def Motifs(Profile, Dna):
    most_prob=[]
    le=len(Profile['A'])
    for i in range(len(Dna)):
        curr= ProfileMostProbableKmer(Dna[i],le,Profile)
        most_prob.append(curr)
    return most_prob

def RandomMotifs(Dna, k, t):
    prob=[]
    for i in range(t):
        num=random.randint(0, len(Dna[i])-k)
        curr=Dna[i][num:num+k]
        prob.append(curr)
    return prob


def Count(Motifs): #generate a count matrix
    count = {} # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = [0] * len(Motifs[0])#matrix of 0
    for i in range(len(Motifs)):# each lane forms a str
        for j in range(len(Motifs[i])): #nt in each str
            symbol = Motifs[i][j]#each nt
            #if symb in count:#similar to dict.get(symbol)
            count[symbol][j] += 1#  add to the index
    return count

def Profile(Motifs):#give a relative frequency distribution of each nt relative to total str of the motifs
    Lis=Count(Motifs)
    profile = {}
    for symbol in"ACGT":
        profile[symbol]=[]
        for nt in Lis[symbol]:#nt is the frequency of each nucleotide, not int of range method
            profile[symbol].append(nt/len(Motifs))
    return profile

def Consensus(Motifs): #the most popular nucleotides in each column of the motif matrix (ties are broken arbitrarily)
    count = Count(Motifs)
    consensus = ""
    for j in range(len(Motifs[0])):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m: # commits to access both from both variables in the dictionary
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs): #The score/number of differences within each string/ k-mers.
    score=0
    base=Count(Motifs)
    ref=Consensus(Motifs)
    for symbol in Count(Motifs):
        for j in range(len(Motifs[0])):
            if ref[j]!=symbol:
                m=base[symbol][j]
                score+=m
    return score

def Pr(Text, Profile):#probabilities calculated from profile
    p=1
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(Text, k, Profile):
    most={}
    for i in range(len(Text)-k+1):
        current=Text[i:i+k]
        most[current]=Pr(current,Profile)
        key_list = list(most.keys())
        val_list = list(most.values())
        position = val_list.index(max(val_list))
    return (key_list[position])

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):#setting first k-mer in each string as bestmotifs for initiation
        BestMotifs.append(Dna[i][0:k]) 
    for i in range((len(Dna[0])-k+1): #comparison repeats until finds best motifs of the cycle
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):#motifs[0] used for profile for 2nd-5th str
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))#return motif with highest probability
        if Score(Motifs) < Score(BestMotifs):#return motif with least difference
            BestMotifs = Motifs
    return BestMotifs

import math

def log2(x):#entropy calculation
    entropy= x*math.log(x ,2) #logarithm base 2 of x
    return entropy (simplified entropy calculation)
    
def CountWithPseudocounts(Motifs): #normalize all possiblity with +1
    count = {} # initializing the count dictionary
    t = len(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)#matrix of 1 instead of 0
    for i in range(len(Motifs)):# each lane forms a str
        for j in range(len(Motifs[i])): #nt in each str
            symbol = Motifs[i][j]#each nt
            #if symb in count:#similar to dict.get(symbol)
            count[symbol][j] += 1#  add to the index
    return count

def ProfileWithPseudocounts(Motifs):#normalized profile
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} 
    Lis=CountWithPseudocounts(Motifs)
    # insert your code here
    for symbol in"ACGT":
        profile[symbol]=[]
        for nt in Lis[symbol]:
            #nt is the frequency of each nucleotide, not int of range method
            profile[symbol].append(nt/(t+4)) # 1 appended, count[symbol] is 4 larger than regular motifs
    return profile# output variable

def Normalize(Probabilities):
    pre_total=0
    for i in Probabilities:
        pre_total+=Probabilities[i]
    fac=1.00/pre_total
    for i in Probabilities:
        Probabilities[i]=Probabilities[i]*fac
    return Probabilities

def WeightedDie(Probabilities):
    kmer = '' # output variable
    p=random.uniform(0, 1)
    odds=0
    for i in Probabilities:
        odds+=Probabilities[i]#continue to add probabilities until greater than p
        if p<odds:
            return i

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    probabilities = {}
    for i in range(len(Text)-k+1):#range over all possible k-mers in Text, computing the probability of each one and placing this probability into a dictionary
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)
    # normalize these probabilities using the Normalize(probabilities) subroutine, and then return the result of rolling a weighted die over this dictionary to produce a k-mer.

def GibbsSampler(Dna, k, t, N):#run N times
    BestMotifs = [] # output variable
    Motifs=RandomMotifs(Dna, k, t)#randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    BestMotifs = Motifs #BestMotifs ← Motifs
    for j in (1, N): # run N times
        i=random.randint(0, t-1)  #select random string in DNA
        temp_motifs=[]
        for a in range(t):
            if a!=i: #Profile ← profile matrix formed from all strings in Motifs except for Motifi
                temp_motifs.append(BestMotifs[a])#set of random motifs from each string
        Profile=ProfileWithPseudocounts(temp_motifs) 
        Motifs[i]= ProfileGeneratedString(Dna[i], Profile, k)#most probable kmers in ith sring
            BestMotifs = Motifs #best motifs =original Motifs, Motifs included most probable kmers from last step.
    return BestMotifs
