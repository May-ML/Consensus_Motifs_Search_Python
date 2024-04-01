def PatternCount(Text, Pattern):##compare pattern by sliding one nucleotide at a time
    count = 0
    for i in range(len(Text)-len(Pattern)+1): #range() exclude the last index
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def FrequencyMap(Text, k): #a dictionary of kmers as key adn number of occurence as values
    freq = {}
    for i in range(len(Text)-k+1):
        kmer = Text[i:i+k]
        if kmer in freq:
            count = freq.get(current)
            freq[current]+=1
        else:
            freq[current] = 1
    return freq

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    for key in freq:#calling keys within the dictionary
        if freq[key]== max(freq.values()): #finding most frequent kmer(s)
            words.append(key)# add each key to words whose corresponding frequency value is equal to m
    return words

def Reverse(Pattern):#complement strings
    for char in Pattern:
        char = Pattern[len(Pattern)-1::-1]
    return char

#def Reverse(Pattern):   # Reverse sequence order function 5'>3' -> 3'>5'
#   return Pattern[::-1]

# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).
def Complement(Pattern):
    # your code here
    base={"A":"T","T":"A","C":"G","G":"C"}
    complement=""
    for i in Pattern:
        complement+=base[i]
    return complement

#def Complement(Pattern): # Complement() function 
#    dict = {'A':'T','G':'C','T':'A','C':'G'}
#    return "".join(dict[i] for i in Pattern)

def ReverseComplement(Pattern):   
    return Complement(Reverse(Pattern))

def PatternMatching(Pattern, Genome): #return position of pattern within the genome
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(str(i))
    return positions

def SymbolArray(Genome, symbol): ## reading the numbers of nucelotides within the frame
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol,ExtendedGenome[i:i+(n//2)])
    return array
