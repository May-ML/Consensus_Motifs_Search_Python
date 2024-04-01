def FasterSymbolArray(Genome, symbol):
    array = {}
    n=len(Genome)
    ExtendedGenome=Genome+Genome[0:n//2]
    array[0] =PatternCount(Genome[0:n//2],symbol)
    for i in range(1,n):
        array[i]=array[i-1]
        if ExtendedGenome[i-1] ==symbol:
            array[i]=array[i]-1
        if ExtendedGenome[i+(n//2)-1]==symbol:
            array[i]=array[i]+1
    return array

# Input:  Strings Text and Pattern
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Text,Pattern):
    count = 0 # output variable
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def SkewArray(Genome):
    skew={}
    array={'A':0,'T':0,'G':1,'C':-1}
    for i in range(len(Genome)):
        skew[0]=0
        if Genome[i] in array:
            skew[i+1]=skew[i]+array.get(Genome[i])
    return skew.values()

##def SkewArray(Genome):
##    skew = {0:0}
##    values = {'A':0, 'T':0, 'C':-1, 'G':1}
##    for i in range(0, len(Genome)):
##        skew[i+1] = skew[i] + values[Genome[i]]
##    return skew.values()

####def SkewArray(Genome):
####    result = [0]
####    for g in Genome:
####        result.append(result[-1] + int(g == 'G') - int(g == 'C'))
####    return result

# find the minimum value of all values in the skew array
# range over the length of the skew array and add all positions achieving the min to positions
def MinimumSkew(Genome):
    positions = [] # output variable
    Range=list(SkewArray(Genome))
    for i in range(len(Genome)):
        if Range[i]==min(Range):
            positions.append(i)
    return positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    count=0
    for i in range(len(p)):
        if p[i]!=q[i]:
            count+=1
    return count

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)):
        current = Text[i:i+len(Pattern)]
        if HammingDistance(current, Pattern)<=d and len(current)==len(Pattern):
            positions.append(i)
    return positions

