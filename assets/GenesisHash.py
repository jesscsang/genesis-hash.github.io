from time import perf_counter
from math import floor


def remappedASCII(ASCII):
    """
    95 Printable characters 32 - 126
    https://en.wikipedia.org/wiki/ASCII
    https://www.cs.cmu.edu/~pattis/15-1XX/common/handouts/ascii.html

    Remaps ASCII value to some unique 12-14 digit value (Tested).
    """
    LEADING_ZERO = "0"

    strASCII = str(ASCII)

    # Use string methods to separate ASCII digits
    if ASCII < 10:
        # Single digit ASCII values
        L_0 = L_1 = LEADING_ZERO
        L_2 = strASCII
    elif ASCII < 100:
        # Two digit ASCII values
        L_0 = LEADING_ZERO
        L_1, L_2 = strASCII
    else:
        # Three digit ASCII values
        L_0, L_1, L_2 = strASCII

    # Convert separate string digits back to int
    L_0 = int(L_0)
    L_1 = int(L_1)
    L_2 = int(L_2)

    """
    Before modulus:
        min = 29124297227663 (14 digits)
        max = 325620017749655509090396059873 (30 digits)
    """
    return floor(floor(((L_0+3)**10 + (L_1+3)**10)/3) /
                 (L_2+3))**3 % (10**14)


def strToRemappedASCII(string):
    resultASCII = ""
    for char in string:
        resultASCII += str(remappedASCII(ord(char)))
    return resultASCII


def appendEvenOrOddIndices(str, type="even", length=-1):
    result = ""
    i = 0
    counter = 0
    lengthReached = False
    while (i < len(str) and not lengthReached):
        if (i % 2 == 0 and type == "even"):
            result += str[i]
            counter += 1
        if (i % 2 != 0 and type != "even"):
            result += str[i]
            counter += 1
        i += 1
        if (length != -1 and counter == length):
            lengthReached = True
    return int(result)


def intronSplicing(decimal, threshold=14):
    """
    5'-GU...AG-3'
    5'-1001...0010-3'
    5'-oeeo...eeoe-3'
    """
    strDecimal = str(decimal)
    dnaLen = len(strDecimal)
    if dnaLen >= 8:
        introns = []
        startRemove = -1
        endRemove = -1
        startFound = False
        endFound = False
        result = ""
        for i in range(dnaLen-7):
            # Start of Intron
            if (not startFound and
                int(strDecimal[i]) % 2 == 1 and
                int(strDecimal[i+1]) % 2 == 0 and
                int(strDecimal[i+2]) % 2 == 0 and
                    int(strDecimal[i+3]) % 2 == 1):
                startRemove = i
                startFound = True
            # End of Intron
            if (startFound and
                int(strDecimal[i+4]) % 2 == 0 and
                int(strDecimal[i+5]) % 2 == 0 and
                int(strDecimal[i+6]) % 2 == 1 and
                    int(strDecimal[i+7]) % 2 == 0):
                endRemove = i+7
                endFound = True
            # Remove Intron
            if endFound:
                startFound = False
                endFound = False
                introns.append(strDecimal[startRemove:endRemove+1])
        for intron in introns:
            intronLen = len(intron)
            if (dnaLen-intronLen >= threshold):
                dnaLen -= intronLen
                strDecimal = strDecimal.replace(intron, "")

    return int(strDecimal)


def decimalToBin(decimal):
    strDecimal = str(decimal)
    binary = ""
    for i in range(len(strDecimal)):
        binary += str(int(strDecimal[i]) % 2)
    return binary


def binToRNA(bin):
    DNA_Result = ""
    bin2DnaMap = {"00": "A", "01": "U", "10": "G", "11": "C"}

    for i in range(0, len(bin)-1, 2):
        binaryPair = bin[i]+bin[i+1]

        for key in bin2DnaMap:
            if (key == binaryPair):
                DNA_Result += bin2DnaMap[key]

    return DNA_Result


def RNA2aminoAcids(RNA):
    aaResult = ""
    RNA_Codons = {
        # 'M' - START, '_' - STOP
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UGU": "C", "UGC": "C",
        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "UUU": "F", "UUC": "F",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAU": "H", "CAC": "H",
        "AUA": "I", "AUU": "I", "AUC": "I",
        "AAA": "K", "AAG": "K",
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUG": "M",
        "AAU": "N", "AAC": "N",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UGG": "W",
        "UAU": "Y", "UAC": "Y",
        "UAA": "_", "UAG": "_", "UGA": "_"
    }

    for i in range(0, len(RNA)-2, 3):
        rnaCodon = RNA[i]+RNA[i+1]+RNA[i+2]
        for key in RNA_Codons:
            if (key == rnaCodon):
                aaResult += RNA_Codons[key]
    return aaResult


def desiredHashLength(_hash, hashLength=256):
    while _hash < 10**(hashLength*2*3*2):
        """
        *2: takes two bits/digits to map to a DNA base 
        *3: takes three DNA bases to map to a RNA codon
        *2: double the length as appendEvenOrOddIndices() will reduce the length by at most half
        """
        _hash = _hash**3

    _hash = appendEvenOrOddIndices(str(_hash), "odd", hashLength*2*3)
    return _hash


def genesisHash(password):
    remappedASCII_PW = strToRemappedASCII(password)
    _hash = appendEvenOrOddIndices(
        remappedASCII_PW, "even") + appendEvenOrOddIndices(remappedASCII_PW, "odd")

    for i in range(3):
        _hash = desiredHashLength(_hash)
        _hash = intronSplicing(_hash)
    _hash = desiredHashLength(_hash)

    binary = decimalToBin(_hash)

    RNA = binToRNA(binary)

    AA = RNA2aminoAcids(RNA)

    return AA


password = "lemon"
# start = perf_counter()
# _hash = genesisHash(password)
# end = perf_counter() - start
# print(f'One 8 character password: {end} sec')
# print(f'All 8 character passwords: {(end*95**3)/(60*60*24)} days')
print(f'hash: {genesisHash(password)}\n')


# counter = 0
# for c1 in (chr(i) for i in range(32, 127)):
#     for c2 in (chr(i) for i in range(32, 127)):
#         for c3 in (chr(i) for i in range(32, 127)):
#             if counter % 100000 == 0:
#                 fname = 'bruteForce3CharPW_'+str(int(counter/100000))+'.txt'
#             counter += 1
#             with open(fname, 'a') as f:
#                 password = c1+c2+c3
#                 f.write(f'{password},{genesisHash(password)}\n')
#                 print(f'{password},{genesisHash(password)}\n')
