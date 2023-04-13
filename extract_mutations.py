from Bio import SeqIO
from tqdm import tqdm
import json
import re
import argparse
from collections import defaultdict
import os

# Define inputs
parser = argparse.ArgumentParser(description="Script for extracting variants of interest")
parser.add_argument('--mutFile', '-m', help='mutation list (csv)', required=True)
parser.add_argument('--vcfFile', '-vcf', help = 'input vcf file', required = True)
parser.add_argument('--prefix', '-out', help = 'output prefix, default is to data/ _mut and haplist.txt')
args = parser.parse_args()

mutFile = args.mutFile
vcfFile = args.vcfFile
hapout = args.prefix

mutList = defaultdict(list)

path = os.getcwd()+'/output/'
# Check whether the specified path exists or not
isExist = os.path.exists(path)

if not isExist:
  
  # Create a new directory because it does not exist
    os.makedirs(path)
print('Created data directory: ' + path)

# Add list of mutations in .csv into mutList()
with open(mutFile, 'r') as file:
    print("Adding mutation list...")
    data = json.load(file)

    for keys in data:
        parts = data[keys].rstrip().split(',')

        m = ()
        for i in range(0, int(len(parts) / 3)):
            m = m + (parts[0 + (i * 3)], parts[2 + (i * 3)])
        mutList[keys].append(m)
        #print(mutList[:5])
    print("Finished adding mutation list")

with open(path + hapout + "_mutlist.txt", "w") as fp:
        json.dump(mutList,fp)
print("Saved mutlist to " + path + hapout + "_mutlist.txt")
        
posColumn = 1  # column of the base position
refColumn = 3  # column of the ref base identity
altColumn = 4  # column of the observed alt bases

### This code makes a default 1:1 mapping of ref_alignment to actual for cases where it's not important
ref2actual = dict()
actual2ref = dict()
for i in range(1,50000):
        ref2actual[i] = i
        actual2ref[i] = i

# Parse through VCF file and extract lines of interest
vcfPosList = list()
for key in mutList:
    for m in mutList[key]:
        for i in range(0,int(len(m)/2)):
            if not ref2actual[int(m[0+(i*2)])] in vcfPosList:
                vcfPosList.append(ref2actual[int(m[0+(i*2)])])
lineDict = dict()

snpLineStart = "1" # pipeline has CHROM=1, so every SNP line will start with 1

with open(vcfFile, 'r') as inFile:
    print("Extracting variants of interest...")
    for line in tqdm(inFile):
        if line.startswith("#CHROM"):
            sampleID_line = line # important to store the line with all sample information so we no what order the samples are in
            lineDict['SAMPLE_ID_LINE'] = sampleID_line

        elif line.startswith(snpLineStart):
            columns = line.split()

            if int(columns[posColumn]) in vcfPosList: # check if this is a position we're interested in, if it is store it
                lineDict[actual2ref[int(columns[posColumn])]] = line
                vcfPosList.remove(int(columns[posColumn])) # once you've found it, remove from the list
    print("Extracted variants of interest")

def extractHap(altString, altHap):
# this function determines how the alt variant is annotated in the vcf line
    if altString == 0:
        return(1)
    else:
        alts = altString.split(',')
        for i in range(0,len(alts)):
            if alts[i] == altHap:
                return(int(i+1))

def categorizeHap(haps, ref, alt):
# this function categorizes the haplotypes of each sample as ref (0), alt (1) or other (-1)
    cat = list()
    for hap in haps:
        if hap == ref:
            cat.append(0)
        elif hap == alt:
            cat.append(1)
        else:
            cat.append(-1)
    return(cat)

def categorizeComboHaps(haps, alts):
    cat = list()
    for i in range(0,len(haps)):
        h = list(haps[i])
        if h == alts:
            cat.append(1)
        else:
            cat.append(0)
    return(cat)
    
sampleList = list() # list of all accession ids in order the appear 
for sample in sampleID_line.rstrip().split()[9:]:
    sampleList.append(sample.split("|")[0])

hapList = list()
hapDict = dict()

print("Categorizing variants...")
for key in mutList:
    for m in tqdm(mutList[key]):
        if len(m) > 2: # if it's a combination mutation treat it differently. Here 1 = all mutatns present, 0 = not all mutants present

            snpColumns = list()
            altHaps = list()

            for i in range(0,int(len(m)/2)):
                mutPos = m[0+(i*2)]
                if int(mutPos) in lineDict:
                    snpColumns.append(lineDict[int(mutPos)].split()[:9] +  [ int(x) for x in lineDict[int(mutPos)].split()[9:] ])
                else:
                    snpColumns.append([0] * len(sampleID_line))

                altHaps.append(extractHap(snpColumns[i][altColumn], m[1+(i*2)]))

                snpHaps = categorizeComboHaps(list(zip(*snpColumns))[9:], altHaps)
        else: # if it's just a single SNP: 1 = mutant present, 0 = reference present, -1 = something else present (other haploytpe or ambiguous base)
            mutPos,mutAlt = m
            if int(mutPos) in lineDict: # check if variant was actually found at location in vcf file
                snpColumns = lineDict[int(mutPos)].split()[:9] + [ int(x) for x in lineDict[int(mutPos)].split()[9:] ]
            else: # if not, set all haplotypes to reference (0)
                snpColumns = [0] * len(sampleID_line)

            refHap = 0 # the reference is always haplotype 0
            altHap = extractHap(snpColumns[altColumn], mutAlt)

            snpHaps = categorizeHap(snpColumns[9:], refHap, altHap)

        for i in range(0, len(sampleList)):
            if sampleList[i] in hapDict:
                hapDict[sampleList[i]].append(snpHaps[i])
            else:
                hapDict[sampleList[i]] = [snpHaps[i]]
with open(path + hapout + "_hapdict.txt","w") as f:
    json.dump(hapDict, f)
print("Saved categorised variants to " + path + hapout + "_hapdict.txt" + "\nDone")
