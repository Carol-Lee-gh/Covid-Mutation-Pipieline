import datetime
import argparse
import json
from tqdm import tqdm
import pandas as pd
import pickle

parser = argparse.ArgumentParser(description = 'Merge metadata with mutation and hap output')
parser.add_argument('--metaFile','-met', help = 'Input metadata file (.json)', required = True)
parser.add_argument('--prefix', '-p', help = 'Hapout prefix from previous script')
parser.add_argument('--out', '-o', help = 'Output prefix')

args = parser.parse_args()

metaFile = args.metaFile
hapout = args.prefix
out = args.out
prefix = args.prefix

###Heatmap thresholds
R2_threshold = 1.2
R2_threshold2 = 0.1
minDays = 14
minSamples = 10
####

###Generate a list of valid dates
start = datetime.datetime.strptime("2019-12-01", "%Y-%m-%d")
#end = '2021-07-13' # manual input for testing
end = str(datetime.date.today()) #generates todays date
end = datetime.datetime.strptime(end, "%Y-%m-%d")
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)]
time=[]
times2days=dict()

### Collecting geographic and time location
#Dates are recorded in the format YYYY-MM-DD. These are converted into integers where Day=1 was the day the first sequence deposited into GISAID was collected.

for date in date_generated:
    times=date.strftime("%Y-%m-%d")
    time.append(times)
#print(time)
def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)
i=1
for times in time:
    times2days[times] = i
    i+=1

### Store metadata as a tuple: AccessionID -> (Country, Time Collected)
metaInformation = dict()

with open(metaFile, 'r') as jsonFile:
    print("Collecting metadata...")
    for line in tqdm(jsonFile):
        sample = json.loads(line)
        if sample['covv_host'] == 'Human': # only collect human samples
            if sample['covv_collection_date'] in times2days: # check if the collection date provided is valid
                metaInformation[sample['covv_accession_id']] = (sample['covv_location'],sample['covv_location'].split("/")[1].strip(),sample['covv_collection_date'],sample['covv_lineage'], sample['covv_gender'], sample['covv_patient_age'],sample['covv_patient_status'],times2days[sample['covv_collection_date']])
            else: # if not, record it as 0000-00-00
                metaInformation[sample['covv_accession_id']] = (sample['covv_location'],sample['covv_location'].split("/")[1].strip(),sample['covv_collection_date'],sample['covv_lineage'], sample['covv_gender'], sample['covv_patient_age'],sample['covv_patient_status'].'0000-00-00')
                
print("Reading hapdict file..." + prefix + "_hapdict.txt)")
with open(prefix + "_hapdict.txt","r") as f:
    hapDict = json.load(f)

mutData = dict()

for key in (hapDict.keys()):
    if key in metaInformation:
        mutData[key] = list(metaInformation[key]) + list(hapDict[key])

###Ratio Analysis
#preprocess and import
print("Reading mutation list...(" + prefix + "_mutlist.txt)")
with open(prefix + "_mutlist.txt", "r") as fp:
        mutList = json.load(fp)

mutations = [''.join([str(u) for u in t]) for t in mutList]
columnNames = ["Location","Country","Date", "Lineage", "Gender", "Age", "Status","Days"] + mutations
data = pd.DataFrame.from_dict(mutData, orient="index", columns = columnNames)
data = data[data.Days != "0000-00-00"]

# remove duplicate columns
data = data.loc[:,~data.columns.duplicated()]
#output collated sampleID, Country, Date, Days, mutation (categorised) to csv for ratio matrix calculation
data.to_csv(prefix + "_ratioData.csv")

#calculate ratios per country and per SNP for Snps that have been present for >14 days, can be used for downstream analyses (see Kuiper et al. 2021, https://doi.org/10.1101/2021.08.04.455042)
countryMutations = {}
countries = data.Country.unique()
print("\nStarting R2 ratio analysis...")
for country in tqdm(countries):
    countryMutations[country] = {}
    countryData = data[data.Country.eq(country)]
    for m in mutations:
        subset = countryData[["Days", m]]
        subset.columns = ["Days", "Mutation"]
        SNP_present = subset[subset.Mutation.eq(1)]
        SNP_not_present = subset[subset.Mutation.eq(0)]
        SNPSamples = SNP_present.shape[0]
        RefSamples = SNP_not_present.shape[0]
        useOldMethod = True
        if SNPSamples > minSamples and RefSamples > minSamples:
            SNPDays = SNP_present["Days"].astype('int').max() - SNP_present["Days"].astype('int').min() + 1 #days between first and last SNP sample
            RefDays = SNP_not_present["Days"].astype('int').max() - SNP_not_present["Days"].astype('int').min() + 1 # days between first and last Ref sample
            if SNPDays >= minDays and RefDays >= minDays:
                useOldMethod = False
                SNPRate = SNPSamples / SNPDays
                RefRate = RefSamples / RefDays
                Ratio = SNPRate / RefRate
                countryMutations[country][m] = Ratio
        if useOldMethod:
            countryMutations[country][m] = 0
print("Calculated R2 ratios")

# Print all countries with any ratio > threshold (default  = 1.2)
#Simply displays all countries of interest with a high ratio for at least one of the selected SNPs. The threshold can be set at the top.

MutationRatios = pd.DataFrame.from_dict(countryMutations, orient="index")

#MutationRatios.to_csv(out + 'countryratios_raw.csv')
output = MutationRatios[(MutationRatios >= R2_threshold).any(1)] 
output2 = MutationRatios[(MutationRatios >= R2_threshold2).any(1)] 

print("writing to output...")
output.to_csv(out + "countryRatios_" + str(R2_threshold) + ".csv")
output2.to_csv(out + "countryRatios_" + str(R2_threshold2) + ".csv")

print("Saved ratio data to " + out + "countryRatios_" + str(R2_threshold) + ".csv"
     + "\n" + out + "countryRatios_" + str(R2_threshold2) + ".csv")
print("Done")