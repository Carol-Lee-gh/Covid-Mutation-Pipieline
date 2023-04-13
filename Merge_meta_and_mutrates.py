import datetime
import argparse
import ujson
import os
from tqdm import tqdm
import pandas as pd
import csv
import bz2
import re

minDays =14
minSamples = 10

parser = argparse.ArgumentParser(description = 'Merge metadata with mutation and hap output, ensure no duplicates present')
parser.add_argument('--metaFile','-m', help = 'Input metadata file (.json or .csv)', required = True)
parser.add_argument('--prefix', '-p', help = 'Hapout prefix from previous script, default is data/countryratios.csv and mutationratios.csv')
parser.add_argument('--out', '-o', help = 'Output prefix, default is data/countryratios.csv and mutationratios.csv')
parser.add_argument('--threshold', '-t', nargs = '?', default = 1.2, help = 'R2 Threshold to filter countries by, default is 1.2. Specify value to change')
args = parser.parse_args()

metaFile = args.metaFile
hapout = args.prefix
R2_threshold = float(args.threshold)
out = args.out
prefix = args.prefix

path = os.getcwd()+ '/output/'

#Get column names
with open(metaFile, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    list_of_column_names = []
    for row in csvreader:
    # loop to iterate through the rows of csv 
        # adding the first row
        list_of_column_names.append(row)
 
        # breaking the loop after the
        # first iteration itself
        break
# printing the result
print('List of column names : ',
      list_of_column_names[0])
      
### Collecting geographic and time location
#Dates are recorded in the format YYYY-MM-DD. These are converted into integers where Day=1 was the day the first sequence deposited into GISAID was collected.

start = datetime.datetime.strptime('2019-12-30', '%Y-%m-%d')
#end = '2021-07-13' # uncomment for manual input 
end = str(datetime.date.today()) #generates todays date
end = datetime.datetime.strptime(end, '%Y-%m-%d')
date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end-start).days)]
time=[]
times2days = dict()
for date in date_generated:
    times = date.strftime("%Y-%m-%d")
    time.append(times)

i = 1
for times in time:
    times2days[times] = i
    i += 1

dict_from_csv = {}

with open(metaFile, mode='r',encoding='utf-8') as inp:
    reader = csv.reader(inp)
    dict_from_csv = {rows[0]:rows[0:] for rows in reader}

print('Reading hapdict file: ' + path + prefix + '_mutlist.txt)')
with open(path + prefix + '_mutlist.txt', 'r') as fp:
    mutList = ujson.load(fp)

print('Reading hapdict file: ' + path + prefix + '_hapdict.txt)')
with open(path + prefix + '_hapdict.txt','r') as f:
    hapDict = ujson.load(f)

mutData = dict()

for key in tqdm(hapDict.keys(), desc = 'Merging metadata and hapdict'):
    if key in dict_from_csv:
        mutData[key] = list(dict_from_csv[key]) + list(hapDict[key])

mutations = list(mutList.keys())
columnNames = list_of_column_names[0] + mutations
data = pd.DataFrame.from_dict(mutData, orient='index', columns=columnNames)
data['Days']= data['Date'].apply(lambda x: [v for k, v in times2days.items() if k in x][0])

data = data[data['Date'].isin(time)]
data['Date'] = pd.to_datetime(data['Date'])
data['Country'] = data['Location'].str.split('/').str[0]
first_column = data.pop('Country')
data.insert(1, 'Country', first_column)

# remove any duplicate columns
data = data.loc[:,~data.columns.duplicated()]

'''
# Check whether the specified path exists or not
isExist = os.path.exists(path)

if not isExist:
  
  # Create a new directory because it does not exist 
    os.makedirs(path)
print('Created data directory: ' + path)
'''
#output raw counts of mutations, refs and alt alleles with sampleID, Country, Date, Days, mutation (categorised),etc to csv for ratio matrix calculation/further analyses
data.set_index('Date', inplace=True)
data.to_csv(path + out + 'ratioData.csv.gz',compression='gzip')

#calculate ratios per country and per SNP for Snps that have been present for >14 days
countryMutations = {}
countries = data.Country.unique()
print('\nStarting R2 ratio analysis...')
for country in tqdm(countries):
    countryMutations[country] = {}
    countryData = data[data.Country.eq(country)]
    for m in mutations:
        subset = countryData[['Days', m]]
        subset.columns = ['Days', 'Mutation']
        SNP_present = subset[subset.Mutation.eq(1)]
        SNP_not_present = subset[subset.Mutation.eq(0)]
        SNPSamples = SNP_present.shape[0]
        RefSamples = SNP_not_present.shape[0]
        useOldMethod = True
        if SNPSamples > minSamples and RefSamples > minSamples:
            SNPDays = SNP_present['Days'].astype('int').max() - SNP_present['Days'].astype('int').min() + 1 #days between first and last SNP sample
            RefDays = SNP_not_present['Days'].astype('int').max() - SNP_not_present['Days'].astype('int').min() + 1 # days between first and last Ref sample
            if SNPDays >= minDays and RefDays >= minDays:
                useOldMethod = False
                SNPRate = SNPSamples / SNPDays
                RefRate = RefSamples / RefDays
                Ratio = SNPRate / RefRate
                countryMutations[country][m] = Ratio
        if useOldMethod:
            countryMutations[country][m] = 0
print('Calculated R2 ratios')

# Print all countries with any ratio > threshold (default  = 1.2)
#Simply displays all countries of interest with a high ratio for at least one of the selected SNPs. The threshold can be set in command.
MutationRatios = pd.DataFrame.from_dict(countryMutations, orient='index')

output = MutationRatios[(MutationRatios >= R2_threshold).any(1)] 

print('writing to output...')
output.to_csv(path + out + 'countryRatios_' + str(R2_threshold) + '.csv')
print('Done')      
