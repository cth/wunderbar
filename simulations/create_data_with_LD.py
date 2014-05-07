#!/usr/bin/env python
import argparse
from random import randrange,seed,random
from math import floor

parser=argparse.ArgumentParser(description='Create simulation data for wunderbar')
parser.add_argument('path', type=str, help='String with the path')
parser.add_argument('--N-ind', type=int, nargs=1, default=[1000], help='Number of individuals. Default 1000')
parser.add_argument('--N-chip-SNPs', type=int, nargs=1, default=[12], help='Number of chip SNPs. Default 12')
parser.add_argument('--N-bar-SNPs', type=int, nargs=1, default=[12], help='Number of bar SNPs. Can not be bigger than N-chip-SNP. Default 12')
parser.add_argument('--F-chip', type=str, nargs=1, default=['chip'], help='Name of plink stem file for the chip. Default=chip')
parser.add_argument('--F-bar', type=str, nargs=1, default=['bar'], help='Name of plink stem file for the barcorde. Default=bar')
parser.add_argument('--Mismatch-rate', type=float, nargs=1, default=[0.01], help='Mismatch rate in barcode. Default 0.01')
parser.add_argument('--Missing', type=bool, nargs=1, default=[False], help='Simulate missing data <True|False>. Default: False')

args = parser.parse_args()
print(args)

path = args.path

with_missing = args.Missing[0]
mismatch_rate = args.Mismatch_rate[0]
number_of_snps = args.N_chip_SNPs[0]
number_of_individuals = args.N_ind[0]
number_of_barsnps = args.N_bar_SNPs[0]

chip_ped = open(path + '/' + args.F_chip[0] + '.ped','w')
chip_map = open(path + '/' + args.F_chip[0] + '.map','w')
bar_ped = open(path + '/' + args.F_bar[0] + '.ped','w')
bar_map = open(path + '/' + args.F_bar[0] + '.map','w')


def recoder(snp):
    #Simply recode 0,1,2,3 into a 1 1, 1 2, 2 2 and 0 0
    recodedsnp = ''
    if(snp==1):
        recodedsnp = '1 1'
    elif(snp==2):
        recodedsnp = '1 2'
    elif(snp==3):
        recodedsnp = '2 2'
    elif(snp==0):
        recodedsnp = '0 0'
                      
    return(recodedsnp)


seed()
chip = []
barcode = []
for i in range(0,number_of_individuals):
    barindividual = []
    individual = []
    snps = []

    #First SNP should never be missing
    snp = randrange(0,3,1)
    snps.append(snp%3+1)



    for j in range(1,number_of_snps):
        snps.append(snps[0]) #Filling out the rest of the snps as the first snp

    #Mismatch to simulate not-perfect LD
    for b in range(0,number_of_snps):
        mismatchsnp = snps[b]
        if(random()<0.1): #10% mismatch rate = ~90%LD
            while mismatchsnp == snps[b]:
                if(with_missing):
                    mismatchsnp = randrange(0,4,1)%4
                    if(mismatchsnp==0): #For 0.25*0.25=0.125 chance for missing data
                        seed()
                        mismatchsnp = randrange(0,4,1)%4
                else:
                    mismatchsnp = randrange(0,3,1)%3+1

            snps[b]=mismatchsnp
    

    IID = 'ind'+str(i)
    FID = 'ind'+str(i)
    MAT = '0'
    PAT = '0'
    PHE = '-9'
    SEX = str(randrange(1,3,1))
    individual.append(IID)
    individual.append(FID)
    individual.append(MAT)
    individual.append(PAT)
    individual.append(PHE)
    individual.append(SEX)
    individual.extend(snps)
    chip.append(individual)

    for b in range(0,number_of_barsnps):
        mismatchsnp = snps[b]
        if(random()<mismatch_rate):
            while mismatchsnp == snps[b]:
                if(with_missing):
                    mismatchsnp = randrange(1,12,1)%4
                else:
                    mismatchsnp = randrange(1,12,1)%3+1

            snps[b]=mismatchsnp

    barindividual.append(IID)
    barindividual.append(FID)
    barindividual.append(MAT)
    barindividual.append(PAT)
    barindividual.append(PHE)
    barindividual.append(SEX)
    barindividual.extend(snps)
    barcode.append(barindividual)
    

for ind in chip:
    formatted_ind = ind[0:6]
    formatted_ind.extend(map(recoder,ind[6:18]))
    for find in formatted_ind:
        chip_ped.write(find + '\t')
    chip_ped.write('\n')

for i in range(1,len(snps)+1):
    chip_map.write(str(i)+'\trs'+str(i)+'\t0\t'+str(i)+'\n')

chip_map.close()
chip_ped.close()

for ind in barcode:
    formatted_ind = ind[0:6]
    formatted_ind.extend(map(recoder,ind[6:18]))
    for find in formatted_ind:
        bar_ped.write(find + '\t')
    bar_ped.write('\n')

for i in range(1,number_of_barsnps+1):
    bar_map.write(str(i)+'\trs'+str(i)+'\t0\t'+str(i)+'\n')

bar_map.close()
bar_ped.close()
