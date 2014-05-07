#!/usr/bin/env python
import argparse
import os
import random
random.seed()

parser=argparse.ArgumentParser(description='Switch a number of individuals in a ped file and copies the map file with the extension .switched. The lists is recorded in list.of.switch')
parser.add_argument('path', type=str, help='String with the path')
parser.add_argument('-N', type=int, nargs=1, default=[1], help='Number of switches. Default 1')
parser.add_argument('-F', type=str, nargs=1, default=['chip'], help='name of plink stem file. Default=chip')
args = parser.parse_args()

path = args.path
numberofswitch = args.N[0]
plink_stem = args.F[0]

nlist=[]
mlist=[]
pedfile = open(path + '/' + plink_stem+'.ped','r').readlines()
for i in range(0,numberofswitch):
    n=0
    m=0
    while n==m or (n in nlist or m in mlist) or (n in mlist or m in nlist):
       n=random.randrange(0,len(pedfile),1)
       m=random.randrange(0,len(pedfile),1)
    
    nlist.append(n)
    mlist.append(m)


switchedpedfile = open(path + '/' + plink_stem + '.switched.ped','w')
listofswitchs = open(path + '/' + plink_stem + '.list.of.switch','w')

switchlist = pedfile
switchlist[len(switchlist)-1] = switchlist[len(switchlist)-1] + '\n'

for l,line in enumerate(switchlist):
    switchlist[l] = line.split()

for k in range(0,numberofswitch):
    switchlist[nlist[k]][6:], switchlist[mlist[k]][6:] = switchlist[mlist[k]][6:],switchlist[nlist[k]][6:]
    listofswitchs.write(pedfile[nlist[k]][0] + " " + pedfile[nlist[k]][1] +";"+pedfile[mlist[k]][0]+ " " + pedfile[mlist[k]][1] + "\n")

for j in switchlist:
    for i in range(0,6):
        switchedpedfile.write(j[i] + '\t')
    for i in range(6,len(j),2):        
        switchedpedfile.write(j[i] + ' ' + j[i+1] + '\t')

    switchedpedfile.write('\n')

os.system('cp '+path + '/' + plink_stem + '.map ' + path + '/' + plink_stem + '.switched.map')
