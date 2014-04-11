     __      __                  .___          ___.                 
    /  \    /  \__ __  ____    __| _/__________\_ |__ _____ _______ 
    \   \/\/   /  |  \/    \  / __ |/ __ \_  __ \ __ \\__  \\_  __ \
     \        /|  |  /   |  \/ /_/ \  ___/|  | \/ \_\ \/ __ \|  | \/
      \__/\  / |____/|___|  /\____ |\___  >__|  |___  (____  /__|   
           \/             \/      \/    \/          \/     \/

About
=====

Wunderbar is a program which is able to identify mislabeled samples in genotype data 
when a few extra independent genotypes are available ( “genetic barcoding“ ). 

Wunderbar calculates probabilities of the likelihood that genotype mismatches have occurred by chance. 
It is capable of reliable and sensitive detection of sample mismatches and swaps even in the presence
of numerous genotyping errors and in the presence of linkage disequibrilium between the individual 
genotypes. It only requires a few SNPs to work well.


Usage
=====

Wunderbar requires files that in plink format, see http://pngu.mgh.harvard.edu/~purcell/plink/ .
Moreover, the files must be in the non-binary, with alleles recoded numerically and tab-separated.
To convert a binary plink file to the required format, use a command similar too:

    plink --bfile my-binary-plink --recode12 --tab --out my-wunder-ready-plink

Wunderbar expects that you have two plink file-sets: A file-set with barcoding genotypes and 
a file-set with the genotypes that you with to do barcoding against. Then to run barcoding
you would use a command like:

    ./wunderbar --barcode-file my-barcodes --file my-test-genotypes --sampling-iterations 10000 --out results

This will create the files ```results.mismatch``` (candidate list of mismatching samples) and ```results.swaps```
(candidate list of mixed up samples).

The full list of options are:

    Usage: wunderbar [options]
        -h, --help                       Display this screen
        -f, --file plinkstem             The genotypes to be tested
        -b, --barcode-file plinkstem     The barcode genotypes
        -o, --out file                   Prefix of the output files
            --ld-adjust iterations       Adjust for population LD structure between barcoding SNPs. Default number of iterations is 0 (Assume no LD). The more iterations the more precise the adjustment.
            --sampling-iterations iterations
                                     Number of iterations to use to calculate p-value from mismatch probability distribution
        --peturb 0..1                A very small pertubation constant to ensure that no sample gets zero mismatch probability (default 0)
        --dot-file file              Create .dot file for GraphViz to visualize swaps


License
=======

    Copyright (C) 2014 Christian Theil Have

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


