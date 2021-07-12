#!/usr/bin/env python

# read the impute output and compute a polygenic risk score based on probabilitic dosage

import sys
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Compute polygenic risk score from impute2 output.')
parser.add_argument('-w', metavar='weight_file', default=None,
                   help='Weights and allele encoding of for each SNP in the polygenic score')
parser.add_argument('-o', metavar='output_file', default=None,
                   help='output file')
parser.add_argument('-d', metavar='delimiter', default="\t",
                   help='delimiter in the weight table')
parser.add_argument('-s', metavar='sample_file', default=None,
                   help='oxford format sample file for the gen files')
parser.add_argument('files', metavar='input_files', help='input files', nargs="*")
args = parser.parse_args()

# cutoffs (currently not used)
cutoff = 0
MAF_cutoff = 0
eff_cutoff = 0

# check that arguments are all present else print help

# read the weights
sys.stderr.write("loading weights: %s\n" % args.w)
weights = pd.read_csv(args.w, index_col="snp_id", sep=args.d)

# read the sample file (if present)
if (args.s):
    sample_ids = pd.read_csv(args.s, sep=" ")
    sys.stderr.write("loading sample ids: %s\n" % args.s)
else:
    sample_ids = None

outfile = sys.stdout
if (args.o):
    outfile = open(args.o, "w")

PRS = None

# file format is id1, id2, position, allele1, allele2, pers1 [P(hom allele1),P(het),P(hom allele2)]
# iterate over all input files
for i in range(len(args.files)):
    sys.stderr.write("processing %s\n" % args.files[i])
    f = open(args.files[i])
    line = f.readline()
    while (line):
        items = line.split(" ")
    
        # look up the weight of the SNP
        snp = items[1]
        if (snp in weights.index):
            weight = weights.loc[snp, 'effect_weight']
        else:
            #sys.stderr.write("No weight for snp: %s\n" % snp)
            line = f.readline()
            continue
    
        # make sure the alleles are correctly encoded
        allele1 = items[3]
        allele2 = items[4]
        effect_allele = weights.loc[snp, 'effect_allele']
        if (allele2 == effect_allele):    
            reverse_coding = False
        elif (allele1 == effect_allele):
            reverse_coding = True
        else:
            sys.stderr.write("%s none of the alleles (%s, %s) matches the effect allele (%s)\n" % (snp, allele1, allele2, effect_allele))
            line = f.readline()
            continue

        nsamples = int((len(items) - 5) / 3)
        if (PRS is None):
            PRS = np.zeros(nsamples)
        
        MAF = 0
        efficacy = 0

        # go through the samples
        for i in range(nsamples):
            # check if any of the three probabilities is above the treshold
            ok = False
            dosage = 0
            for gt in range(3):
                posterior = float(items[5 + gt + 3 * i])
                dosage = dosage + gt * posterior
                if (posterior >= cutoff):
                    efficacy = efficacy + 1
                    ok = True
                    MAF = MAF + dosage
                # if (not ok):
                #     dosage = "NA"
            if (reverse_coding):
                dosage = 2 - dosage
            PRS[i] = PRS[i] + dosage * weight
            efficacy = float(efficacy) / nsamples
            MAF = float(MAF) / (2 * nsamples)
            # print "E = %s MAF = %s" % (efficacy, MAF)
        if (efficacy >= eff_cutoff and MAF >= MAF_cutoff):
            pass
        line = f.readline()
        if (not line):
            break

if (not sample_ids is None):
    if (sample_ids.shape[0] != nsamples + 1):
        sys.stderr.write("warning: the number of samples in the gen file and the sample files differ! (genfile: %s sample file: %s\nThe sample file will not be used!\n" % (nsamples, sample_ids.shape[0] - 1))
        sample_ids = None

# finally write the output
outfile.write("id\tPRS\n")

for i in range(nsamples):
    sid = i
    if (not sample_ids is None):
        sid = sample_ids.loc[i + 1, "ID_1"]
    outfile.write("%s\t%f\n" % (sid, PRS[i]))

