#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='Phase allelic reads by chromsome names')
parser.add_argument('--in', type=str, dest="inpairs", help='input pairs.gz')
parser.add_argument('--out', type=str, dest="outpfx", help='output prefix')
parser.add_argument('--phase1', type=str, dest="phase1", help='phase1 suffix')
parser.add_argument('--phase2', type=str, dest="phase2", help='phase2 suffix')

args = parser.parse_args()
inf = args.inpairs
outp = args.outpfx
nphase1 = args.phase1
nphase2 = args.phase2

import os
import numpy as np
import gzip
from time import perf_counter as pc

def run():
    start_time = pc()
    print("Phasing pairs file...")
    split_pairs(inf, outp, nphase1, nphase2)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def split_pairs(inf, outp, nphase1, nphase2):
    ### v3: directly write to output without storing in memory
    outbfname = outp + "_allele1.pairs"
    pairs1 = open(outbfname, "a")
    outbfname = outp + "_allele2.pairs"
    pairs2 = open(outbfname, "a")
    outbfname = outp + "_mix.pairs"
    mix_pairs = open(outbfname, "a")
    outbfname = outp + "_bad.pairs"
    bad_pairs = open(outbfname, "a")
    phase_suffixes = [nphase1, nphase2]
    phase_suffix_lens = [len(suffix) for suffix in phase_suffixes]

    with gzip.open(inf, 'rt') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#"):
                if line.startswith("#columns"):
                    fields = line.split(' ')
                    fields.extend(["phase1", "phase2"])
                    line = ' '.join(fields)
                pairs1.write(line + '\n')
                pairs2.write(line+ '\n')
                mix_pairs.write(line + '\n')
                bad_pairs.write(line + '\n')
                continue
            
            fields = line.split('\t')
            readID, chrom1, pos1, chrom2, pos2, strand1, strand2, pair_type, XA1, XA2, NM1, NM2, AS1, AS2, XS1, XS2 = fields ### XA mode
            # readID, chrom1, pos1, chrom2, pos2, strand1, strand2, pair_type, XB1, XB2, AS1, AS2, XS1, XS2 = line.strip().split('\t') ### XB mode, not supported for now

            ### start phasing
            lphase, lchrom, lpos = phase_side_XA(chrom1, pos1, XA1, AS1, XS1, NM1, phase_suffixes, phase_suffix_lens)
            rphase, rchrom, rpos = phase_side_XA(chrom2, pos2, XA2, AS2, XS2, NM2, phase_suffixes, phase_suffix_lens)
            fields.extend([lphase, rphase])
            
            ### determine bad pairs:
            if (lphase == "!") or (rphase == "!"): 
                bad_pairs.write('\t'.join(fields) + '\n')
            else: 
            ### both left right can be mapped
                lpfx = [qry.split("_")[1] for qry in lchrom] 
                rpfx = [qry.split("_")[1] for qry in rchrom]

                ### determined paired chrom and pos:
                dpfx = compare_prefix(lpfx, rpfx) 

                if (lphase == ".") and (rphase == '.'): 
                ### non specific: write each entries twice?
                    mix_pairs.write('\t'.join(fields) + '\n')
                elif (lphase != ".") and (rphase != '.'): 
                ###  both left and right have specific mapping to one allele
                    if lphase == "0" and rphase == "0":
                        pairs1.write('\t'.join(fields) + '\n')
                    elif lphase == "1" and rphase == "1":
                        pairs2.write('\t'.join(fields) + '\n')
                    else:
                        bad_pairs.write('\t'.join(fields) + '\n')
                elif (lphase != ".") and (rphase == '.'):
                ### only left has specific mapping
                    pair = dpfx[0]
                    fields[1] = lchrom[pair[0]]
                    fields[5], fields[2] = get_sign_value(lpos[pair[0]])
                    fields[3] = rchrom[pair[1]]
                    fields[6], fields[4] = get_sign_value(rpos[pair[1]])
                    if lphase == "0":
                        pairs1.write('\t'.join(fields) + '\n')
                    elif lphase == "1":
                        pairs2.write('\t'.join(fields) + '\n')
                elif (lphase == ".") and (rphase != '.'): 
                ### only right has specific mapping
                    pair = dpfx[0]
                    fields[1] = lchrom[pair[0]]
                    fields[5], fields[2] = get_sign_value(lpos[pair[0]])
                    fields[3] = rchrom[pair[1]]
                    fields[6], fields[4] = get_sign_value(rpos[pair[1]])
                    if rphase == "0":
                        pairs1.write('\t'.join(fields) + '\n')
                    elif rphase == "1":
                        pairs2.write('\t'.join(fields) + '\n')

    pairs1.close()
    pairs2.close()
    mix_pairs.close()
    bad_pairs.close()

def get_chrom_phase(chrom, phase_suffixes, phase_suffix_lens):
    for suffix, length in zip(phase_suffixes, phase_suffix_lens):
        if chrom.endswith(suffix):
            return str(phase_suffixes.index(suffix)), chrom, chrom[: -length]
    return "!", chrom, chrom

### to del with Ming's data, we need to consider the chrom name when 
### selecting mapped chromsomes without SNP to prevent extend beyond chrom end / start
def phase_side_XA(chrom, pos, XA, AS, XS, NM, phase_suffixes, phase_suffix_lens):
    phase, chrom_all, chrom_base = get_chrom_phase(chrom, phase_suffixes, phase_suffix_lens)
    XAs = [i for i in XA.split(";") if len(i.strip()) > 0]
    if len(XAs) >= 1:
        alt_chrom, alt_pos, alt_CIGAR, alt_NM = XAs[0].split(",")
        M1 = pos
        M2 = alt_pos
    else:
        M1 = pos
        M2 = pos
    if AS > XS:  
    ### Primary hit has higher score than the secondary
        return phase, [chrom_all, chrom_all], [M1, M2]
    elif len(XAs) >= 1:
        if len(XAs) >= 2:
            alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM = XAs[1].split(",")
            M3 = alt2_pos
            if int(alt2_NM) == int(alt_NM) == NM:
                return "!", ["!", "!"], [M1, M2]

        alt_phase, alt_chrom_all, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes, phase_suffix_lens)
        alt_is_homologue = (chrom_base == alt_chrom_base) and (
            ((phase == "0") and (alt_phase == "1"))
            or ((phase == "1") and (alt_phase == "0"))
        )

        if alt_is_homologue: 
            return ".", [chrom_all, alt_chrom_all], [M1, M2]

    return "!", ["!", "!"], [M1, M2]

def compare_prefix(array1, array2):
    return [(i, j) for i, a in enumerate(array1) for j, b in enumerate(array2) if a == b]

def get_sign_value(number):
    number = int(number)
    sign = "+" if number >= 0 else "-"
    absolute_value = str(abs(number))
    return sign, absolute_value

if __name__ == "__main__":
    run()


