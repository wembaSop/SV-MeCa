#!/usr/bin/env python3

## Author: Rudel Christian Nkouamedjo Fankep
## date : 10.11.2022
"""
This script take the ctx output of breakdancer and parse it to a vcf file 

"""

from datetime import datetime
from os import path
import argparse
import time

def my_parser():
    parser = argparse.ArgumentParser(description="Programm to transform breakdancer output(.ctx) in to a VCF(v4.2)")
    parser.add_argument("reference", metavar="FASTA", type=str, help="Reference file (.fa/.fasta) used for the alignment")
    parser.add_argument("breakdancer", type=str, help="output file (.ctx) from brankdancer-max" )
    parser.add_argument("-o","--output", type=str, help="name of the output file (VCF) [default=prefix fom breakdancer CTX output]")
    return parser.parse_args()

def create_output_file(name,content):
    with open(name,"w") as f:
        for line in content:
            f.write(line+"\n")
        
def read_ctx_file(ctx):
    
    with open(ctx,"r") as f:
        content = f.readlines()
    return content

def make_metadata(fasta):
    pair={}
    with open(fasta,"r") as f:
        flag = False
        for line in f.readlines():
            if line.startswith(">"):
                flag = True
                ident=line.strip().split()[0][1:]
                pair[ident]=0
            else:
                if flag:
                    pair[ident]+=len(line.strip())
                else:
                    continue
    contigs = [(i,pair[i]) for i in pair.keys()]    
    mandatory_line = [
        "##fileformat=VCFv4.2",
        "##fileDate={}".format(datetime.today().strftime('%Y%m%d')),
        "##software=breakdancer-max",
        "##source={}".format(path.basename(__file__)),
        "##reference={}".format(path.basename(fasta))
        ]
    ## Add INFO currently CHROM2 END SVLEN SVTYPE
    mandatory_line.append('##INFO=<ID=CHROM2,Number=1,Type=String,Description="Identifier of the chromosome where the SV end">')
    mandatory_line.append('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the SV">')
    mandatory_line.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">')
    mandatory_line.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">')
    ## add FORMAT currently PR, GT wil be ./.
    mandatory_line.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    mandatory_line.append('##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Total number of supporting read pairs>')
    ## Add Contigs informations
    for i in contigs:
        mandatory_line.append("##contig=<ID={},length={}>".format(i[0],i[1]))
    ## add columns header
    mandatory_line.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown')

    return mandatory_line



def main():
    st = time.time()
    args = my_parser()
    ctx_file = read_ctx_file(args.breakdancer)
    
    header = make_metadata(args.reference)
    for line in ctx_file:
        if line.startswith("#"):
            continue
        else:
            line = line.split()
            header.append("{chrom}\t{pos}\t{id}\t{ref}\t<{alt}>\t{qual}\t{filter}\tCHROM2={chrom2};END={end};SVTYPE={svtype};SVLEN={svlen}\tGT:PR\t{gt}:{dp}".format(
                chrom=line[0],
                pos=line[1],
                id=".",
                ref=".",
                gt="./.",
                filter=".",
                chrom2=line[3],
                end=line[4],
                alt=line[6],
                svtype=line[6],
                svlen=abs(int(line[7])),
                qual=line[8],
                dp=line[9]))
    if args.output:
        output_name = args.output
    else:
        output_name = path.basename(args.breakdancer).replace("ctx","vcf")

    create_output_file(output_name,header)

    ##print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(time.time()-st)))
    ##print("done!")



if __name__ == "__main__":
    main()

