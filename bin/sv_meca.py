#!/usr/bin/env python3

## Authors: Rudel Christian Nkouamedjo Fankep, Corinna Ernst
## date : 16.04.2024

import pandas as pd
import pickle
import xgboost as xgb
import argparse
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description='Your script description here')

    # Mandatory arguments
    parser.add_argument('model_ins', type=str, help='pickle model for insertion')
    parser.add_argument('model_del', type=str, help='pickle model for deletions')

    # Optional arguments with default values
    parser.add_argument('--stats', '-s', type=str, required=True, help='txt file with some informations about the data (coverage, read length)')
    parser.add_argument('--breakdancer', '-bd', type=str, required=True, help='Breakdancer VCF')
    parser.add_argument('--delly', '-dl', type=str, required=True, help='Delly VCF')
    parser.add_argument('--insurveyor', '-is', type=str, required=True, help='INSurVeyor VCF')
    parser.add_argument('--lumpy', '-lp', type=str, required=True, help='Lumpy VCF')
    parser.add_argument('--manta', '-mt', type=str, required=True, help='Manta VCF')
    parser.add_argument('--pindel', '-pd', type=str, required=True, help='Pindel VCF')
    parser.add_argument('--tardis', '-td', type=str, required=True, help='TARDIS VCF')
    parser.add_argument('--survivor', '-su', type=str, required=True, help='SURVIVOR VCF')
    parser.add_argument('--output', '-o', type=str, required=True, help='Output file name')
    args = parser.parse_args()
    return args

def parse_stats(filename):
    infos={}
    with open(filename) as f:
        for line in f.readlines():
            tmp_line = line.strip().split("=")
            infos[tmp_line[0]] = tmp_line[1]

    return int(float(infos["coverage"])), int(float(infos["readlen"]))


def main(my_args):
    COLS = ["SURVIVOR_ID", "CHROM", "POS", "SVLEN", "END", "SVTYPE", "BD_ID", "DELLY_ID", "INSV_ID", "LUMPY_ID", "MANTA_ID", "PINDEL_ID", "TARDIS_ID"]
    CONTIGS=[]
    FINAL_SET = pd.DataFrame(columns=COLS)
    FINAL_SET["SVLEN"] = 0
    FINAL_SET["SVTYPE"] = ""
    MCOV , RLEN = parse_stats(args.stats) # RLEN not needed for DEL

    print("... parsing survivor")

    with __builtins__.open(my_args.survivor,'r') as infile:
        for line in infile:
    #        print(line)#
            if 'contig=<ID=' in line:
                CONTIGS.append(line.strip())
            if not line.startswith('#') and not line.startswith("X") and not line.startswith("Y"):
                ll = line.split('\t')
                ll[0]=ll[0].replace("chr","")
                SVTYPE = [_ for _ in ll[7].split(';') if _.startswith("SVTYPE")][0].split("=")[1]
                if SVTYPE in ["DEL", "INS"] and ll[0].isdigit():
                    DATA_LINE = {"SURVIVOR_ID" : ll[2], "SVTYPE": SVTYPE, "CHROM": int(ll[0]), "POS": int(ll[1]), "REF":str(ll[3]),"ALT":str(ll[4])}
                    END = abs(int([_ for _ in ll[7].split(';') if _.startswith("END")][0].split("=")[1]))
                    SVLEN = abs(int([_ for _ in ll[7].split(';') if _.startswith("SVLEN")][0].split("=")[1]))
                    DATA_LINE["SVLEN"] = SVLEN
                    DATA_LINE["END"] = END


                    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BD_HG00514	DL_HG00514	HG00514	LP_HG00514	MT_HG00514	PD_HG00514	TD_HG00514
                    # get standalone IDs
                    ID_IND = ll[8].split(':').index('ID')
                    SUPP_VEC = [_ for _ in ll[7].split(';') if _.startswith("SUPP_VEC")][0].split("=")[1]

                    if SUPP_VEC[0] == "1":
                        #BreakDancer
                        DATA_LINE["BD_ID"] = ll[9].split(':')[ID_IND]
                    if SUPP_VEC[1] == "1":
                        #Delly
                        DATA_LINE["DELLY_ID"] = ll[10].split(':')[ID_IND]
                    if SUPP_VEC[2] == "1":
                        #InSurveyor
                        DATA_LINE["INSV_ID"] = ll[11].split(':')[ID_IND]
                    if SUPP_VEC[3] == "1":
                        #InSurveyor
                        DATA_LINE["LUMPY_ID"] = ll[12].split(':')[ID_IND]
                    if SUPP_VEC[4] == "1":
                        #Manta
                        DATA_LINE["MANTA_ID"] = ll[13].split(':')[ID_IND].replace("_", ":")
                    if SUPP_VEC[5] == "1":
                        #Pindel
                        DATA_LINE["PINDEL_ID"] = ll[14].split(':')[ID_IND]
                    if SUPP_VEC[6] == "1":
                        #TARDIS
                        DATA_LINE["TARDIS_ID"] = ll[15].split(':')[ID_IND]
                    DATA_LINE_DF = pd.DataFrame([DATA_LINE]) 
                    FINAL_SET = pd.concat([FINAL_SET, DATA_LINE_DF], ignore_index=True)
                    
    FINAL_SET = FINAL_SET.sort_values(by=['CHROM', 'POS'])

    #DEBUG
    #print(FINAL_SET)
    #FINAL_SET.to_csv("hg002_52_basic.tsv", sep='\t', index=False)
    #DATA = pd.read_csv("hg002_basic.tsv", sep='\t')

    # DF DEL
    DATA = FINAL_SET.loc[FINAL_SET["SVTYPE"] == "DEL"]
    DATA = DATA[['SURVIVOR_ID', "CHROM", "POS", "REF", "ALT", "END",'SVLEN', 'BD_ID', 'DELLY_ID', 'LUMPY_ID', 'MANTA_ID', 'PINDEL_ID', 'TARDIS_ID']]
    DATA.drop(DATA.loc[DATA.BD_ID.isnull() & DATA.DELLY_ID.isnull() & DATA.LUMPY_ID.isnull() & DATA.MANTA_ID.isnull() & DATA.PINDEL_ID.isnull() & DATA.TARDIS_ID.isnull() ].index, inplace=True)

    # DF INS
    DATA_INS = FINAL_SET.loc[FINAL_SET["SVTYPE"] == "INS"]
    DATA_INS = DATA_INS[['SURVIVOR_ID', "CHROM", "POS", "REF", "ALT", "END", 'SVLEN', 'DELLY_ID', 'INSV_ID', 'MANTA_ID', 'PINDEL_ID', 'TARDIS_ID']]
    DATA_INS.drop(DATA_INS.loc[DATA_INS.DELLY_ID.isnull() & DATA_INS.INSV_ID.isnull() & DATA_INS.MANTA_ID.isnull() & DATA_INS.PINDEL_ID.isnull() & DATA_INS.TARDIS_ID.isnull() ].index, inplace=True)

    print("... parsing BreakDancer")## Only DEL

    with open(my_args.breakdancer) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                ## Only Deletions
                if ll[2] in DATA.BD_ID.values:
                    _ind = ll[8].split(':').index('PR')
                    DATA.loc[DATA['BD_ID'] == ll[2], "BD_RD"] = int(ll[9].split(':')[_ind])/MCOV
                    _ind = ll[8].split(':').index('SBP')
                    DATA.loc[DATA['BD_ID'] == ll[2], "BD_SB1"] = ll[9].split(':')[_ind]
                    _ind = ll[8].split(':').index('SBE')
                    DATA.loc[DATA['BD_ID'] == ll[2], "BD_SB2"] = ll[9].split(':')[_ind]


    #########################
                    
    print("... parsing Delly")## DEL & INS

    with open(my_args.delly) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                ## For Deletions
                if ll[2] in DATA.DELLY_ID.values:
                    #DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_QUAL"] = float(ll[5]) ### XXX convert to normalized rank score!
                    
                    if ll[6] == "PASS":
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_PASS"] = int(1)
                    elif ll[6] == "LowQual":
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_PASS"] = int(0)
                    
                    ### DELLY_VF: RV/(RR+RV) for PRECISE calls; DV/(DR+DV) for IMPRECISE calls
                    ### DELLY_HQSUPPORT: Estimator of local high quality read depth: RR+RV for PRECISE calls, DR+DV for IMPRECISE calls, normalized by mean coverage
                    if "IMPRECISE" in ll[7]:
                        # IMPRECISE
                        _ind = ll[8].split(':').index('DR')
                        DR = int(ll[9].split(':')[_ind])
                        _ind = ll[8].split(':').index('DV')
                        DV = int(ll[9].split(':')[_ind])
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_HQSUPPORT"] = DR+DV
                    else:
                        #PRECISE
                        _ind = ll[8].split(':').index('RV')
                        RV = int(ll[9].split(':')[_ind])
                        _ind = ll[8].split(':').index('RR')
                        RR = int(ll[9].split(':')[_ind])
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_HQSUPPORT"] = RR+RV
                    
                    #DELLY_PE	Delly (VCF)	0-Inf	yes	Count of supporting paired reads (PE), normalized by mean coverage
                    PE = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("PE=")][0]
                    DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_PE"] = int(PE)/MCOV

                    #DELLY_SR	Delly (VCF)	0-Inf	yes	Count of supporting split reads (SR), normalized by mean coverage
                    if "SR=" in ll[7]:
                        SR = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("SR=")][0]
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_SR"] = int(SR)/MCOV

                    if not "IMPRECISE" in ll[7]:
                        SRQ = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("SRQ=")][0]
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_SRQ"] = SRQ
                        CE = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("CE=")][0]
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_CE"] = CE
                    
                    if "HOMLEN" in ll[7]:
                        homlen = [_ for _ in ll[7].split(';') if _.startswith("HOMLEN=")][0].split("=")[1]
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_HOMLEN"] = int(homlen)
                    else:
                        DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_HOMLEN"] = 0
                    
                    #DELLY_GQ	Delly (VCF)	0-10000	yes	Genotype quality
                    gq_ind = ll[8].split(':').index('GQ')
                    DATA.loc[DATA['DELLY_ID'] == ll[2], "DELLY_GQ"] = ll[9].split(':')[gq_ind]
                ## For Insertions
                elif ll[2] in DATA_INS.DELLY_ID.values:
                    if ll[6] == "PASS":
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_PASS"] = int(1)
                    elif ll[6] == "LowQual":
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_PASS"] = int(0)

                    #DELLY_MAPQ	Delly (VCF)	0-60	yes	Median mapping quality of paired-end reads
                    MAPQ = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("MAPQ=")][0]
                    DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_MAPQ"] = MAPQ

                    if not "IMPRECISE" in ll[7]:
                        SRQ = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("SRQ=")][0]
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_SRQ"] = SRQ
                        CE = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith("CE=")][0]
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_CE"] = CE

                    #DELLY_FT_PASS	Delly (VCF)	0 or 1	yes	Binary flag indicating that the per-sample genotype filter (FT) was passed
                    ft_ind = ll[8].split(':').index('FT')
                    FT = ll[9].split(':')[ft_ind]
                    if FT == 'PASS':
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_FT_PASS"] = 1
                    elif FT == "LowQual":
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_FT_PASS"] = 0
                    else:
                        DATA_INS.loc[DATA_INS['DELLY_ID'] == ll[2], "DELLY_FT_PASS"] = ""

    #DATA.to_csv("hg00733_basic_second.tsv", sep='\t', index=False)
    #########################

    print("... parsing InSurVeyor") #Only INS


    with __builtins__.open(my_args.insurveyor) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                if ll[2] in DATA_INS["INSV_ID"].values:

                    # INSV_ANOMALOUS_DEPTH
                    if "ANOMALOUS_DEPTH" in ll[6]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_ANOMALOUS_DEPTH"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_ANOMALOUS_DEPTH"] = 0
                    
                    # INSV_ALT_SHORTER_THAN_REF
                    if "ALT_SHORTER_THAN_REF" in ll[6]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_ALT_SHORTER_THAN_REF"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_ALT_SHORTER_THAN_REF"] = 0
                    
                    # INSV_LOW_SUPPORT
                    if "LOW_SUPPORT" in ll[6]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_LOW_SUPPORT"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_LOW_SUPPORT"] = 0

                    # NO_DISC_SUPPORT
                    if "NO_DISC_SUPPORT" in ll[6]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_NO_DISC_SUPPORT"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_NO_DISC_SUPPORT"] = 0
                    
                    # HOMOPOLYMER_INSSEQ
                    if "HOMOPOLYMER_INSSEQ" in ll[6]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_HOMOPOLYMER_INSSEQ"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_HOMOPOLYMER_INSSEQ"] = 0

                    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles."> #XXX considered in global SVLEN?
                    ##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description="Length of the inserted sequence."> 
                    #COLS = COLS + ["INSV_DIFF_SVLEN"] # XXX Difference between SVINSLEN and SVLEN
                    SVLEN = int([_ for _ in ll[7].split(';') if _.startswith("SVLEN=")][0].split('=')[1])
                    SVINSLEN = int([_ for _ in ll[7].split(';') if _.startswith("SVINSLEN=")][0].split('=')[1])
                    DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_DIFF_SVLEN"] = SVINSLEN - SVLEN #XXX always >=0

                    ##INFO=<ID=DISCORDANT,Number=2,Type=Integer,Description="Discordant pairs supporting the left and right breakpoints of this insertion."> ### TODO has to be normalized , LEFT & RIGHT
                    #COLS = COLS + ["INSV_DSCORDANT_LEFT", "INSV_DSCORDANT_RIGHT"]
                    if "DISCORDANT=" in ll[7]:
                        DISCORDANT = [_ for _ in ll[7].split(';') if _.startswith("DISCORDANT=")][0]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_DISCORDANT_LEFT"] = int(DISCORDANT.split("=")[1].split(',')[0])/MCOV
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_DISCORDANT_RIGHT"] = int(DISCORDANT.split("=")[1].split(',')[1])/MCOV

                    ##INFO=<ID=SPLIT_READS,Number=2,Type=Integer,Description="Split reads supporting the left and right breakpoints of this insertion.">
                    ##INFO=<ID=FWD_SPLIT_READS,Number=2,Type=Integer,Description="Forward split reads supporting the left and right breakpoints of this insertion.">
                    ##INFO=<ID=REV_SPLIT_READS,Number=2,Type=Integer,Description="Reverse split reads supporting the left and right breakpoints of this insertion.">
                    # XXX SPLIT_READS is given by FWD and REV values
                    # XXX use (normalized) SPLIT_READS left and right &
                    # XXX use ratios (strand bias) left and right (by dividing lower count by the higer one): 
                    #       min(FWD_SPLIT_READS[0], REV_SPLIT_READS[0])/max(FWD_SPLIT_READS[0], REV_SPLIT_READS[0])
                    #       min(FWD_SPLIT_READS[1], REV_SPLIT_READS[1])/max(FWD_SPLIT_READS[1], REV_SPLIT_READS[1])
                    #COLS = COLS + ["INSV_SPLIT_READS_LEFT", "INSV_SPLIT_READS_RIGHT"]
                    #COLS = COLS + ["INSV_SPLIT_RATIO_LEFT", "INSV_SPLIT_RATIO_RIGHT"]
                    if "SPLIT_READS" in ll[7]:
                        SPLIT_READS = [_ for _ in ll[7].split(';') if _.startswith("SPLIT_READS=")][0].split('=')[1]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_SPLIT_READS_LEFT"] = int(SPLIT_READS.split(',')[0])/MCOV

                        FWD_SPLIT_READS = [_ for _ in ll[7].split(';') if _.startswith("FWD_SPLIT_READS=")][0].split('=')[1]
                        FWD_SPLIT_READS = FWD_SPLIT_READS.split(',')
                        FWD_SPLIT_READS = [int(FWD_SPLIT_READS[0]), int(FWD_SPLIT_READS[1])]
                        REV_SPLIT_READS = [_ for _ in ll[7].split(';') if _.startswith("REV_SPLIT_READS=")][0].split('=')[1]
                        REV_SPLIT_READS = REV_SPLIT_READS.split(',')
                        REV_SPLIT_READS = [int(REV_SPLIT_READS[0]), int(REV_SPLIT_READS[1])]
                        
                        #FWD_SPLIT_READS=0,40;REV_SPLIT_READS=0,4
                        if FWD_SPLIT_READS[0] + REV_SPLIT_READS[0]:
                            DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_SPLIT_RATIO_LEFT"] = min(FWD_SPLIT_READS[0], REV_SPLIT_READS[0])/max(FWD_SPLIT_READS[0], REV_SPLIT_READS[0])
                        else:
                            DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_SPLIT_RATIO_LEFT"] = 0
                            
                    ##INFO=<ID=SPANNING_READS,Number=2,Type=Integer,Description="Negative evidence."> # XXX only A_ and T_INS
                    ## XXX normalize use left and right
                    #COLS = COLS + ["INSV_SPANNING_LEFT", "INSV_SPANNING_RIGHT"]
                    if "SPANNING_READS" in ll[7]:
                        SPANNING_READS = [_ for _ in ll[7].split(';') if _.startswith("SPANNING_READS=")][0].split('=')[1]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_SPANNING_RIGHT"] = int(SPANNING_READS.split(',')[1])/MCOV

                    ##INFO=<ID=TRANS_QUERY_COV,Number=2,Type=Integer,Description="Length of the prefix and suffix of the transposed sequence that was covered by reads."> #XXX use unnormalized , 2 values
                    #COLS = COLS + ["INSV_TRANS_QUERY_COV_LEFT", "INSV_TRANS_QUERY_COV_RIGHT"]
                    if "TRANS_QUERY_COV" in ll[7]:
                        TRANS_QUERY_COV = [_ for _ in ll[7].split(';') if _.startswith("TRANS_QUERY_COV=")][0].split('=')[1]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_TRANS_QUERY_COV_LEFT"] = int(TRANS_QUERY_COV.split(',')[0])
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_TRANS_QUERY_COV_RIGHT"] = int(TRANS_QUERY_COV.split(',')[1])

                    ##INFO=<ID=STABLE_DEPTHS,Number=2,Type=Integer,Description="Depths of the stable regions (in practice, the regions left and right of the insertion site)."> # XXX normalize, 2 values
                    #COLS = COLS + ["INSV_STABLE_DEPTHS_LEFT", "INSV_STABLE_DEPTHS_RIGHT"]
                    if "STABLE_DEPTHS" in ll[7]:
                        STABLE_DEPTHS = [_ for _ in ll[7].split(';') if _.startswith("STABLE_DEPTHS=")][0].split('=')[1]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_STABLE_DEPTHS_LEFT"] = int(STABLE_DEPTHS.split(',')[0])/MCOV
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_STABLE_DEPTHS_RIGHT"] = int(STABLE_DEPTHS.split(',')[1])/MCOV

                    ##INFO=<ID=AVG_STABLE_NM,Number=2,Type=Integer,Description="Average NM of stable reads."> # XXX ???? only A_INS and T_INS
                    ### COLS = COLS + ["INSV_AVG_STABLE_NM_LEFT", "INSV_AVG_STABLE_NM_RIGHT"]
                    if "AVG_STABLE_NM" in ll[7]:
                        AVG_STABLE_NM = [_ for _ in ll[7].split(';') if _.startswith("AVG_STABLE_NM=")][0].split('=')[1]
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_AVG_STABLE_NM_LEFT"] = int(AVG_STABLE_NM.split(',')[0])/RLEN
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_AVG_STABLE_NM_RIGHT"] = int(AVG_STABLE_NM.split(',')[1])/RLEN

                    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Breakpoints and/or inserted sequence are imprecise."> # XXX only for T_INS 
                    #COLS = COLS + ["INSV_IMPRECISE"]
                    if "IMPRECISE" in  ll[7]:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_IMPRECISE"] = 1
                    else:
                        DATA_INS.loc[DATA_INS['INSV_ID'] == ll[2], "INSV_IMPRECISE"] = 0

    #########################                
    print("... parsing Lumpy") #Only DEL


    with open(my_args.lumpy) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                if int(ll[2]) in DATA.LUMPY_ID.values:
                    ID = int(ll[2]) ## XXX numeric ID for Lumpy!!!
                    
                    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
                    if "IMPRECISE" in ll[7]:
                        DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_PRECISE"] = 0
                    else:
                        DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_PRECISE"] = 1
                    
                    #Genotype quality
                    _ind = ll[8].split(':').index('GQ')
                    DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_GQ"] = ll[9].split(':')[_ind]
                    
                    
                    #LUMPY_RPFR	LUMPY/smoove (VCF)	0-1	yes	no	no	RP/RO; fraction of reference allele PR observation count and reference allele observation count, with partial observations recorded fractionally
                    ro_ind = ll[8].split(':').index('RO')
                    ro = int(ll[9].split(':')[ro_ind])
                    rp_ind = ll[8].split(':').index('RP')
                    rp = int(ll[9].split(':')[rp_ind])
                    try:
                        DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_RPFR"] = rp/ro
                    except:
                        DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_RPFR"] = 0
                        
                    #LUMPY_QAFR	LUMPY/smoove (VCF)	0-1	yes	no	no	QA/AO; fraction of sum of quality of alternate observations and alternate allele observations, with partial observations recorded fractionally
                    #LUMPY_ASCFR	LUMPY/smoove (VCF)	0-1	yes	no	no	ASC/AO; fraction of  alternate allele clipped-read observation count and alternate allele observations, with partial observations recorded fractionally
                    #LUMPY_APFR	LUMPY/smoove (VCF)	0-1	yes	no	no	AP/AO; fraction of alternate allele PR observation count and alternate allele observations, with partial observations recorded fractionally
                    qa_ind = ll[8].split(':').index('QA')
                    qa = int(ll[9].split(':')[qa_ind])
                    ao_ind = ll[8].split(':').index('AO')
                    ao = int(ll[9].split(':')[ao_ind])
                    asc_ind = ll[8].split(':').index('ASC')
                    asc = int(ll[9].split(':')[asc_ind])
                    ap_ind = ll[8].split(':').index('AP')
                    ap = int(ll[9].split(':')[ap_ind])
                    DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_QAFR"] = qa/ao
                    DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_ASCFR"] = asc/ao
                    DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_APFR"] = ap/ao
                    
                    #LUMPY_AB	LUMPY/smoove (VCF)	0-1	yes	no	no	Allele balance, fraction of observations from alternate allele, QA/(QR+QA)
                    qr_ind = ll[8].split(':').index('QR')
                    qr = int(ll[9].split(':')[qr_ind])
                    DATA.loc[DATA['LUMPY_ID'] == ID, "LUMPY_AB"] = qa/(qr+qa)

    #########################
                    
    print("... parsing Manta") # DEL & INS


    with open(my_args.manta) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                # DEL 
                if ll[2] in DATA.MANTA_ID.values:
                    DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_QUAL"] = float(ll[5])/MCOV
                    
                    ##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints">
                    if "HOMLEN=" in ll[7]:
                        HOMLEN = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith('HOMLEN=')][0]
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_HOMLEN"] = float(HOMLEN)
                    else:
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_HOMLEN"] = 0
                    
                    pr_ref, pr_alt = 0, 0
                    if "PR" in ll[8]:
                        _ind = ll[8].split(':').index('PR')
                        pr_ref = int(ll[9].split(':')[_ind].split(',')[0])
                        pr_alt = int(ll[9].split(':')[_ind].split(',')[1])
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_PR"] = float(ll[9].split(':')[_ind].split(',')[1])/MCOV
                    else: DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_PR"] = 0

                    if "SampleFT" in ll[6]:
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_SampleFT"] = 1
                    else:
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_SampleFT"] = 0
                    if "NoPairSupport" in ll[6]:
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_NoPairSupport"] = 1
                    else:
                        DATA.loc[DATA['MANTA_ID'] == ll[2], "MANTA_NoPairSupport"] = 0
                elif ll[2] in DATA_INS.MANTA_ID.values:

                    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
                    if "IMPRECISE" in ll[7]:
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_PRECISE"] = 0
                    else:
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_PRECISE"] = 1
                    
                    ##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints">
                    if "HOMLEN=" in ll[7]:
                        HOMLEN = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith('HOMLEN=')][0]
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_HOMLEN"] = float(HOMLEN)
                    else:
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_HOMLEN"] = 0
                    
                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
                    _ind = ll[8].split(':').index('GQ')
                    DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_GQ"] = ll[9].split(':')[_ind]
                    
                    sr_ref, sr_alt = 0, 0
                    if "SR" in ll[8]:
                        _ind = ll[8].split(':').index('SR')
                        sr_ref = int(ll[9].split(':')[_ind].split(',')[0])
                        sr_alt = int(ll[9].split(':')[_ind].split(',')[1])
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_SR"] = float(ll[9].split(':')[_ind].split(',')[1])/MCOV
                    else: DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_SR"] = 0
                    
                    pr_ref, pr_alt = 0, 0
                    if "PR" in ll[8]:
                        _ind = ll[8].split(':').index('PR')
                        pr_ref = int(ll[9].split(':')[_ind].split(',')[0])
                        pr_alt = int(ll[9].split(':')[_ind].split(',')[1])
                        DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_PR"] = float(ll[9].split(':')[_ind].split(',')[1])/MCOV
                    else: DATA_INS.loc[DATA_INS['MANTA_ID'] == ll[2], "MANTA_PR"] = 0


    #########################
                        
    print("... parsing Pindel") # DEL & INS


    with open(my_args.pindel) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                if ll[2] in DATA.PINDEL_ID.values:
                
                    #PINDEL_DP	Pindel (VCF)	0-Inf	yes	yes	yes	Local coverage estimate: ratio of sum of AD values and mean coverage
                    ad_ind = ll[8].split(':').index('AD')
                    ad_sum = int(ll[9].split(':')[ad_ind].split(',')[0]) + int(ll[9].split(':')[ad_ind].split(',')[1])
                    DATA.loc[DATA['PINDEL_ID'] == ll[2], "PINDEL_DP"] = ad_sum/MCOV
                    
                    #PINDEL_VF	Pindel (VCF)	0-1	yes	yes	yes	Variant fraction estimate: Number of alternate allele supporting reads divided by sum of AD values
                    DATA.loc[DATA['PINDEL_ID'] == ll[2], "PINDEL_VF"] = int(ll[9].split(':')[ad_ind].split(',')[1])/ad_sum
                    
                    #PINDEL_HOMLEN	Pindel (VCF)	0-Inf	yes	yes	yes	Length of base pair identical micro-homology at event breakpoints, set to zero if missing
                    if "HOMLEN" in ll[7]:
                        HOMLEN = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith('HOMLEN=')][0]
                        DATA.loc[DATA['PINDEL_ID'] == ll[2], "PINDEL_HOMLEN"] = HOMLEN
                    else:
                        DATA.loc[DATA['PINDEL_ID'] == ll[2], "PINDEL_HOMLEN"] = 0
                elif ll[2] in DATA_INS.PINDEL_ID.values:
            
                    #PINDEL_DP	Pindel (VCF)	0-Inf	yes	yes	yes	Local coverage estimate: ratio of sum of AD values and mean coverage
                    

                    ad_ind = ll[8].split(':').index('AD')
                    ad_sum = int(ll[9].split(':')[ad_ind].split(',')[0]) + int(ll[9].split(':')[ad_ind].split(',')[1])
                    DATA_INS.loc[DATA_INS['PINDEL_ID'] == ll[2], "PINDEL_DP"] = ad_sum/MCOV
                    
                    #PINDEL_VF	Pindel (VCF)	0-1	yes	yes	yes	Variant fraction estimate: Number of alternate allele supporting reads divided by sum of AD values
                    DATA_INS.loc[DATA_INS['PINDEL_ID'] == ll[2], "PINDEL_VF"] = int(ll[9].split(':')[ad_ind].split(',')[1])/ad_sum
                    
                    #PINDEL_HOMLEN	Pindel (VCF)	0-Inf	yes	yes	yes	Length of base pair identical micro-homology at event breakpoints, set to zero if missing
                    if "HOMLEN" in ll[7]:
                        HOMLEN = [_.split('=')[1] for _ in ll[7].split(';') if _.startswith('HOMLEN=')][0]
                        DATA_INS.loc[DATA_INS['PINDEL_ID'] == ll[2], "PINDEL_HOMLEN"] = HOMLEN
                    else:
                        DATA_INS.loc[DATA_INS['PINDEL_ID'] == ll[2], "PINDEL_HOMLEN"] = 0

    #########################
                        
    print("... parsing TARDIS") # DEL & INS 


    with open(my_args.tardis) as infile:
        for line in infile:
            if not line.startswith('#'):
                ll = line.rstrip().split('\t')
                if ll[2] in DATA.TARDIS_ID.values:
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                    gt = ll[9].split(':')[0]
                    if gt in ['0/0', './.']:
                        DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_GT"] = 0
                    elif gt in ['0/1', '1/0']:
                        DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_GT"] = 0.5
                    elif gt == '1/1':
                        DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_GT"] = 1
                    
                    ##FORMAT=<ID=CNVL,Number=1,Type=Float,Description="CNV Likelihood">
                    _ind = ll[8].split(':').index('CNVL')
                    DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_CNVL"] = float(ll[9].split(':')[_ind])
                    
                    ##INFO=<ID=RPSUP,Number=1,Type=Integer,Description="Number of supporting read pairs">
                    ##INFO=<ID=SRSUP,Number=1,Type=Integer,Description="Number of supporting split reads">
                    ##FORMAT=<ID=RP,Number=1,Type=Integer,Description="Read Pair Support">
                    ##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Split Read Support">
                    # ...identical for single sample runs
                    _ind = ll[8].split(':').index('RP')
                    DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_RPSUP"] = float(ll[9].split(':')[_ind])/MCOV
                    _ind = ll[8].split(':').index('SR')
                    DATA.loc[DATA['TARDIS_ID'] == ll[2], "TARDIS_SRSUP"] = float(ll[9].split(':')[_ind])/MCOV
                elif ll[2] in DATA_INS.TARDIS_ID.values:
                    ##FORMAT=<ID=CNVL,Number=1,Type=Float,Description="CNV Likelihood">
                    _ind = ll[8].split(':').index('CNVL')
                    DATA_INS.loc[DATA_INS['TARDIS_ID'] == ll[2], "TARDIS_CNVL"] = float(ll[9].split(':')[_ind])
                    
                    ##INFO=<ID=RPSUP,Number=1,Type=Integer,Description="Number of supporting read pairs">
                    ##INFO=<ID=SRSUP,Number=1,Type=Integer,Description="Number of supporting split reads">
                    ##FORMAT=<ID=RP,Number=1,Type=Integer,Description="Read Pair Support">
                    ##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Split Read Support">
                    # ...identical for single sample runs
                    _ind = ll[8].split(':').index('SR')
                    DATA_INS.loc[DATA_INS['TARDIS_ID'] == ll[2], "TARDIS_SRSUP"] = float(ll[9].split(':')[_ind])/MCOV

    ### XXX ensure feature order DEL
    BASE_COLS = ["SURVIVOR_ID", "CHROM", "POS", "REF", "ALT", "END", "SVLEN", "BD_ID", "DELLY_ID", "LUMPY_ID", "MANTA_ID", "PINDEL_ID", "TARDIS_ID"]
    BD_COLS = ["BD_RD", "BD_SB1", "BD_SB2"]
    DELLY_COLS = ["DELLY_PASS", "DELLY_HQSUPPORT", "DELLY_PE", "DELLY_SR", "DELLY_SRQ", "DELLY_CE", "DELLY_HOMLEN", "DELLY_GQ"] 
    LUMPY_COLS = ["LUMPY_PRECISE", "LUMPY_GQ", "LUMPY_RPFR", "LUMPY_QAFR", "LUMPY_ASCFR", "LUMPY_APFR", "LUMPY_AB"]
    MANTA_COLS = ["MANTA_QUAL", "MANTA_HOMLEN", "MANTA_PR", "MANTA_SampleFT", "MANTA_NoPairSupport"]
    PINDEL_COLS = ["PINDEL_DP", "PINDEL_VF", "PINDEL_HOMLEN"]
    TARDIS_COLS = ["TARDIS_GT", "TARDIS_CNVL", "TARDIS_RPSUP", "TARDIS_SRSUP"]
    print(len(BD_COLS), len(DELLY_COLS), len(LUMPY_COLS), len(MANTA_COLS), len(PINDEL_COLS), len(TARDIS_COLS)) ### DEBUG
    print("Sum of Covariates DEL:",sum([len(BD_COLS), len(DELLY_COLS), len(LUMPY_COLS), len(MANTA_COLS), len(PINDEL_COLS), len(TARDIS_COLS)]))
    col_order = BASE_COLS + BD_COLS + DELLY_COLS + LUMPY_COLS + MANTA_COLS + PINDEL_COLS + TARDIS_COLS
    # convert to basic type to avoid object passing to xgboost
    for column in DATA.columns[12:]:
        DATA[column] = pd.to_numeric(DATA[column], errors='coerce')

    DATA = DATA.reindex(columns=col_order)
    #DATA.to_csv("hg002_52_del.tsv", sep='\t', index=False)

    ### XXX ensure feature order INS
    BASE_COLS_INS = ["SURVIVOR_ID", "CHROM", "POS", "REF", "ALT", "END", "SVLEN", "DELLY_ID", "INSV_ID", "MANTA_ID", "PINDEL_ID", "TARDIS_ID"]
    DELLY_COLS_INS = ["DELLY_PASS", "DELLY_MAPQ", "DELLY_SRQ", "DELLY_CE", "DELLY_FT_PASS"] 
    INSV_COLS = ["INSV_ANOMALOUS_DEPTH", "INSV_ALT_SHORTER_THAN_REF", "INSV_LOW_SUPPORT", "INSV_NO_DISC_SUPPORT", "INSV_HOMOPOLYMER_INSSEQ", "INSV_DIFF_SVLEN"]
    INSV_COLS = INSV_COLS + ["INSV_DISCORDANT_LEFT", "INSV_DISCORDANT_RIGHT", "INSV_SPLIT_READS_LEFT", "INSV_SPLIT_RATIO_LEFT"]
    INSV_COLS = INSV_COLS + ["INSV_SPANNING_RIGHT", "INSV_TRANS_QUERY_COV_LEFT", "INSV_TRANS_QUERY_COV_RIGHT"]
    INSV_COLS = INSV_COLS + ["INSV_STABLE_DEPTHS_LEFT", "INSV_STABLE_DEPTHS_RIGHT", "INSV_AVG_STABLE_NM_LEFT", "INSV_AVG_STABLE_NM_RIGHT", "INSV_IMPRECISE"]
    MANTA_COLS_INS = ["MANTA_PRECISE", "MANTA_HOMLEN", "MANTA_GQ", "MANTA_SR", "MANTA_PR"]
    PINDEL_COLS_INS = ["PINDEL_DP", "PINDEL_VF", "PINDEL_HOMLEN"]
    TARDIS_COLS_INS = ["TARDIS_CNVL", "TARDIS_SRSUP"]
    #DATA.to_csv("hg00733_basic_second.tsv", sep='\t', index=False)
    print(len(DELLY_COLS_INS), len(INSV_COLS), len(MANTA_COLS_INS), len(PINDEL_COLS_INS), len(TARDIS_COLS_INS)) ### DEBUG
    print("Sum of Covariates INS:",sum([len(DELLY_COLS_INS), len(INSV_COLS), len(MANTA_COLS_INS), len(PINDEL_COLS_INS), len(TARDIS_COLS_INS)]))
    col_order = BASE_COLS_INS + DELLY_COLS_INS + INSV_COLS + MANTA_COLS_INS + PINDEL_COLS_INS + TARDIS_COLS_INS
    # convert to basic type to avoid object passing to xgboost
    for column in DATA_INS.columns[11:]:
        DATA_INS[column] = pd.to_numeric(DATA_INS[column], errors='coerce')

    DATA_INS = DATA_INS.reindex(columns=col_order)
    #DATA_INS.to_csv("hg002_52_ins.tsv", sep='\t', index=False)
    #DATA.to_csv("hg002_ins_full.tsv", sep='\t', index=False)


    # Model DEL
    xgbc = pickle.load(open(my_args.model_del, "rb"))

    OUT_COLS = ['SURVIVOR_ID', "CHROM", "POS", "REF", "ALT", "END", 'BD_ID', 'DELLY_ID', 'LUMPY_ID', 'MANTA_ID', 'PINDEL_ID', 'TARDIS_ID']

    #xgbc.predict(DATA.loc[:, ~DATA.columns.isin(OUT_COLS)])

    #df = pd.DataFrame(zip(xgbc.get_booster().feature_names, xgbc.feature_importances_), columns=['feature', 'importance'])
    #print(df.sort_values(by=['importance'], ascending=False))

    PREDS = xgbc.predict(DATA.loc[:, ~DATA.columns.isin(OUT_COLS)])
    PROBS = xgbc.predict_proba(DATA.loc[:, ~DATA.columns.isin(OUT_COLS)])
    DATA["Prediction"] = PREDS
    #print(PROBS)
    DATA["Probs"] = PROBS[:,1]
    DATA.to_csv("hg002_52_del.pred.tsv", sep='\t', index=False)


    #Model INS
    xgbc = pickle.load(open(my_args.model_ins, "rb"))

    OUT_COLS = ["SURVIVOR_ID", "CHROM", "POS", "REF", "ALT", "END", "DELLY_ID", "INSV_ID", "MANTA_ID", "PINDEL_ID", "TARDIS_ID"]

    #xgbc.predict(DATA_INS.loc[:, ~DATA_INS.columns.isin(OUT_COLS)])

    #df = pd.DataFrame(zip(xgbc.get_booster().feature_names, xgbc.feature_importances_), columns=['feature', 'importance'])
    #print(df.sort_values(by=['importance'], ascending=False))

    PREDS_INS = xgbc.predict(DATA_INS.loc[:, ~DATA_INS.columns.isin(OUT_COLS)])
    PROBS_INS = xgbc.predict_proba(DATA_INS.loc[:, ~DATA_INS.columns.isin(OUT_COLS)])
    DATA_INS["Prediction"] = PREDS_INS
    #print(PROBS)
    DATA_INS["Probs"] = PROBS_INS[:,1]
    #DATA_INS.to_csv("hg002_52_ins.pred.tsv", sep='\t', index=False)
    return DATA, DATA_INS, CONTIGS

def create_header(the_infile, the_outname, the_contigs):
    the_outfile = open(the_outname,"w")
    header=["##fileformat=VCFv4.1", 
            "##source=SV-MeCav1.0", 
            f'##fileDate={str(datetime.now()).split()[0].replace("-","")}',
            "##FILTER=<ID=LowProb,Description=\"Low probability to represent a real Structural Variant (SV) based on quality metrics from SV callers\">",
            "##INFO=<ID=IDS,Number=.,Type=String,Description=\"IDs of initial variant calls from SV callers supporting the variant. See SUPP_VEC relate about the order\">",
            "##INFO=<ID=SUPP,Number=1,Type=Integer,Description=\"Number of standalone SV callers supporting the variant call\">",
            "##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Vector of supporting SV callers in the following order: BreakDancer, Delly, InSurVeyor, Lumpy, Manta, Pindel, TARDIS. 1=support the variant, 0=do not support.\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">",
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">",
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">"
            ]
    for line in header:
        the_outfile.write(line+"\n")
    for line in the_contigs:
        the_outfile.write(line+"\n")
    the_outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"+"\n")
    return the_outfile

def build_vcf(vcf_file, my_df, svtype, columns, filter):
    # Iterate over each row in the DataFrame
    for index, row in my_df.iterrows():
        # Extract the selected columns for the current row
        #extracted_values = [str(row[column]) for column in columns]
        idx=columns.index("BD_ID")
        SUPP=0
        SUPP_VEC=""
        IDS=""

        for column in columns[idx:-1]:
            try:
                if str(row[column])=="nan":
                    SUPP_VEC+="0"
                else:
                    SUPP+=1
                    SUPP_VEC+="1"
                    IDS+=row[column]+","
            except KeyError:
                SUPP_VEC+="0"
        # remove "," at the end of IDS
        if IDS[-1]==",":
            IDS=IDS[:-1]
        line=[str(row["CHROM"]),str(row["POS"]),row["SURVIVOR_ID"], row["REF"], row["ALT"],str(int(row["Probs"]*1000)),f'{"PASS" if row["Probs"]>filter else "LowProb"}',f'IDS={IDS[:-1]};SUPP={SUPP};SUPP_VEC={SUPP_VEC};SVLEN={row["SVLEN"]};SVTYPE={svtype};END={row["END"]}']
        vcf_file.write("\t".join(line)+"\n")

if __name__ == "__main__":
    args = parse_arguments()

    DATA, DATA_INS, mycontigs = main(args)
    
    outfile = create_header(args.survivor, args.output, mycontigs)
    build_vcf(outfile, DATA, "DEL", ['SURVIVOR_ID', "CHROM", "POS", "REF", "ALT", "END", 'BD_ID', 'DELLY_ID', "INSV_ID",'LUMPY_ID', 'MANTA_ID', 'PINDEL_ID', 'TARDIS_ID','Probs'],0.5)
    build_vcf(outfile, DATA_INS, "INS", ['SURVIVOR_ID', "CHROM", "POS", "REF", "ALT", "END", 'BD_ID', 'DELLY_ID', "INSV_ID",'LUMPY_ID', 'MANTA_ID', 'PINDEL_ID', 'TARDIS_ID','Probs'],0.5)
    outfile.close()