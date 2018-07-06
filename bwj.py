# module load python/3.2.2
# python /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/bwj.py -p /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/integrated_call_male_samples_v3.20130502.ALL.panel -o /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/output/chr21 -m /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/chr21.1kg.phase3.v5a.markers -b /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/chr21.1kg.phase3.v5a.bgl -t /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/Prins_28887542_hdl_chr21.txt --pop EUR --bin /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/ImpG-master/ImpG-Bins --maf 0.0001 --lambd 0.1

import sys
import math
import argparse
import os
import bisect
import subprocess

# Specify arguments for the python utility
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-m","--marker",help="Input map file for all markers.",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("-p","--panel",help="Input panel file with sample ID and population of 1000 Genome database.",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("-b","--bgl",help="Input bgl file as a Beagle file. Note: please unzip the file if ever compressed.",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("-t","--typed",help="Z-scores of typed SNPs",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("-o","--output",help="output directory to store all results.",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("--bin",help="Directory to store ImpG-Summary-GenBeta binary code, perhaps need 'make' first. The directory should contain 'ImpG-Summary' and 'ImpG-Summary-GenBeta'.",
                    required=True,default=argparse.SUPPRESS)
parser.add_argument("--pop",help="Specify a population in 1000 Genome database",
                    choices=["EUR","AFR","AMR","EAS","SAS"],default="ALL")
parser.add_argument("--windows",help="window width",default=1000000,type=int)
parser.add_argument("--buffer",help="buffer size",default=500000,type=int)
parser.add_argument("--maf",help="maf threshold passing to ImpG-Summary-GenBeta function",default=0.0001,type=float)
parser.add_argument("--lambd",help="lambda value passing to ImpG-Summary-GenBeta function",default=0.1,type=float)

args=parser.parse_args()


# Construct the directories for the whole process
if not os.path.exists(args.output):
    os.makedirs(args.output)

dir_maps = args.output + "/maps/"
if not os.path.exists(dir_maps):
    os.makedirs(dir_maps)

dir_haps = args.output + "/haps/"
if not os.path.exists(dir_haps):
    os.makedirs(dir_haps)

dir_typed = args.output + "/typed/"
if not os.path.exists(dir_typed):
    os.makedirs(dir_typed)

dir_betas = args.output + "/betas/"
if not os.path.exists(dir_betas):
    os.makedirs(dir_betas)

dir_imp = args.output + "/imp/"
if not os.path.exists(dir_imp):
    os.makedirs(dir_imp)

# Specify some common variables
pop_panel = args.output + "/phase3.1KG." + args.pop + ".txt"
file_map = args.output + "/maps.names"
file_typed_map = args.output + "/typed.maps.names"
oprefix = dir_maps + "phase3.1KG." + args.pop
window = int(args.windows)
buf = int(args.buffer/2)
binary_ImpZ = args.bin + "/ImpG-Summary"
binary_GenBeta = args.bin + "/ImpG-Summary-GenBeta"

# Subfunctions for the whole process
def ref_pop(input_panel,pop_panel):
    f = open(input_panel,mode="r")
    w = open(pop_panel,mode="w")
    for line in f:
        line = line[:-1]
        cols = line.split()
        indv = cols[0]
        if indv == "sample":
            continue
        pop = cols[2]
        if pop == args.pop or args.pop == "ALL":
            w.write(indv+"\n")
    w.close()
def split_map(fname,oprefix,window,buf,sequence_nm):
    f = open(fname, 'r')
    all_snps = []
    positions = []
    for line in f:
        line = line[:-1]
        cols = line.split()
        # ignore lines that have more than 4 fields
        if(len(cols) != 4):
            continue
        # parse everything out and check for format
        snp_id = cols[0]
        prefix = snp_id[:2] # should be "rs" for snps
        suffix = snp_id[2:] # should be a number for snps
        # snp_id should be "rs" followed by a number, otherwise, ignore it
        if(prefix != "rs" or (not suffix.isdigit())):
            continue
        # snp_pos should be a number
        snp_pos = cols[1]
        if(not snp_pos.isdigit()):
            continue
        # should be A, T, C, G
        ref = cols[2]
        if(ref != "A" and ref != "T" and ref != "C" and ref != "G"):
            continue
        alt = cols[3]
        # should be A, T, C, G
        if(alt != "A" and alt != "T" and alt != "C" and alt != "G"):
            continue
        all_snps.append((snp_id, ref+" "+alt))
        positions.append(int(snp_pos))
    f.close()

    min_pos = positions[0]
    max_pos = positions[len(all_snps)-1]

    min_win_id = (min_pos-1)//window
    max_win_id = (max_pos-1)//window
    # iterate through windows
    total_count = 0
    s = open(sequence_nm, "w")
    for i in range(min_win_id, max_win_id+1):
        # compute window range
        low = i*window-buf
        high = (i+1)*window+buf
        if(i == 0):
            low = 1
        if(i == max_win_id):
            high = (i+1)*window
        # index of the snp, with position greater than low
        low_idx =  bisect.bisect(positions, low)
        # index of the snp, with position greater than high
        high_idx = bisect.bisect(positions, high)
        if(high_idx - low_idx > 0):
            name = oprefix+"."+str(low)+"-"+str(high)+".map"
            name_fields = name.split("/")
            s.write(name_fields[len(name_fields)-1]+"\n")
            w = open(name, "w")
            w.write("SNP_ID SNP_Pos Ref_Allele Alt_Allele\n")
            window_count = 0
            start_idx = low_idx-1
            if(start_idx < 0):
                start_idx = 0
            for j in range(start_idx, high_idx):
                line = all_snps[j][0]+" "+str(positions[j])+" "+all_snps[j][1]
                w.write(line+"\n")
                window_count = window_count + 1
            total_count = total_count + window_count
            w.close()
    s.close()
def split_hap(sequencefile_nm,hapfile_nm,popfile_nm,output,directory):
    # load the set of individual sample ids of the population we choose
    indvs = set()
    popfile = open(popfile_nm, "r")
    for line in popfile:
        line = line[:-1]
        indvs.add(line)
    popfile.close()

    # find the columns corresponding to the individual sample id
    indv_line = "SNP_id "
    indv_idx = list()
    hapfile = open(hapfile_nm, "r")
    line = hapfile.readline()
    line = line[:-1]
    cols = line.split()
    for i in range(2, len(cols)):
        if (cols[i] in indvs):
            indv_line = indv_line + cols[i] + " "
            indv_idx.append(i)

    # read in the file that contains all the haplotype
    all_haps = dict()
    for line in hapfile:
        line = line[:-1]
        cols = line.split()
        snp_id = cols[1]
        if (snp_id[0:2] != "rs"):
            continue
        hap_line = ""
        for idx in indv_idx:
            hap_line = hap_line + cols[idx]
        all_haps[snp_id] = hap_line
    hapfile.close()

    print("All " + str(len(all_haps)) + " loaded...")

    # iterate through all the map files
    sequencefile = open(sequencefile_nm, "r")
    for mapfile_nm in sequencefile:
        print(1)
        print(mapfile_nm)
        mapfile_nm = mapfile_nm[:-1]
        print(2)
        print(mapfile_nm)
        print("Extracting haplotypes for windows in " + mapfile_nm)
        mapfile_nm = directory + mapfile_nm

        # iterate through one map file
        of = open(output + mapfile_nm.split("/")[-1] + ".haps", "w")
        of.write(indv_line + "\n")
        mapfile = open(mapfile_nm, "r")
        first_line_skipped = False
        for line in mapfile:
            if (not first_line_skipped):
                first_line_skipped = True
                continue
            line = line[:-1]
            cols = line.split()
            name = cols[0]
            pos = cols[1]
            refalt = cols[2] + cols[3]
            if (name not in all_haps):
                continue
            hap_line = all_haps[name]
            of.write(name + " ")
            for i in range(0, len(hap_line)):
                if (hap_line[i] == refalt[0]):
                    of.write("1")
                else:
                    of.write("2")
            of.write("\n")
        mapfile.close()
        of.close()
def split_type(sequencefile_nm,typedfile_nm,output,directory,order_nm,th):
    # read the snps on the array
    typed_snps = dict()
    typedfile = open(typedfile_nm, "r")
    for line in typedfile:
        line = line[:-1]
        cols = line.split()
        snp_id = cols[0]
        if(snp_id[0:2] == "rs"):
            typed_snps[snp_id] = (cols[1], cols[2], cols[3], cols[4])
    typedfile.close()

    # iterate through all the map files
    first_line = "SNP_id SNP_pos Ref Alt Z-score"
    order = open(order_nm, "w")
    sequencefile = open(sequencefile_nm, "r")
    for mapfile_nm in sequencefile:
        count = 0
        mapfile_nm = mapfile_nm[:-1]
        mapfile_nm = directory+mapfile_nm
        of = open(output+mapfile_nm.split("/")[-1]+".typed", "w")
        of.write(first_line+"\n")
        mapfile = open(mapfile_nm, "r")
        for line in mapfile:
            line = line[:-1]
            cols = line.split()
            snp_id = cols[0]
            if(snp_id in typed_snps):
                snp_info = typed_snps[snp_id]
                if(snp_info[1]==cols[2] and snp_info[2]==cols[3]):
                    # of.write(snp_id+" "+snp_info[0]+" "+snp_info[1]+" "+snp_info[2]+" "+snp_info[3]+"\n")
                    of.write(snp_id+" "+cols[1]+" "+snp_info[1]+" "+snp_info[2]+" "+snp_info[3]+"\n")
                elif(snp_info[1]==cols[3] and snp_info[3]==cols[2]):
                    # of.write(snp_id+" "+snp_info[0]+" "+snp_info[2]+" "+snp_info[1]+" "+str((-1)*float(snp_info[3]))+"\n")
                    of.write(snp_id+" "+cols[1]+" "+snp_info[2]+" "+snp_info[1]+" "+str((-1)*float(snp_info[3]))+"\n")
                else:
                    continue
                count = count + 1
        mapfile.close()
        of.close()
        if(count > th):
            order.write(mapfile_nm.split("/")[-1]+"\n")
    sequencefile.close()
    order.close()
def gen_beta(sequencefile_nm,directory,binary,maf,lambd):


    # iterate through all the map files
    sequencefile = open(sequencefile_nm, "r")
    for mapfile_nm in sequencefile:
        mapfile_nm = mapfile_nm[:-1]
        cmd_line = binary + " -f " + str(maf) + " -l " + str(lambd)
        cmd_line = cmd_line + " -h " + directory + "/haps/" + mapfile_nm + ".haps"
        cmd_line = cmd_line + " -t " + directory + "/typed/" + mapfile_nm + ".typed"
        cmd_line = cmd_line + " -m " + directory + "/maps/" + mapfile_nm
        cmd_line = cmd_line + " -p " + directory + "/betas/" + mapfile_nm + ".step1"
        subprocess.call(cmd_line.split())
    sequencefile.close()
def imp_z(sequencefile_nm,directory,binary):

    # iterate through all the map files
    sequencefile = open(sequencefile_nm, "r")
    for mapfile_nm in sequencefile:
        mapfile_nm = mapfile_nm[:-1]
        cmd_line = binary + " "
        cmd_line = cmd_line + " -t " + directory + "/typed/" + mapfile_nm + ".typed"
        cmd_line = cmd_line + " -m " + directory + "/maps/" + mapfile_nm
        cmd_line = cmd_line + " -p " + directory + "/betas/" + mapfile_nm + ".step1"
        cmd_line = cmd_line + " -o " + directory + "/imp/" + mapfile_nm + ".impz"
        subprocess.call(cmd_line.split())
    sequencefile.close()
def merge_z(directory,window,buf,sequence_nm):

    output = directory + "/all.impz"
    # adjust start position and end position
    seqfile = open(sequence_nm, "r")
    start_end = list()
    for line in seqfile:
        line = line[:-1]
        print(line)
        cols = line.split(".")
        low_high = cols[-2].split("-")
        low = int(low_high[0])
        high = int(low_high[1])
        start_end.append((low + buf + 1, high - buf))
    start_end[0] = (1, start_end[0][1])
    nmap = len(start_end)
    start_end[nmap - 1] = (start_end[nmap - 1][0], start_end[nmap - 1][1] + buf)
    seqfile.close()

    # iterate through all z-score file
    out = open(output, "w")
    out.write("SNP_id SNP_pos Ref_allele Alt_Allele Z-score Var\n")
    idx = 0
    seqfile = open(sequence_nm, "r")
    for line in seqfile:
        line = line[:-1]
        zscfile = open(directory + "/imp/" + line + ".impz", "r")
        first_line_skipped = False
        for zline in zscfile:
            if (first_line_skipped == False):
                first_line_skipped = True
                continue
            zline = zline[:-1]
            cols = zline.split()
            pos = int(cols[1])
            if (pos >= start_end[idx][0] and pos <= start_end[idx][1]):
                out.write(zline + "\n")
        idx = idx + 1
    seqfile.close()
    out.close()


########## The below is the code for the whole process


# Step 1: Output a reference sample file based on specific population
ref_pop(args.panel,pop_panel)

# Step 2: Generate SNP mapping files for partitions of a chromosome
split_map(args.marker,oprefix,window,buf,file_map)

# Step 3: Generate haplotype reference panel files for partitions of a chromosome
split_hap(file_map,args.bgl,pop_panel,dir_haps,dir_maps)

# Step 4: Generate z-scores of typed SNPs for partitions of a chromosome (check the REF and ALT allele at the same time)
th = 10   # a threshold to specify the minimal SNP number in one SNPs partition, not so important
split_type(file_map,args.typed,dir_typed,dir_maps,file_typed_map,th)

# Step 5: Generate beta files for all the partitions of one chromosome
gen_beta(file_map,args.output,binary_GenBeta,args.maf,args.lambd)

# Step 6: Calculate imputation z-scores
imp_z(file_map,args.output,binary_ImpZ)

# Step 7: Merge the Z-scores files computed by ImpG-Summary into one file
merge_z(args.output,window,buf,file_map)
