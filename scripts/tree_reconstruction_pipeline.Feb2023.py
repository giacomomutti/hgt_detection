#!/usr/bin/env python

"""
    Tree reconstruction pipeline - automated phylogenetic reconstruction pipeline - it resembles the
    steps followed by a phylogenetist to build a gene family tree with error-control
    of every step

    Copyright (C) 2019 - Marina Marcet-Houben, Toni Gabaldon

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
"""

import argparse
import os
import sys
import subprocess as sp
import pandas as pd
import glob
import time

blast_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/blastp"
muscle_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/muscle"
mafft_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/mafft"
kalign_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/kalign"
tcoffee_path = "t_coffee"
trimal_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/trimal"
iqtree_path = "/gpfs/projects/bsc40/project/pipelines/tree_pipeline/bin/iqtree"

###############################################################################
# Basic modules
###############################################################################

def load_sequences(contigFile,delimiter):
    """ This function is needed to load a set of fasta sequences into memory. 
    Params
    --------
    
    contigFile --> contains the name of the fasta file you want to load into memory
    delimiter --> Leave empty in case you don't want to split the header, set the separator you want to use if you want to split the fasta header
        
    Return
    -------
        
    This function will return a dictionary which contained the headers as keys and the sequences as values.
        
    """
    seqs = {}
    name = ""
    s = []
    for line in open(contigFile):
        line = line.strip()
        if ">" in line:
            if name != "":
                seqs[name] = "".join(s)
            if delimiter == "":
                name = line.replace(">","")
            else:
                name = line.replace(">","").split(delimiter)[0]
            s = []
        else:
            s.append(line.upper())
    if len(s) != 0:
        seqs[name] = "".join(s)
    return seqs

def print_sequence(code,sequence,outfile):
    """ Prints sequence into a file
    Param
    -------
        
    code --> Header the sequence should be assigned
    sequence --> The aminoacid or nucleotide sequence to be printed
    outfile --> the handle to the file where the sequence should be printed
        
    Return 
    -------
    No return is expected from this function
        
    """
    print(">"+code,file=outfile)
    i = 0
    while i < len(sequence):
        print(sequence[i:i+60],file=outfile)
        i += 60

#This function will execute a command in bash
def run_command(cmd,ommit):
    """ Runs a command in bash, the ommit option allows the script to
    continue even if the command fails. Not recommended to set it to
    True"""
        
    if ommit:
        try: process = sp.Popen(cmd,shell=True)
        except: pass
        process.communicate("Y\n")
        if process.wait() != 0: print("Error ocurred, but you chose to ommit it")
    else:
        try: process = sp.Popen(cmd,shell=True)
        except OSError: sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0: sys.exit("ERROR: Execution cmd failed")


def create_folder(name):
    """ Creates a folder if it's not present """
    if not os.path.exists(name):
        cmd = "mkdir "+name
        try:
            run_command(cmd,False)
        except:
            print("Unable to create directory ",name)

def remove_file(fileName):
    """ Deletes files """
    if "*" in fileName:
        files = [f for f in glob.glob(fileName)]
        if len(files) != 0:
            cmd = "rm -rf "+fileName
            run_command(cmd,True)
    else:
        if os.path.exists(fileName):
            cmd = "rm -rf "+fileName
            run_command(cmd,False)

###############################################################################
# Blast related modules
###############################################################################

def run_blast(inFile,dbFile,outH,outLog):
    """ Performs a homology search	"""
    cmd = blast_path+" -query "+inFile+" -db "+dbFile+" -outfmt 6 -out "+ \
    outH+" -dbsize 100000 -evalue 1 -max_target_seqs 5000"
    print(cmd,file=outLog)
    run_command(cmd,False)

def filter_blast_results(infile, outfile, outfileSeqs, dbFile, evalue, overlap, numHits,seqsLen,outLog):
    """ Filters blast results based on float,evalue and numHits """
    if os.stat(infile).st_size == 0:
        print("No blast results were found",file=outLog)
        exit()
    df = pd.read_csv(infile,sep="\t",header=None)
    seqs = load_sequences(dbFile," ")
    len_seq = len(seqs[df[0][0]])
    df[12] = (df[7] - df[6] + 1) / len_seq
    df_filtered = df[(df[10] < evalue) & (df[12] > overlap)]
    df_filtered.to_csv(outfile,sep="\t",header=False,index=False)
    num_seqs = 0
    with open(outfileSeqs,"w") as outF:
        for index,row in df_filtered.iterrows():
            sequence = str(seqs[row[1]]).replace("*","")
            if seqsLen:
                if len(sequence) > seqsLen:
                    pass
                else:
                    if num_seqs <= numHits:
                        outF.write(">"+row[1]+"\n"+sequence+"\n")
                        num_seqs += 1
            else:
                if num_seqs <= numHits:
                    outF.write(">"+row[1]+"\n"+sequence+"\n")
                    num_seqs += 1

###############################################################################
# Functions related to alignments
###############################################################################

def clean_seqs(outSeqs,outSeqsR):
    """ Removes non-regular aminoacids from sequence and replaces them with X.
    It also returns the file with the reverse sequence"""
    seqs = load_sequences(outSeqs,"")
    with open(outSeqs,"w") as outfile:
        with open(outSeqsR,"w") as outfileR:
            for code in seqs:
                s = seqs[code].upper()
                s = s.replace("J","X").replace("B","X").replace("O","X").replace("U","X").replace("Z","X")
                print_sequence(code,s,outfile)
                print_sequence(code,s[::-1],outfileR)

def get_cds_seqs(outSeqs,cdsFile,outCDS):
    """ If a CDS file was provided, the program will return a .seqs_cds file"""
    seqs = load_sequences(outSeqs,"")
    seqsCDS = load_sequences(cdsFile,"")
    with open(outCDS,"w") as outfile:
        for code in seqs:
            if code in seqsCDS:
                print_sequence(code,seqsCDS[code],outfile)
            else:
                sys.exit("Sequence "+code+" was not found in the CDS file")

def run_muscle(outSeqs,outSeqsR,outAlg,completed_algs,outLog):
    """ Runs the muscle alignments """
    if not os.path.exists(outAlg+".F.muscle") or args.replace:
        try:
            cmd = muscle_path+" -in "+outSeqs+" -out "+outAlg+".F.muscle"
            print(cmd,file=outLog)
            run_command(cmd,False)
            completed_algs.append(outAlg+".F.muscle")
            print("Muscle forward alignment completed",file=outLog)
        except:
            print("Muscle forward alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".F.muscle")
    if not os.path.exists(outAlg+".RF.muscle") or args.replace:
        try:
            cmd = muscle_path+" -in "+outSeqsR+" -out "+outAlg+".R.muscle"
            print(cmd,file=outLog)
            run_command(cmd,False)
            turn_alignment(outAlg+".R.muscle",outAlg+".RF.muscle")
            completed_algs.append(outAlg+".RF.muscle")
            print("Muscle Reverse alignment completed",file=outLog)
        except:
            print("Muscle Reverse alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".RF.muscle")
    return completed_algs

def turn_alignment(inFile,outFile):
    """ Given a fasta files it turns all the sequences arround """
    seqs = load_sequences(inFile,"")
    with open(outFile,"w") as outfile:
        for code in seqs:
            print_sequence(code,seqs[code][::-1],outfile)
        
def run_mafft(outSeqs,outSeqsR,outAlg,completed_algs,outLog):
    """ Runs the mafft alignments """
    if not os.path.exists(outAlg+".F.mafft") or args.replace:
        try:
            cmd = mafft_path+" --auto "+outSeqs+" > "+outAlg+".F.mafft"
            print(cmd,file=outLog)
            run_command(cmd,False)
            completed_algs.append(outAlg+".F.mafft")
            print("Mafft forward alignment completed",file=outLog)
        except:
            print("Mafft forward alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".F.mafft")
    if not os.path.exists(outAlg+".RF.mafft") or args.replace:
        try:
            cmd = mafft_path+" --auto "+outSeqsR+" > "+outAlg+".R.mafft"
            print(cmd,file=outLog)
            run_command(cmd,False)
            turn_alignment(outAlg+".R.mafft",outAlg+".RF.mafft")
            completed_algs.append(outAlg+".RF.mafft")
            print("Mafft Reverse alignment completed",file=outLog)
        except:
            print("Mafft Reverse alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".RF.mafft")
    return completed_algs
            
def run_kalign(outSeqs,outSeqsR,outAlg,completed_algs,outLog):
    """ Runs the kalign alignments """
    if not os.path.exists(outAlg+".F.kalign") or args.replace:
        try:
            cmd = kalign_path+" -f fasta -i "+outSeqs+" -o "+outAlg+".F.kalign"
            print(cmd,file=outLog)
            run_command(cmd,False)
            completed_algs.append(outAlg+".F.kalign")
            print("Kalign Forward alignment completed",file=outLog)
        except:
            print("Kalign Forward alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".F.kalign")
    if not os.path.exists(outAlg+".RF.kalign") or args.replace:
        try:
            cmd = kalign_path+" -f fasta -i "+outSeqsR+" -o "+outAlg+".R.kalign"
            print(cmd,file=outLog)
            run_command(cmd,False)
            turn_alignment(outAlg+".R.kalign",outAlg+".RF.kalign")
            completed_algs.append(outAlg+".RF.kalign")
            print("Kalign Reverse alignment completed",file=outLog)
        except:
            print("Kalign Reverse alignment not completed",file=outLog)
    else:
        completed_algs.append(outAlg+".RF.kalign")
    return completed_algs
        
def run_mcoffee(outSeqs,completed_algs,outAlg,outLog):
    """ Performs a consensus alignment with M-coffee """
    cmd = tcoffee_path+" "+outSeqs+" -n_core 1 -output fasta -quiet -aln "+ \
          ",".join(completed_algs)+" -outfile "+outAlg+".metalig"
    print(cmd,file=outLog)
    run_command(cmd,False)
    print("M-coffee completed",file=outLog)

def run_trimal(outAlg,completed_algs,outPaths,outLog):
    """ This function will turn the metalign into phylip and trim it"""
    cmd = trimal_path+" -in "+outAlg+".metalig -out "+outAlg+".metalig -phylip"
    print(cmd,file=outLog)
    run_command(cmd,False)
    with open(outPaths,"w") as outfile:
        print("\n".join(completed_algs),file=outfile)
    ct_thr = 1 / len(completed_algs)
    cmd = trimal_path+" -compareset "+outPaths+" -forceselect "+outAlg+ \
        ".metalig -out "+outAlg+".interm.clean -phylip -ct "+str(ct_thr)+" -cons 30"
    print(cmd,file=outLog)
    try:
        run_command(cmd,False)
    except:
        print("Trimal failed",file=outLog)
        exit()
    cmd = trimal_path+" -in "+outAlg+".interm.clean -gt 0.1 -cons 30 -out " +\
        outAlg+".clean"
    print(cmd,file=outLog)
    try:
        run_command(cmd,False)
    except:
        print("Trimal failed",file=outLog)
        exit()
    print("Trimming completed",file=outLog)	

def run_trimal_cds(outAlg,completed_algs,outPaths,cdsFile,outLog):
    """ Performs the backtranslation when a CDS file is provided """
    cmd = trimal_path+" -backtrans "+cdsFile+ \
        " -out "+outAlg+".interm.clean_cds -compareset "+outPaths+ " -forceselect "+outAlg+".metalig -ct 0.1667 -splitbystopcodon"
    print(cmd,file=outLog)
    try:
        run_command(cmd,False)
    except:
        print("Trimal failed",file=outLog)
        exit()
    cmd = trimal_path+" -in "+outAlg+".interm.clean_cds -out "+outAlg+".clean_cds -phylip\
        -gt 0.1 -cons 30"
    print(cmd,file=outLog)
    try:
        run_command(cmd,False)
    except:
        print("Trimal failed", file=outLog)
        exit()
    print("Conversion to CDS alignment completed",file=outLog)
        
###############################################################################
# Tree reconstruction function
###############################################################################

def reconstruct_tree(algFile,treeFile,mlFile,modelList,numseqs,treeThreads,outLog):
    """ Reconstructs phylogenetic tree using iqtree and then processes results"""
    cmd = iqtree_path+" -nt "+treeThreads+" -quiet -mem 4G" \
        " -cmin 4 -cmax 10 -s "+algFile
    if numseqs > 4:
        cmd += " -bb 1000"
    if len(modelList) != 0:
        cmd += " -mset "+modelList
    print(cmd,file=outLog)
    run_command(cmd,False)
    print("Tree reconstructed",file=outLog)
    process_tree_files(algFile,treeFile,mlFile)

def reconstruct_tree_cds(algFile,treeFile,mlFile,numseqs,treeThreads,outLog):
    """ Reconstructs phylogenetic tree using a codon model """
    cmd = iqtree_path+" -nt "+treeThreads+" -st CODON -quiet -mem 4G -cmin 4 -cmax 10 -s "+algFile+" -redo"
    if numseqs > 4:
        cmd += " -bb 1000"
    print(cmd,file=outLog)
    run_command(cmd,False)
    process_tree_files(algFile,treeFile,mlFile)
    print("Tree reconstructed",file=outLog)

def process_tree_files(algFile,treeFile,mlFile):
    """ Processes the files created by iqtree and only saves relevant info"""
    cmd = "mv "+algFile+".treefile "+treeFile
    run_command(cmd,False)
    model = [line.strip().split(": ")[1] for line in open(algFile+".iqtree") if "Model of substitution:" in line]
    lk = [line.strip().split() for line in open(algFile+".iqtree") if model[0] in line and ":" not in line]
    with open(mlFile,"w") as outfile:
        print(model[0]+"\t"+lk[0][1],file=outfile)
    if not args.keep:
        cmd = "rm -rf "+algFile+".*"
        run_command(cmd,False)

desc = ("Pipeline used to reconstruct trees for phylomeDB\n")
parser = argparse.ArgumentParser(description = desc, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-i", "--in", dest = "inFile", type = str, default = None, help = "Input file containing the query sequence/s or alignment")
parser.add_argument("-d", "--db", dest = "dbFile", type = str, default = None, help = "Input file containing the target sequence database")
parser.add_argument("--cds", dest = "cdsFile", type = str, default = None, help = "Input file containing CDS corresponding to input protein seqs")
parser.add_argument("-o", "--out", dest = "outFolder", type = str, default = ".", help = "Output folder where all generated files will be dumped")
parser.add_argument("-s", "--steps", dest = "steps", type = str, choices=["hat","h","a","t","ha","at"], default = "hat", help = "Indicates which parts of the pipeline will be executed. h = homology, a = alignment, t = tree reconstruction")
parser.add_argument("-r", "--replace", dest = "replace", default = False, action = "store_true", help = "Over-write any previously generated file")
parser.add_argument("-R", "--hard_replace", dest = "hard_replace", default = False, action = "store_true", help = "First delete all the files previously produced and then start analysis")
parser.add_argument("-k","--keep",dest = "keep", default = False, action = "store_true", help = "Keeps intermediate files of the analysis. Default is to delete them.")
parser.add_argument("--overlap",dest="overlap", type = float, default = 0.33, help = "Homology step: Overlap by which the blast hits will be filtered")
parser.add_argument("--evalue",dest="evalue", type = float, default = 1e-05, help = "Homology step: E-value by which the blast hits will be filtered")
parser.add_argument("--numHits",dest="numHits", type = int, default = 150, help = "Homology step: Maximum number of hits that will be selected")
parser.add_argument("--minNumAlgs",dest="minNumAlgs", type = int, default = 2, help = "Minimum number of alignments to be completed to build tree")
parser.add_argument("--models",dest="models", default = "DCmut,JTTDCMut,LG,WAG,VT", help = "List of models for iqtree, separated by a ,")
parser.add_argument("--treeThreads",dest="treeThreads", default = "1", help = "Number of threads to run iqtree")
parser.add_argument("--seqLen",dest="seqLen", default = None, help = "Maximum sequence length allowed")
args = parser.parse_args()

#Check the presence of the input files and create the output folder if it's not
# present
inFile = args.inFile

if not os.path.exists(inFile):
    sys.exit("Provided input file does not exist.")

if "h" in args.steps:
    if not os.path.exists(args.dbFile):
        sys.exit("DB file does not exist.")

    if not os.path.exists(args.dbFile+".phr"):
        sys.exit("DB has not be formatted. Please use formatdb to format it")

if args.cdsFile:
    if not os.path.exists(args.cdsFile):
        sys.exit("You provided a CDS file, but it cannot be found.")

if not os.path.exists(args.outFolder):
    create_folder(args.outFolder)
       
if args.hard_replace:
    fileNames = [fileName for fileName in glob.glob(args.outFolder+"/*")]
    for fileName in fileNames:
        if ".fasta" not in fileName:
            cmd = "rm -rf "+fileName
            run_command(cmd,False)

os.chdir(args.outFolder)
tag = args.inFile.split("/")[-1].split(".")[0]
outH = args.outFolder+"/"+tag+".homology.blast.out"
outHF = args.outFolder+"/"+tag+".homology.blast.filter"
outSeqs = args.outFolder+"/"+tag+".seqs"
outSeqsR = args.outFolder+"/"+tag+".reverse.seqs"
outAlg = args.outFolder+"/"+tag+".alg"
outPaths = args.outFolder+"/"+tag+".paths"
treeThreads = args.treeThreads
if args.seqLen:
    seqsLen = int(args.seqLen)
else:
    seqsLen = args.seqLen

#Start log file:
if not os.path.exists(args.outFolder+"/"+tag+".log"):
    outLog = open(args.outFolder+"/"+tag+".log","w")
    print("Starting log:",time.strftime('%X %x'),file=outLog)
else:
    outLog = open(args.outFolder+"/"+tag+".log","a")
    print("Re-starting log:",time.strftime('%X %x'),file=outLog)

if not args.cdsFile:
    outTree = args.outFolder+"/"+tag+".tree.ml.iqtree.nw"
    outML = args.outFolder+"/"+tag+".rank.ml.txt"
    finalAlg = outAlg+".clean"
else:
    outTree = args.outFolder+"/"+tag+".cds.tree.ml.iqtree.nw"
    outCDS_seqs = args.outFolder+"/"+tag+".seqs_cds"
    finalAlg = outAlg+".clean_cds"
    outML = args.outFolder+"/"+tag+".cds.rank.ml.txt"

if "h" not in args.steps:
    outSeqs = args.inFile

if "h" not in args.steps and "a" not in args.steps:
    finalAlg = args.inFile
        
if "h" in args.steps:
    print("Searching for homologs - ",time.strftime('%X %x'),file=outLog)
    if not os.path.exists(outH) or args.replace:
        run_blast(args.inFile, args.dbFile, outH,outLog)
    if not os.path.exists(outHF) or args.replace:
        filter_blast_results(outH,outHF,outSeqs,args.dbFile,args.evalue,args.overlap,args.numHits,seqsLen,outLog)

if "a" in args.steps:
    print("Building alignments - ",time.strftime('%X %x'),file=outLog)
    numseqs = len([line.strip() for line in open(outSeqs) if ">" in line])
    if numseqs >= 3:
        if not os.path.exists(outAlg+".clean") or args.replace:
            if not os.path.exists(outSeqsR) or args.replace:
                clean_seqs(outSeqs,outSeqsR)
            if args.cdsFile:
                if not os.path.exists(outCDS_seqs) or args.replace:
                    get_cds_seqs(outSeqs,args.cdsFile,outCDS_seqs)
            completed_algs = []
            completed_algs = run_muscle(outSeqs,outSeqsR,outAlg,completed_algs,outLog)
            completed_algs = run_mafft(outSeqs,outSeqsR,outAlg,completed_algs,outLog)
            completed_algs = run_kalign(outSeqs,outSeqsR,outAlg,completed_algs,outLog)
            if len(completed_algs) >= args.minNumAlgs:
                if not os.path.exists(outAlg+".metalig") or args.replace:
                    run_mcoffee(outSeqs,completed_algs,outAlg,outLog)
                if not os.path.exists(outAlg+".clean") or args.replace:
                    run_trimal(outAlg,completed_algs,outPaths,outLog)
                if args.cdsFile:
                    if not os.path.exists(outAlg+".clean_cds") or args.replace:
                        run_trimal_cds(outAlg,completed_algs,outPaths,outCDS_seqs,outLog)

if "t" in args.steps:
    print("Building tree - ",time.strftime('%X %x'),file=outLog)
    if "a" not in args.steps:
        try:
            numseqs = int([line.split() for line in open(finalAlg)][0][0])
        except:
            sys.exit("Alignment needs to be in phylip format")
    else:
        seqs = load_sequences(outSeqs,"")
        numseqs = len(set([s for s in seqs.values()]))
        # ~ numseqs = len([line.strip() for line in open(outSeqs) if ">" in line])
    if numseqs >= 3:
        if not args.cdsFile:
            if not os.path.exists(outTree) or args.replace:
                reconstruct_tree(finalAlg,outTree,outML,args.models,numseqs,treeThreads,outLog)
        else:
            if not os.path.exists(outTree) or args.replace:
                reconstruct_tree_cds(finalAlg,outTree,outML,numseqs,treeThreads,outLog)

outLog.close()

if not args.keep:
    remove_file(args.outFolder+"/"+tag+".dnd")
    remove_file(outSeqsR)
    remove_file(outAlg+".F*")
    remove_file(outAlg+".R*")
    remove_file(outPaths)
    remove_file(outAlg+".clean.*")

