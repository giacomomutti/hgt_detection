#!/usr/bin/env python3
from glob import glob
import pandas as pd
import string
import argparse
import os
import subprocess as sp

def parse_args():
    parser=argparse.ArgumentParser(description="Get a new taxdump including new close proteomes")
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
                    help='Input file: mnemo\ttaxonomy'),
    parser.add_argument('-d', '--data_dir', dest='taxdump', required=True, type=str,
                help='taxdump directory'),
    parser.add_argument('-t', '--taxids', dest='taxids', required=True, type=str,
            help='taxonomy to lineage already present'),
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, type=str,
            help='output new taxdump directory'),
    parser.add_argument('-s', '--seqdir', dest='seqs', required=True, type=str,
        help='directory where raw sequences are stored. Name of file must be the same as mnemo')
    args=parser.parse_args()
    return args

# taxdump_edit.pl -names names.dmp -nodes nodes.dmp -taxa NAME -parent XXX -rank NAME -division X

alph = string.ascii_uppercase

def clean_seq(seq):
    oseq = seq.replace('J', 'X').replace('B', 'X').replace('O', 'X')
    oseq = oseq.replace('U', 'X').replace('Z', 'X').replace('*', '')

    return oseq


def read_fasta(fafile):
    seqs = {}
    seqid = ''

    for line in open(fafile):
        line = line.replace('\n', '')
        if '>' in line:
            if seqid == '':
                seqid = line.replace('>', '').split(' ', 1)[0]
                s = []
            else:
                seqs[seqid] = clean_seq(''.join(s))
                seqid = line.replace('>', '')
                s = []
        else:
            s.append(line)

    if len(s) > 0:
        seqs[seqid] = clean_seq(''.join(s))

    return(seqs)


def write_fasta(seqs, rename=False, taxid=None, taxiddic = None, virus=False, filename=None,
                seqlen=60, ofilenm=None, taxidmap=None):
    ostr = ''
    if rename:
        i = 0
        j = 0
        k = 0

    for seq in seqs:
        if rename:
            if taxid is not None:
                protid = alph[i] + alph[j] + str(k).zfill(len(str(len(seqs))))
                ostr += ('>%s_%s_%s\n' % (taxid, filename, protid))
            elif virus and taxiddic is not None:
                protid = alph[i] + alph[j] + str(k).zfill(len(str(len(seqs))))
                virus = seq.split('|')[4]
                ostr += ('>%s_%s_%s\n' % (taxiddic[virus], filename, protid))


            j += 1
            if i > len(alph) - 1:
                i = 0
                if j > len(alph) - 1:
                    j = 0
            elif j > len(alph) - 1:
                j = 0
                i += 1
                if i > len(alph) - 1:
                    i = 0
            k += 1

        else:
            ostr += ('>%s\n' % seq)
        for l in range(0, len(seqs[seq]), seqlen):
            ostr += seqs[seq][l:l+seqlen] + '\n'

        if taxidmap is not None:
            ofile = open(taxidmap, 'a')
            ofile.write(('%s_%s_%s\t%s\n' % (taxid, filename, protid, taxid)))
            ofile.close()
    
    if ofilenm is not None:
        ofile = open(ofilenm, 'w')
        ofile.write(ostr)
        ofile.close()

    return(ostr)

keydict = {'p':'phylum','c':'class','o':'order', 
           'f': 'family', 'g': 'genus', 's': 'species'}


if __name__ == '__main__':
    inputs=parse_args()
    dir_path = inputs.outdir+'/taxdump/'
    
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    
    eukprot_tax = dir_path+'/new_eukprot_tax.tsv'
    new_ids_txt = dir_path+'/new_ids.tsv'

    taxids = pd.read_csv(inputs.taxids, sep="\t", names=['taxid', 'k', 'p', 'c', 'o', 'f', 'g', 's'])

    sp.run('cp %s/nodes.dmp %s/names.dmp %s' % (inputs.taxdump, inputs.taxdump, dir_path), shell=True)
    sp.run('Rscript scripts/get_eukprot_tax.R -i %s -o %s' % (inputs.input, eukprot_tax), shell=True)

    cmdl = list()
    dict_species = {}

    with open(eukprot_tax) as toadd:
        for line in toadd:
            mnemo = line.strip().split('\t')[0]
            line = line.strip().split('\t', 1)[1].split(';')

            for i, content in enumerate(line):
                rank = content.split('__', 1)[0]
                if rank in keydict.keys():
                    name = content.split('__', 1)[1]

                    # if rank == 's':
                    #     if name in taxids['s'].values:
                    #         dict_species[mnemo] = taxids.loc[taxids['s']==name].taxid.item()

                    if name in taxids[rank].values:
                        # print('%s already found' % (name))
                        continue

                    rank = keydict[rank]
                    parent = line[i - 1].split('__', 1)[1]
                    cmd = ('perl scripts/taxdump_edit.pl -names %s/names.dmp -nodes %s/nodes.dmp '
                        '-taxa \"%s\" -parent $(echo %s | taxonkit name2taxid --data-dir %s | cut -f 2) '
                        '-rank %s -division 8' % (dir_path, dir_path, name, parent, dir_path, rank))
                    cmdl.append(cmd)

    cmdl = list(set(cmdl))

    for rank in keydict.values():
        for item in cmdl:
            if rank in item:
                sp.run(item, shell=True)

    cmd = ('sed \'s/\\t.*s__/\\t/g\' %s | taxonkit name2taxid --data-dir %s -i 2 > %s' % (eukprot_tax, dir_path, new_ids_txt))
    sp.run(cmd, shell=True)

    df = pd.read_csv(new_ids_txt, sep="\t", names=['acc', 'sp', 'taxid'])
    df['genomeid'] = [x.split('.')[0].replace('_', '') for x in df['acc']]
    taxids = dict(zip(df['genomeid'], df['taxid']))

    # taxid_map_file=inputs.taxdump+"taxid.map"
    # taxid_map = pd.read_csv(taxid_map_file, names=['mnemo', 'taxid'], sep='\t')
    # new_taxid_map = pd.DataFrame(taxids.items(), columns=['mnemo', 'taxid'])
    # new_taxid_map['taxid']= new_taxid_map['taxid'].map(str)
    # pd.concat([taxid_map, new_taxid_map]).to_csv(out_taxid_map_file, header=False, sep='\t', index=False)

    files = glob(inputs.seqs+'*')
    odir = inputs.outdir+'/parsed/'

    if not os.path.exists(odir):
        os.makedirs(odir)

    out_taxid_map_file=dir_path+"taxid.map"

    for file in files:
        print('Processing: %s' % file, flush=True)
        filenm = ''.join(file.rsplit('/', 1)[1].split('_')[0:2]).split('.')[0]
        # print(filenm)
        ofilenm = '%s/%s_protein_parsed.fa' % (odir, filenm)
        # print(ofilenm)
        seqs = read_fasta(file)
        for key in seqs.keys():
            header = key.split(' ')[0]
            # new_map.append(('%s\t%s') % (header, taxids[filenm]))

        write_fasta(seqs, True, taxid=taxids[filenm],
                    filename=filenm, ofilenm=ofilenm, taxidmap=out_taxid_map_file)

    sp.run('cat %s/* > %s/all_seqs.fasta' % (odir, odir), shell=True)

    bdb_odir = inputs.outdir+'/blastdb/'
    if not os.path.exists(bdb_odir):
        os.makedirs(bdb_odir)

    sp.run('makeblastdb -parse_seqids -dbtype prot -in %s/all_seqs.fasta -taxid_map %s -out %s/bdb' % (odir, out_taxid_map_file, bdb_odir), shell=True)
