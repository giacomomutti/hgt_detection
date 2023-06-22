import json

configfile: "data/configs/HGT_config.yml"


proteome=config["proteome"]
blast=config["blast"]
db=config["db"]
foldseek_db=config["foldseek_db"]
# self_db=config["self_db"]
kingdom=config["kingdom"]
taxdump=config["taxdump"]

# lineage=config["lineage"] # this is if one day we want to do it at a finer scale

# SAMPLES WITH TAXID 
samples = json.load(open(config["samples"]))
species = list(samples.keys())

# HGTECTOR FILTERS
max_hits=config["max_hits"]

# BLAST FILTERING PARAMS
evalue=config["evalue"]
len_ratio=config["len_ratio"]
min_cov=config["min_cov"]

# HGT FILTERING PARAMS
min_AI=config["min_AI"]
min_AI2=config["min_AI2"]
min_outpct=config["min_outpct"]

# TREE RECONSTRUCTION PIPELINE PARAMS
n_seqs=config["n_seqs"]

# output directories
output=config["output"]


# SAMPLES WITH SEED IDs
samples = json.load(open(config["samples"]))
species = list(samples.keys())


outfiles = []

for sp, seed_file in samples.items():

    id_file = output+f"{sp}/ids/{sp}_all.txt"
    outfiles.append(id_file)

    # this should be optional as some proteomes are not in uniprot and do not have pdbs
    # outfiles.append(output+f"{sp}/pdbs")
    outfiles.append(output+f"{sp}/hgt_trees/{sp}_trees.nwk")


def get_trees(wildcards):
    checkpoint_output = checkpoints.merge_IDs.get(**wildcards).output.all_ids
    with open(checkpoint_output) as all_genes:
        seeds = [gn.strip() for gn in all_genes]
    return expand(output+"{sp}/rooted_tree/rooted_{i}.nwk", sp=wildcards.sp, i=seeds)


def get_tree_plots(wildcards):
    checkpoint_output = checkpoints.merge_IDs.get(**wildcards).output.all_ids
    with open(checkpoint_output) as all_genes:
        seeds = [gn.strip() for gn in all_genes]
    return expand(output+"{sp}/tree_plots/{i}.png", sp=wildcards.sp, i=seeds)


rule hgt_hits:
    input:
        outfiles


# Selfblast to compute self bitscore and eventually use it to define gene families
# with diamond the e-values were different from blast!
rule selfblast:
    input:
        proteome+"{sp}.fa"
    output:
        output+"{sp}/selfblast/{sp}_self.blast",
    threads: 4
    shell:
        """
makeblastdb -parse_seqids -dbtype prot -in {input}
blastp -query {input} -db {input} -evalue 1e-10 \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads} > {output}
"""

# get taxonomy to annotate blast results
rule get_taxonomy:
    input:
        blast+"{sp}.blast"
    output:
        output+"{sp}/taxonomy/{sp}_tax.txt"
    params:
        taxdump=taxdump+"/{sp}/taxdump"
    # conda:
    #     "data"
    shell:
        """
cut -f2 {input} | cut -f1 -d\'_\' | sort | uniq | \
taxonkit reformat --data-dir {params.taxdump} -I1  \
-f \'{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}\' > {output}
"""

# Parse blast to compute new indexes and filter it
rule parse_blast:
    input:
        prot=proteome+"{sp}.fa",
        blast=blast+"{sp}.blast",
        tax=output+"{sp}/taxonomy/{sp}_tax.txt",
        selfblast=output+"{sp}/selfblast/{sp}_self.blast"
    output:
        output+"{sp}/blast/{sp}_filtered.blast"
    params:
        seed_id=lambda wcs: samples[wcs.sp][0]
    shell:
        """
Rscript ./scripts/parse_blast.R -i {input.blast} -t {input.tax}\
 -e {evalue} -c {min_cov} -l {len_ratio}\
 -s {params.seed_id} -a {input.selfblast} -o {output}
"""


# Compute alien indexes
rule compute_HGT_index:
    input:
        output+"{sp}/blast/{sp}_filtered.blast"
    output:
        output+"{sp}/blast/{sp}_filtered_AI.tsv"
    shell:
        "Rscript ./scripts/compute_HGT_index.R -i {input} -k {kingdom} -o {output}"


# Run hgtector
rule HGTector2:
    input:
        prot=proteome+"{sp}.fa",
        blast=output+"{sp}/blast/{sp}_filtered.blast"
    output:
        hgt_blast=output+"{sp}/hgtector/{sp}.blast",
        search=directory(output+"{sp}/hgtector/{sp}_search"),
        analysis=directory(output+"{sp}/hgtector/{sp}_analysis")
    params:
        log=output+"{sp}/hgtector/{sp}_hgtector.log",
        seed_id=lambda wcs: samples[wcs.sp][0],
        tdump=taxdump+"/{sp}/taxdump"
    # conda:
    #     "HGT"
    shell:
        """
awk 'NR>1' {input.blast} | cut -f1,2,3,11,12,13,17 > {output.hgt_blast}
hgtector search -i {input.prot} -m precomp -s {output.hgt_blast} \
-o {output.search} -t {params.tdump} --evalue {evalue} | tee {params.log}

hgtector analyze -i {output.search} -o {output.analysis} \
-t {params.tdump} -k {max_hits} --input-tax {params.seed_id} --outliers none | tee -a {params.log}
"""


# Merge ids from hgtector and AI
checkpoint merge_IDs:
    input:
        hgtector=output+"{sp}/hgtector/{sp}_analysis",
        AI=output+"{sp}/blast/{sp}_filtered_AI.tsv",
    output:
        all_ids=output+"{sp}/ids/{sp}_all.txt",
        both=output+"{sp}/ids/{sp}_both.txt",
        hgtector=output+"{sp}/ids/{sp}_hgtector.txt",
        indexes=output+"{sp}/ids/{sp}_AIs.txt",
        a=temp(output+"{sp}/ids/a_hgt.txt"),
        b=temp(output+"{sp}/ids/b_hgt.txt"),
        plot=directory(output+"{sp}/plots")
    params:
        blast=output+"{sp}/blast/{sp}_filtered.blast"
    shell:
        """
mkdir -p {output.plot}

cut -f1 {input.hgtector}/hgts/* | sort > {output.a}

awk 'NR>1' {input.AI} | awk -v min_AI={min_AI} -v min_AI2={min_AI2} -v min_outpct={min_outpct} \
'$12>min_AI || $17>min_AI2 || $5>min_outpct' | cut -f1 | sort > {output.b}

comm -12 {output.a} {output.b} > {output.both}
comm -23 {output.a} {output.b} > {output.hgtector}
comm -13 {output.a} {output.b} > {output.indexes}

echo "There are $(wc -l < {output.both}) hits in common, $(wc -l < {output.hgtector}) only found by hgtector\
 and $(wc -l < {output.indexes}) only found with the AIs"

cat {output.both} {output.hgtector} {output.indexes} > {output.all_ids}

Rscript scripts/plot_results_HGT.R -i {input.hgtector}/scores.tsv -a {input.AI} -b {params.blast} \
--both {output.both} --onlyai {output.indexes} --onlyhgtector {output.hgtector} -o {output.plot}
"""

rule foldseek:
    input: 
        output+"{sp}/ids/{sp}_all.txt"
    output: 
        pdbs=directory(output+"{sp}/pdbs"),
        foldseek_tmp=temp(output+"{sp}/pdbs/{sp}_foldseek_out.tmp"),
        foldseek_out=output+"{sp}/pdbs/{sp}_foldseek_out.tsv",
        foldseek_html=output+"{sp}/pdbs/{sp}_foldseek_out.html"
    # params:
    #     tmpdir=temp(output+"{sp}/pdbs/tmp")
    threads:
        4
    shell:
        """
python scripts/get_pdb.py -i {input} -o {output.pdbs}

foldseek easy-search {output.pdbs} {foldseek_db} {output.foldseek_tmp} {resources.tmpdir} --threads {threads} \
--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,taxid,lddt,alntmscore,rmsd,prob

cat {output.foldseek_tmp} | taxonkit reformat -I 15 -f \'{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}\' > {output.foldseek_out}

foldseek easy-search {output.pdbs} {foldseek_db} {output.foldseek_html} {resources.tmpdir} --threads {threads}\
 --format-mode 3
"""

rule get_inparalogs:
    input:
        proteome+"{sp}.fa"
    output:
        po_dir=temp(directory(output+"{sp}/queries/{sp}_po")),
        inparalogs=output+"{sp}/queries/{sp}_inparalogs.tsv"
    threads:
        2
    shell:
        """
mkdir -p {output.po_dir}
input=$(realpath {input})
inpar=$(basename {output.inparalogs})
home_dir=$PWD

cd {output.po_dir}

proteinortho -selfblast -clean -cpus={threads} -project={wildcards.sp} $input
awk 'NR>1' {wildcards.sp}.proteinortho.tsv | cut -f4 > ../$inpar

cd $home_dir
"""


# get fasta of filtered blast hits
rule get_fastas:
    input:
        all_ids=output+"{sp}/ids/{sp}_all.txt",
        prot=proteome+"{sp}.fa",
        paralogs=output+"{sp}/queries/{sp}_inparalogs.tsv",
        # dummy=output+"{sp}/dummy/{seed}.txt"
    output:
        ids=output+"{sp}/seed_ids/{seed}.txt",
        fa=output+"{sp}/sequence/{seed}.fa",
        a=output+"{sp}/seed_ids/{seed}_self.txt",
        b=output+"{sp}/seed_ids/{seed}_close.txt",
        c=output+"{sp}/seed_ids/{seed}_others.txt",
        # tmp=temp(output+"{sp}/seed_ids/{sp}_tmpids.txt")
    params:
        # sdb=self_db+"{sp}",
        blast=output+"{sp}/blast/{sp}_filtered.blast",
        sub_db=lambda wcs: samples[wcs.sp][1],
        taxid=lambda wcs: samples[wcs.sp][0]
    shell:
        """
set +o pipefail;
mkdir -p $(dirname {output.ids})
mkdir -p $(dirname {output.fa})

grep "{wildcards.seed}" {input.paralogs} | tr ',' '\\n' | sort | uniq > {output.ids}
echo "{wildcards.seed}" >> {output.ids}

sort -u {output.ids} -o {output.a}
seqkit grep -f {output.a} {input.prot} > {output.fa}

while read gene; do
    awk -v gene="$gene" 'NR>1 && $1==gene && ( $17!="{params.taxid}" && $3!=100 )' {params.blast} \
| awk -v seqs={n_seqs} 'NR<=seqs' | cut -f2 | sort | uniq >> {output.ids}
done < {output.a}

sort -u {output.ids} -o {output.ids}

seqkit grep -f {output.ids} {params.sub_db} >> {output.fa}
seqkit seq -n -i {output.fa} | grep -v -f {output.a} | cat > {output.b}

grep -v -f {output.b} {output.ids} | grep -v -f {output.a} | sort -u > {output.c}
blastdbcmd -db {db} -entry_batch {output.c} >> {output.fa}
"""
# awk '$1=="{wildcards.seed}" && $11<1e-10' {input.selfblast} | cut -f2 \
# blastdbcmd -db {params.sdb} -entry_batch {output.a} > {output.fa}
# blastdbcmd -db {params.sdb} -entry {wildcards.seed} > {output.fa}


# align
rule famsa:
    input:
        output+"{sp}/sequence/{seed}.fa"
    output:
        output+"{sp}/aligned/{seed}/{seed}_aln_famsa.fa"
    threads:
        2
    shell:
        "famsa -t {threads} {input} {output}"

rule mafft:
    input:
        output+"{sp}/sequence/{seed}.fa"
    output:
        output+"{sp}/aligned/{seed}/{seed}_aln_mafft.fa"
    shell:
        "mafft --auto {input} > {output}"

rule muscle:
    input:
        output+"{sp}/sequence/{seed}.fa"
    output:
        output+"{sp}/aligned/{seed}/{seed}_aln_muscle.fa"
    threads:
        4
    shell:
        "muscle -align {input} -output {output} -threads {threads}"

rule kalign:
    input:
        output+"{sp}/sequence/{seed}.fa"
    output:
        output+"{sp}/aligned/{seed}/{seed}_aln_kalign.fa"
    shell:
        "kalign -i {input} -o {output}"

# alns_method=['mafft', 'muscle', 'kalign']
alns_method=['famsa']


# rule mergealign:
#     input:
#         lambda wildcards: expand(output+"{sp}/aligned/{seed}/{seed}_aln_{method}.fa", 
#         sp=wildcards.sp, seed=wildcards.seed, method=alns_method)
#     output:
#         output+"{sp}/aligned/{seed}/{seed}_consensus.fa"
#     shell:
#         """
# input_dir=$(dirname {output})
# cp {input} {output}
# """
# python2 scripts/MergeAlign.py -a $input_dir -f {output}


# trim module
rule trim_aln:
    input:
        output+"{sp}/aligned/{seed}/{seed}_aln_famsa.fa"
    output:
        output+"{sp}/aligned/{seed}/{seed}_aln_trimmed.fa"
    shell:
        "trimal -gappyout -in {input} -out {output}"

# tree reconstruction
rule tree:
    input:
        output+"{sp}/aligned/{seed}/{seed}_aln_trimmed.fa"
    output:
        output+"{sp}/tree/{seed}.nwk"
    shell:
        "fasttree {input} > {output}"


rule root_tree:
    input:
        output+"{sp}/tree/{seed}.nwk"
    output:
        output+"{sp}/rooted_tree/rooted_{seed}.nwk"
    shell:
        """
nleaves=$(tr -cd ',' < {input} | wc -c)
if [ $nleaves -eq 0 ]; then
    cat {input} > {output}
else
    FastRoot.py -i {input} > {output}
fi
"""

# QC tree and alignment

# plot tree HGT
rule plot_tree:
    input:
        tree=output+"{sp}/rooted_tree/rooted_{seed}.nwk",
        aln=output+"{sp}/aligned/{seed}/{seed}_aln_famsa.fa",
    output:
        output+"{sp}/tree_plots/{seed}.png"
    params:
        tax=output+"{sp}/taxonomy/{sp}_tax.txt",
        self_seeds=output+"{sp}/seed_ids/{seed}_self.txt",
        taxid=lambda wcs: samples[wcs.sp][0]
    shell:
        """
Rscript scripts/plot_tree.R -i {input.tree} -a {input.aln} -o {output} \
 -t {params.tax} --selfid {params.taxid} -s {params.self_seeds}
"""

rule summarize_results:
    input:
        trees=get_trees,
        plots=get_tree_plots
    output:
        temp(output+"{sp}/hgt_trees/{sp}_trees.nwk")
    shell:
        "cat {input.trees} {input.plots} > {output}"
