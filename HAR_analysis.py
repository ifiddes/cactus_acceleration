# HAR analysis
# using the outgroup human alignment as a positive control for HAR identification

# select appris transcripts
hgsql -Ne 'select wgEncodeGencodeAttrsV24lift37.transcriptId, geneId, tag from wgEncodeGencodeTagV24lift37 join wgEncodeGencodeAttrsV24lift37 on wgEncodeGencodeTagV24lift37.transcriptId = wgEncodeGencodeAttrsV24lift37.transcriptId where tag like "appris%"' hg19 > appris.txt

# select basic genes
hgsql -Ne 'select name from wgEncodeGencodeBasicV24lift37' hg19 > basic_transcripts.txt

hgsql -Ne 'select geneId,geneName,geneType,transcriptId,transcriptType from wgEncodeGencodeAttrsV24lift37' hg19 > attrs.tsv

hgsql -Ne 'select * from wgEncodeGencodeBasicV24lift37' hg19 | cut -f 2- > wgEncodeGencodeAttrsV24lift37.gp

pyfasta flatten hg19.fa

# begin filtering appris transcripts
# first, pull down only name of protein coding genes + transcripts

from pycbio.bio.bio import *
from pycbio.bio.transcripts import *
from pycbio.bio.psl import *
from comparativeAnnotator.comp_lib.annotation_utils  import *

def get_coding_genes(attrs):
    return {x.split()[0] for x in open(attrs) if x.split()[2] == 'protein_coding'}


def get_coding_transcripts(attrs):
    return {x.split()[-2] for x in open(attrs) if x.split()[-1] == 'protein_coding'}


basic_transcripts = {x.rstrip() for x in open('basic_transcripts.txt')}
coding_genes = get_coding_genes('attrs.tsv')
coding_transcripts = get_coding_transcripts('attrs.tsv')
appris = {x.split()[0]: x.split() for x in open('appris.txt')}

# now filter appris set for tx and gene only coding
appris = {x: y for x, y in appris.iteritems() if x in coding_transcripts and y[1] in coding_genes and x in basic_transcripts}

gp = 'wgEncodeGencodeAttrsV24lift37.gp'
tx_dict = get_transcript_dict(gp)

# invert the appris dict so that it maps genes to all transcripts
# filter out any transcript with incomplete ends or not having proper starts/stops

stops = {'TGA', 'TAA', 'TAG'}
start = 'ATG'
seq_dict = get_sequence_dict('hg19.fa')

appris_by_gene = defaultdict(list)
for tx_id, (_, gene_id, tag) in appris.iteritems():
    if tx_id not in tx_dict:
        continue
    tx = tx_dict[tx_id]
    cds = tx.get_cds(seq_dict)
    if cds[:3] == start and cds[-3:] in stops and tx.cds_start_stat == 'cmpl' and tx.cds_end_stat == 'cmpl':
        appris_by_gene[gene_id].append([tx_id, tag])


good_tags = {'appris_principal', 'appris_principal_1', 'appris_principal_2', 'appris_principal_3', 'appris_principal_4', 'appris_principal_5'}
ordering = ['appris_principal', 'appris_principal_1', 'appris_principal_2', 'appris_principal_3', 'appris_principal_4', 'appris_principal_5', 'appris_alternative_1', 'appris_alternative_2']
ordering = {x: i for i, x in enumerate(ordering)}
one_to_one = {}
for gene_id, vals in appris_by_gene.iteritems():
    # remove the bad tags
    vals = [[tx_id, tag] for tx_id, tag in vals if tag in good_tags]
    if len(vals) == 0:
        continue
    tx_ids, tags = zip(*vals)
    best_tag = sorted(tags, key=lambda x: ordering[x])[0]
    best_txs = [tx_id for tx_id, tag in vals if tag == best_tag]
    if len(best_txs) > 1:
        # pick based on which is longest
        tx_recs = [tx_dict[x] for x in best_txs if x in tx_dict]
        # for human we have to deal with chrY, which was filtered out of tx dict
        if len(tx_recs) == 0:
            continue
        best_tx = sorted(tx_recs, key=lambda x: -len(x))[0].name
    else:
        best_tx = best_txs[0]
    one_to_one[gene_id] = best_tx

# write out the BEDs

with open('cleaned_cds.bed', 'w') as outf:
    for tx in one_to_one.itervalues():
        tx = tx_dict[tx]
        outf.write('\t'.join(map(str, tx.get_bed())) + '\n')


# bash code, making use of the human outgroup alignment for this.

genePredToBed wgEncodeGencodeAttrsV24lift37.gp wgEncodeGencodeAttrsV24lift37.bed

workDir=$PWD
halFile=/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1509/cactus/1509_outgroups.hal

# genomes
srcOrg=hg19
targetGenomes="SPRET_EiJ PWK_PhJ CAST_EiJ WSB_EiJ A_J CAROLI_EiJ Pahari_EiJ Rattus micOch1 jacJac1 oryCun2 hg19 panTro4 gorGor3 ponAbe2 rheMac3 oviAri3 bosTau8 canFam3 felCat8 loxAfr3"
acceleratedGenomes="hg19 panTro4 gorGor3 ponAbe2"

cdsBed=cleaned_cds.bed
fourDSites=$workDir/4d_sites.bed


export PYTHONPATH=/hive/users/ifiddes/comparativeAnnotator

# target metrics provided by David Thybert
targetCoverage=0.3
expectedLength=45


# output conserved BED
humanConservedBed=${workDir}/human_conserved_regions.bed
# output conserved wiggle
humanConservedWig=${workDir}/human_conserved_scores.wig

# conserved BED after merging
humanMergedConservedBed=${workDir}/human_merged_conserved_regions.bed
singleCopyScoredBed=${workDir}/human_single_copy_scored_regions.bed

# output conserved model to modify
consModel=${workDir}/conserved.mod
nonconsModel=${workDir}/non_conserved.mod

# LRT output BED
rawLrtBed=${workDir}/lrt_results.bed
# filtered BED for being uniquely mappable and having a positive LRT
filteredLrtResults=${workDir}/filtered_lrt_results.bed
# hits within 10kb of genes
lrtGeneHits=${workDir}/lrt_gene_hits.bed

# we must subset the HAL on juggernaut/kolossus because of a weird hdf5 error. We will save the extracted alignments
# in case we need to run this again.
subsetDir=${workDir}/subset_alignments

# jobTree needs a place to store the jobs. We can delete this once it is completed.
subsetTreeDir=${workDir}/subset_jobTree
phastTreeDir=${workDir}/phast_jobTree
filterTreeDir=${workDir}/single_copy_jobTree
lrtTreeDir=${workDir}/lrt_jobTree


codeDir=/hive/users/ifiddes/comparativeAnnotator

# extract 4d sites
hal4dExtract ${halFile} ${srcOrg} ${cdsBed} ${fourDSites}

# train neutral model
neutralModel=${workDir}/neutral_model.mod
python ${codeDir}/hal/phyloP/halPhyloPTrain.py --numProc 32 --noAncestors --no4d ${halFile} ${srcOrg} ${fourDSites} ${neutralModel} \
--targetGenomes SPRET_EiJ PWK_PhJ CAST_EiJ WSB_EiJ C57B6J A_J CAROLI_EiJ Pahari_EiJ Rattus micOch1 jacJac1 oryCun2 hg19 panTro4 gorGor3 ponAbe2 rheMac3 oviAri3 bosTau8 canFam3 felCat8 loxAfr3 \
--tree "((((((((((SPRET_EiJ:0.002,(PWK_PhJ:0.001,(CAST_EiJ:0.001,(WSB_EiJ:1e-05,(C57B6J:3e-06,A_J:4e-06)1:2e-06)1:0.0001)1:1e-06)1:1e-06)1:0.015,CAROLI_EiJ:0.02)1:0.01,Pahari_EiJ:0.03)1:0.02,Rattus:0.013)1:0.065,micOch1:0.15)1:0.117,jacJac1:0.1859)1:0.07,oryCun2:0.21)1:0.01,((((hg19:0.00642915,panTro4:0.00638042)1:0.00217637,gorGor3:0.00882142)1:0.00935116,ponAbe2:0.0185056)1:0.00440069,rheMac3:0.007)1:0.1)1:0.02,((oviAri3:0.019,bosTau8:0.0506)1:0.17,(canFam3:0.11,felCat8:0.08)1:0.06)1:0.02)1:0.02,loxAfr3:0.15);"


nice -n 19 python ${codeDir}/phast/phast_subset.py ${halFile} ${srcOrg} ${neutralModel} ${subsetDir} \
--target-genomes ${targetGenomes} --batchSystem singleMachine --maxThreads 40 --jobTree ${subsetTreeDir}

# use the output of this subset as input to the phastCons pipeline
# these can be cluster jobs, 8gb of RAM is fine
# limit to 100 cpus because these jobs can be >1hr
python ${codeDir}/phast/phastcons.py ${halFile} ${srcOrg} ${neutralModel} ${humanConservedBed} ${humanConservedWig} \
--output-conserved-model ${consModel} --output-nonconserved-model ${nonconsModel} --pre-extracted ${subsetDir} \
--target-genomes ${targetGenomes} --target-coverage ${targetCoverage} --expected-length ${expectedLength} \
--jobTree ${phastTreeDir} --batchSystem parasol --parasolCommand /cluster/home/ifiddes/bin/remparasol \
--defaultMemory 8589934591 --maxCpu 100

# merge nearby conserved elements
bedtools merge -d 50 -i ${humanConservedBed} | awk '($3-$2) >= 400' > ${humanMergedConservedBed}

# score these for overlapping with single copy regions
python ${codeDir}/phast/find_single_copy_regions.py --jobTree ${filterTreeDir} --batchSystem parasol --parasolCommand ~/bin/remparasol \
--hal ${halFile} --ref_genome ${srcOrg} --conserved_bed ${humanMergedConservedBed} --out_bed ${singleCopyScoredBed}

# run the likelihood ratio test
# has to be singleMachine due to hal2maf problems
python ${codeDir}/phast/run_acceleration_tests.py --jobTree ${lrtTreeDir} --batchSystem singleMachine --maxThreads 40 \
--accelerated_genomes ${acceleratedGenomes} --target_genomes ${targetGenomes} \
--hal ${halFile} --ref_genome ${srcOrg} --conserved_bed ${singleCopyScoredBed} --conserved_model ${consModel} \
--non_conserved_model ${nonconsModel} --out_bed ${rawLrtBed}

# extract the highest branch for each, requiring a score difference of at least 10
python ${codeDir}/phast/filter_acceleration_results.py --cutoff 10 ${rawLrtBed} > ${filteredLrtResults}
bedSort ${filteredLrtResults} ${filteredLrtResults}

# find if any of these are within 10kb of any reference annotations
referenceBed=/hive/users/ifiddes/ihategit/pipeline/gencode_V24/human_with_pseudogenes.bed
bedtools closest -a ${filteredLrtResults} -b ${referenceBed} -d | awk '$NF <= 10000' | cut -d $'\t' -f 1,2,3,4,5 | sort | uniq > ${lrtGeneHits}

# we received a HAR locus BED from Katie Pollard. hg19 coordinates!
bedtools intersect -a ${filteredLrtResults} -b nchaes_merged_hg19.bed -u > overlap.bed
# not very good, 77 hits