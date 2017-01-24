"""
Run acceleration tests on each region in the conserved_bed.
"""

import os
import argparse
import itertools
from ete3 import Tree
from collections import OrderedDict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, TempFileTree, popenCatch, setLoggingFromOptions, logger


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hal', help='HAL alignment file.', required=True)
    parser.add_argument('--ref_genome', help='Reference genome.', required=True)
    parser.add_argument('--conserved_bed', help='Conserved BED file.', required=True)
    parser.add_argument('--conserved_model', help='Original conserved model', required=True)
    parser.add_argument('--non_conserved_model', help='Original nonconserved model.', required=True)
    parser.add_argument('--out_bed', help='output BED', required=True)
    parser.add_argument('--target_genomes', nargs='+', required=True, help='target genomes')
    parser.add_argument('--accelerated_genomes', nargs='+', required=True, help='target genomes')
    parser.add_argument('--single_copy_percent_cutoff', default=0.98, type=float,
                        help='Percent of region that must be single copy before we test it for acceleration')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


#### functions


def format_ratio(numerator, denominator):
    """
    Convenience function that converts two numbers, integer or no, to a ratio
    """
    if denominator == 0:
        return float("nan")
    return float(numerator) / denominator


def grouper(iterable, size):
    """
    http://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
    """
    it = iter(iterable)
    chunk = tuple(itertools.islice(it, size))
    while chunk:
        yield chunk
        chunk = tuple(itertools.islice(it, size))


def rename_model(target, model, accelerated_genomes):
    """Iteratively rename each ancestor of accelerated_genomes, walking down the tree to each ancestor"""
    new_model = os.path.join(target.getGlobalTempDir(), 'region_specific_conserved_subtree.mod')
    lines = open(model).readlines()
    t = Tree(lines[-1].split('TREE: ')[1], format=1)
    # this model may not have all of the genomes, if they were not aligned in this region
    accelerated_genomes = list(set(t.get_leaf_names()) & set(accelerated_genomes))
    if len(accelerated_genomes) > 1:
        anc = t.get_common_ancestor(accelerated_genomes)
        nodes = anc.get_descendants()
        leaves = anc.get_leaves()
        internal_nodes = [x for x in nodes if x not in leaves]
        for n in [anc] + internal_nodes:
            oldest_name = [x.name for x in n.get_children() if x.name != '1']
            if len(oldest_name) == 1:
                n.name = oldest_name[0] + '_Anc'
            else:
                n.name = '_'.join(oldest_name)
            with open(new_model, 'w') as outf:
                for l in lines[:-1]:
                    outf.write(l)
                outf.write('TREE: ' + t.write(format=1) + '\n')
            yield n.name, new_model
            n.name = '1'
    else: # only one accelerated genome here -- get common ancestor above will return root node
        with open(new_model, 'w') as outf:
            for l in lines[:-1]:
                outf.write(l)
            outf.write('TREE: ' + t.write(format=1) + '\n')
        yield accelerated_genomes[0], new_model


def conserved_model_contains_sufficient_outgroups(model, accelerated_genomes, outgroup_genomes, percent_outgroups=0.5):
    """makes sure that this region, when extracted, has at least percent_outgroups present"""
    lines = open(model).readlines()
    t = Tree(lines[-1].split('TREE: ')[1], format=1)
    outgroup_nodes = set(t.get_leaf_names()) - set(accelerated_genomes)
    return format_ratio(len(outgroup_nodes), len(outgroup_genomes)) >= percent_outgroups


#### main pipeline


def extract_maf_wrapper(target, args):
    """
    Main pipeline wrapper. Calls out to hal2maf once for each region in args.conserved_bed
    """
    accelerated_genomes = set(args.accelerated_genomes + [args.ref_genome])
    outgroup_genomes = set(args.target_genomes) - accelerated_genomes
    bed_recs = [x.split() for x in open(args.conserved_bed)]
    result_dir = target.getGlobalTempDir()
    result_tree = TempFileTree(result_dir)
    for chunk in grouper(bed_recs, 50):
        result_path = result_tree.getTempFile()
        target.addChildTargetFn(extract_and_calculate, args=(args, chunk, accelerated_genomes, outgroup_genomes, result_path))
    target.setFollowOnTargetFn(cat_results, args=(args, result_tree.listFiles()))


def extract_conserved_region(target, args, chrom, start, stop):
    """extract maf for a specific region and the target genomes"""
    maf_path = os.path.join(target.getGlobalTempDir(), 'region.maf')
    size = int(stop) - int(start)
    cmd = 'hal2maf --noAncestors --targetGenomes {} --refGenome {} --refSequence {} --start {} --length {} {} {}'
    cmd = cmd.format(','.join(args.target_genomes), args.ref_genome, chrom, start, size, args.hal, maf_path)
    system(cmd)
    target.logToMaster('extracted MAF for region {}:{}-{}'.format(chrom, start, stop))
    return maf_path


def write_region_bed(target, chrom, start, stop):
    region_bed = os.path.join(target.getGlobalTempDir(), 'region.bed')
    with open(region_bed, 'w') as outf:
        outf.write('\t'.join([chrom, start, stop]) + '\n')
    return region_bed


def calculate_region_specific_model(target, args, maf_path):
    """calculate the initial region-specific conserved model"""
    # phyloFit writes to the file *.mod in whatever the path set by --out-root is
    region_specific_conserved = os.path.join(target.getGlobalTempDir(), 'region_specific_conserved')
    cmd = 'phyloFit --init-model {} --scale-only --out-root {} {}'
    cmd = cmd.format(args.conserved_model, region_specific_conserved, maf_path)
    system(cmd)
    region_specific_conserved += '.mod'
    return region_specific_conserved


def test_ancestral_nodes(target, region_specific_conserved, accelerated_genomes, maf_path, region_bed, outf_handle):
    """Run LRT test on each ancestral node below the common ancestor of the accelerated genomes"""
    for anc_name, branch_model in rename_model(target, region_specific_conserved, accelerated_genomes):
        assert len(anc_name) > 0, (anc_name, branch_model, region_specific_conserved, accelerated_genomes)
        region_specific_accelerated = os.path.join(target.getGlobalTempDir(), 'region_specific_accelerated')
        cmd = 'phyloFit --init-model {} --scale-subtree {}:loss --out-root {} {}'
        cmd = cmd.format(branch_model, anc_name, region_specific_accelerated, maf_path)
        region_specific_accelerated += '.mod'
        system(cmd)
        cmd = 'phastOdds --output-bed --features {} --background-mods {} --feature-mods {} {}'
        cmd = cmd.format(region_bed, branch_model, region_specific_accelerated, maf_path)
        r = popenCatch(cmd)
        l = r.split()
        # discard the result if the test is not positive
        if int(l[-1]) > 0:
            l[-2] = anc_name
            outf_handle.write('\t'.join(l) + '\n')


def extract_and_calculate(target, args, chunk, accelerated_genomes, outgroup_genomes, result_path):
    """
    Runs the fitting, testing each ancestor that is under accelerated_genomes
    """
    with open(result_path, 'w') as outf_handle:
        for chrom, start, stop, score in chunk:
            # do not test a region if it is not at least single_copy_percent_cutoff
            if float(score) <= args.single_copy_percent_cutoff:
                target.logToMaster('locus {}:{}-{} has insufficient single copy percentage {}'.format(chrom, start, stop, score))
                continue
            # extract region MAF
            maf_path = extract_conserved_region(target, args, chrom, start, stop)
            # write the region BED needed by phastOdds
            region_bed = write_region_bed(target, chrom, start, stop)
            # calculate the initial region-specific conserved model
            region_specific_conserved = calculate_region_specific_model(target, args, maf_path)
            # should we continue?
            if not conserved_model_contains_sufficient_outgroups(region_specific_conserved, accelerated_genomes,
                                                                 outgroup_genomes):
                target.logToMaster('locus {}:{}-{} has insufficient outgroups'.format(chrom, start, stop))
                continue
            test_ancestral_nodes(target, region_specific_conserved, accelerated_genomes, maf_path, region_bed,
                                 outf_handle)


def cat_results(target, args, paths):
    """
    Concatenates final phastcons output into one gff file.
    """
    fofn = os.path.join(target.getGlobalTempDir(), 'fofn')
    with open(fofn, 'w') as outf:
        for p in paths:
            outf.write(p + '\n')
    tmp_p = os.path.join(target.getGlobalTempDir(), os.path.basename(args.out_bed) + '.tmp')
    cat_cmd = 'cat {} | xargs -n 50 cat > {}'.format(fofn, tmp_p)
    system(cat_cmd)
    os.rename(tmp_p, args.out_bed)


if __name__ == '__main__':
    from phast.run_acceleration_tests import *
    args = parse_args()
    setLoggingFromOptions(args)
    s = Stack(Target.makeTargetFn(extract_maf_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
