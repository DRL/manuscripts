#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
parse SnpEffVCF
"""

from __future__ import division
from itertools import izip
from os.path import isfile
import sys
import gzip
import itertools



def get_columns_from_vcf(infile):
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Nfld    Ro1_Mar149      Ro1_QC  Ro1_Ref Ro1_Ro19        SCRIR03 SCRIR04 SCRIR05
    columns = []
    with gzip.open(infile, 'rb') as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                columns = line.lstrip("#").rstrip("\n").split("\t")
                break
    return columns

def read_vcf(infile):
    with gzip.open(infile, 'rb') as fh:
        header, seqs = '', []
        i = 0
        for line in fh:
            if not line.startswith('#'):
                i += 1
                #if i == 1000:
                #    sys.exit('Done')
                fields = line.rstrip("\n").split("\t")
                record = {k : v for k, v in izip(columns, fields)}
                record['INFO'] = {field.split('=')[0] : field.split('=')[1] for field in record['INFO'].split(";") }
                if 'ANN' in record['INFO']:
                    record['INFO']['ANN'] = [value.split("|") for value in record['INFO'].get('ANN', '').split(",")]
                    record['ALLELES'] = {str(idx) : allele for idx, allele in enumerate((record['REF'] + "," + record['ALT']).split(','))}
                    for sample in samples:
                        record[sample] = {k : v for k, v in izip(record['FORMAT'].split(":"), record[sample].split(":")) }
                    yield record

def get_info(sample, genotype, annotations, gene_dict):
    for allele in genotype:
        if allele in annotations:
            gene = annotations[allele][3]
            if annotations[allele][1] == 'missense_variant':
                gene_dict[gene][sample]['N'] = gene_dict[gene][sample].get('N', 0) + 1
            elif annotations[allele][1] == 'synonymous_variant':
                gene_dict[gene][sample]['S'] = gene_dict[gene][sample].get('S', 0) + 1
            else:
                pass
            if 'ann' not in gene_dict[gene][sample]:
                gene_dict[gene][sample]['ann'] = []
            gene_dict[gene][sample]['ann'].append(annotations[allele][1])
    return gene_dict





if __name__ == '__main__':

    ALLOWED_ANN = [
        'coding_sequence_variant',
        'chromosome',
        'inframe_insertion',
        'disruptive_inframe_insertion',
        'inframe_deletion',
        'disruptive_inframe_deletion',
        'exon_variant',
        'exon_loss_variant',
        'frameshift_variant',
        'gene_variant',
        'intron_variant',
        'intragenic_variant',
        'missense_variant',
        'initiator_codon_variant',
        'stop_retained_variant',
        'rare_amino_acid_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'splice_region_variant',
        'stop_lost',
        'start_lost',
        'stop_gained',
        'synonymous_variant',
        'start_retained',
        'stop_retained_variant',
        'transcript_variant'
    ]
    ALL_ANN = {
        'coding_sequence_variant',
        'chromosome',
        'coding_sequence_variant',
        'inframe_insertion',
        'disruptive_inframe_insertion',
        'inframe_deletion',
        'disruptive_inframe_deletion',
        'downstream_gene_variant',
        'exon_variant',
        'exon_loss_variant',
        'frameshift_variant',
        'gene_variant',
        'intergenic_region',
        'conserved_intergenic_variant',
        'intragenic_variant',
        'intron_variant',
        'conserved_intron_variant',
        'miRNA',
        'missense_variant',
        'initiator_codon_variant',
        'stop_retained_variant',
        'rare_amino_acid_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'splice_region_variant',
        'splice_region_variant',
        'splice_region_variant',
        'stop_lost',
        '5_prime_UTR_premature',
        'start_lost',
        'stop_gained',
        'synonymous_variant',
        'start_retained',
        'stop_retained_variant',
        'transcript_variant',
        'regulatory_region_variant',
        'upstream_gene_variant',
        '3_prime_UTR_variant' ,
        '3_prime_UTR_truncation + exon_loss',
        '5_prime_UTR_variant' ,
        '5_prime_UTR_truncation + exon_loss_variant'  ,
        'sequence_feature + exon_loss_variant'
    }

    try:
        vcf_f = sys.argv[1]
    except:
        sys.exit("Usage : ./parse_snpeff.py VCF.gz > FILE.txt")

    if not isfile(vcf_f):
        sys.exit("Usage : ./parse_snpeff.py VCF.gz > FILE.txt")

    columns = get_columns_from_vcf(vcf_f)
    samples = columns[9:]

    singleton_sets = {}
    for ann in ALLOWED_ANN:
        singleton_sets[ann] = {}
        for set_size in range(1, len(samples) + 1):
            for subset in itertools.combinations(samples, set_size):
                singleton_sets[ann][frozenset(subset)] = 0
    type_count = {}
    for sample in samples:
        type_count[sample] = {'total' : {}, 'coding' : {}}
        type_count['all'] = {'total' : {}, 'coding' : {}}
    otherGtChar = {
        '.' : '.',
        '/' : '/',
        '|' : '|'
        }

    '''
    #TYPE    : The type of allele, either snp, mnp, ins, del, or complex.
    #DPB     : Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype"
    #MQM     : Mean mapping quality of observed alternate alleles
    #MQMR    : Mean mapping quality of observed reference alleles
    '''
    gene_dict = {}
    print "writing data.txt"
    with open(vcf_f + ".data.txt", 'w') as fh:
        fh.write("#CHROMOSOME\tPOS\tREF\tALT\tTYPE\tDPB\tMQM\tMQMR\tANNs\t%s\tGENE\tANN" % ("\t".join(sorted(samples))))
        fh.write("\n")
        for record in read_vcf(vcf_f):
            filtered_ANNs = {}
            filtered_ANNs_l = []
            annotation_strings = []
            for annotation_info in record['INFO']['ANN']:
                annotation_type = annotation_info[1]
                if '&' in annotation_type: # means there are two effects for the same allele
                    multiple_annotations_types = annotation_type.split("&")
                    allowed_annotation_types = "&".join([ann for ann in multiple_annotations_types if ann in ALLOWED_ANN])
                    if allowed_annotation_types:
                        if not allowed_annotation_types in singleton_sets:
                            singleton_sets[allowed_annotation_types] = {}
                            for set_size in range(1, len(samples) + 1):
                                for subset in itertools.combinations(samples, set_size):
                                    singleton_sets[allowed_annotation_types][frozenset(subset)] = 0
                        annotation_info[1] = allowed_annotation_types
                        filtered_ANNs[annotation_info[0]] = annotation_info
                        filtered_ANNs_l.append(annotation_info)
                else:
                    if annotation_info[1] in ALLOWED_ANN:
                        filtered_ANNs[annotation_info[0]] = annotation_info
                        filtered_ANNs_l.append(annotation_info)

                ann_set = {ann[1] : [] for ann in filtered_ANNs_l}
                genes =  {ann[3] for ann in filtered_ANNs_l}
                first = 1
                types = {allele : var_type for allele, var_type in izip(record['ALT'], record["INFO"]["TYPE"].split(","))}

            for sample in sorted(samples):
                genotype = {record['ALLELES'].get(char, ",") for char in record[sample]['GT']}
                for allele in genotype:
                    if allele in types:
                        type_count[sample]['total'][types[allele]] = type_count[sample]['total'].get(types[allele], 0) + 1
                        type_count['all']['total'][types[allele]] = type_count['all']['total'].get(types[allele], 0) + 1
            for annotation in filtered_ANNs.values():
                ann = annotation[1]
                gene = annotation[3]
                if gene not in gene_dict:
                    gene_dict[gene] = {sample : {"N" : 0, "S" : 0, "ann" : []} for sample in samples}
                if ann not in ann_set:
                    ann_set[ann] = []
                if (first):
                    fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (record["CHROM"], record["POS"], record["REF"], record["ALT"], record["INFO"]["TYPE"], record["INFO"]["DPB"], record["INFO"]["MQM"], record["INFO"]["MQMR"], len(filtered_ANNs)))
                    for sample in sorted(samples):

                        fh.write("\t%s" % record[sample]['GT'])
                        genotype = {record['ALLELES'].get(char, ",") for char in record[sample]['GT']}

                        for allele in genotype:
                            if allele in types:

                                type_count[sample]['coding'][types[allele]] = type_count[sample]['coding'].get(types[allele], 0) + 1
                                type_count['all']['coding'][types[allele]] = type_count['all']['coding'].get(types[allele], 0) + 1
                            if allele in filtered_ANNs:
                                ann_set[filtered_ANNs[allele][1]].append(sample)
                        gene_dict = get_info(sample, genotype, filtered_ANNs, gene_dict)
                    fh.write("\t%s" % ",".join(genes))
                    first = 0

                fh.write("\t%s" % "|".join([annotation[0], annotation[1], annotation[3], annotation[8], annotation[9], annotation[10]]))
            for ann_string in ann_set:
                if (frozenset(ann_set[ann_string])):
                    singleton_sets[ann_string][frozenset(ann_set[ann_string])] += 1
            if (filtered_ANNs):
                fh.write("\n")

    print "writing effect.txt"
    with open(vcf_f + ".gene_effect.txt", 'w') as fh:
        fh.write("#GENE\t%s" % "\t".join(sorted(samples)))
        fh.write("\n")
        for gene in sorted(gene_dict):
            fh.write(gene)
            for sample in sorted(gene_dict[gene]):
                fh.write("\t%s" % gene_dict[gene][sample]['ann'])
            fh.write("\n")

    print "writing misssense_synonymous.txt"
    with open(vcf_f + ".missesense_synonymous.txt", 'w') as fh:
        fh.write("#GENE\t%s" % "\t".join(sorted(samples)))
        fh.write("\n")
        for gene in sorted(gene_dict):
            fh.write(gene)
            for sample in sorted(gene_dict[gene]):
                dNdS = gene_dict[gene][sample]['N'] / gene_dict[gene][sample]['S'] if gene_dict[gene][sample]['S'] > 0 else 0.0
                fh.write("\t%s" % round(dNdS,2))
            fh.write("\n")

    print "writing private sets"
    with open(vcf_f + ".all_private_sets.txt", 'w') as fh:
        fh.write("#ANNOTATION\tCOUNT\tSETSIZE\tSUBSET\n")
        for annotation in sorted(singleton_sets):
            for subset in sorted(singleton_sets[annotation]):
                fh.write("%s" % annotation)
                fh.write("\t%s" % singleton_sets[annotation][subset])
                fh.write("\t%s" % len(subset))
                fh.write("\t%s" % (",".join(sorted(subset))))
                fh.write("\n")

    print "writing type counts"
    with open(vcf_f + ".type_counts.txt", 'w') as fh:
        fh.write("#SAMPLE\t")
        header = set()
        for sample in sorted(type_count):
            for category in sorted(type_count[sample]):
                for var_type in sorted(type_count[sample][category]):
                    header.add(category + "_" + var_type)
        fh.write("%s" % "\t".join(sorted(header)))
        fh.write("\n")
        for sample in sorted(type_count):
            fh.write("%s" % sample)
            for category in sorted(type_count[sample]):
                for var_type in sorted(type_count[sample][category]):
                    fh.write("\t%s" % type_count[sample][category][var_type])
            fh.write("\n")




