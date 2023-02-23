#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @File            : gtf.py
# @Date            : 2023/02/17 23:17:42
# @Author          : Baoting Nong (nong55@foxmail.com)
# @Link            : https://github.com/nongbaoting
# @Version         : 1.0.0
# @Description     : 

import os, sys, fire, re, gzip
from collections import defaultdict

def gtf_attr(attr):
    attributes = defaultdict(str)
    ats = attr.strip(' ').strip(';').split('; ')
    for at in ats:
        k, v = at.split(' ')[0:2]
        attributes[k] = v.strip('"')
    return attributes


class EXON:
    def __init__(self, cell):
        self.chrom = cell[0]
        self.source = cell[1]
        self.feature = cell[2]
        self.start = int(float(cell[3]))
        self.end = int(float(cell[4]))
        self.score = cell[5]
        self.strand = cell[6]
        self.frame = cell[7]
        self.attr = cell[8]
        self.gene_type = 'NA'
        attrDict = gtf_attr(self.attr)
        self.gene_id = attrDict['gene_id']
        self.gene_name = self.gene_id
        if 'gene_name' in attrDict:
            self.gene_name = attrDict['gene_name']
        else:
            self.gene_name = self.gene_id

        if 'transcript_id' in attrDict:
            self.transcript_id = attrDict['transcript_id']
        if 'gene_type' in attrDict:
            self.gene_type = attrDict['gene_type']
        if 'transcript_type' in attrDict:
            self.transcript_type = attrDict['transcript_type']
        else:
            self.transcript_type = 'NA'

    def to_gtf_line(self):
        line = '\t'.join([self.chrom, self.source, self.feature, str(self.start),
                        str(self.end), str(self.score), self.strand, self.frame, self.attr])
        return line

    def gtf_attrs_cat(self,attrsDict):
        it = []
        for key, val in attrsDict.items():
            it.append('{} "{}"'.format(key, val))
        attr = '; '.join(it) + ';'
        return attr


class Main:
    def to_gene_gtf(sef, gtf,outFi):
        geneDict = {}
        with open(gtf, 'r') as f:
            for line in f:
                if re.match('#', line):
                    continue
                if re.match('track', line):
                    continue
                cell = line.strip().split('\t')
                exon = EXON(cell)
                exon.feature = "gene"
                exon.attr = 'gene_id "' + exon.gene_id + '";'
                if exon.gene_id in geneDict:
                    gene = geneDict[exon.gene_id]
                    gene.start = min(gene.start, exon.start )
                    gene.end = max(gene.end, exon.end)
                    
                    geneDict[exon.gene_id] = gene
                else:
                    geneDict[exon.gene_id] = exon

        with open(outFi, 'w') as fo:
            for exon in geneDict.values():
                fo.write(exon.to_gtf_line() +'\n')

if __name__ == '__main__':
    fire.Fire(Main)