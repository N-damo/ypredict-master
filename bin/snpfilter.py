#!/usr/bin/env python
#coding: utf-8
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import subprocess
import argparse, os, gzip, json
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description = 'filter snp')
    parser.add_argument('-snp', '--snp', help = 'input snp file')
    parser.add_argument('-y', '--y', help = 'input Y fasta')
    args = parser.parse_args()
    return args



class snpfilter(object):

    def __init__(self, snp):
        self.snp = snp
        self.filter = self.filter_snp()
        
    def filter_snp(self):
        df = pd.read_csv(self.snp,header = 1)
        df.columns = ['name','haplogroup','other_name','rs','b37','b38','mutation']
        df = df[['haplogroup','b38','mutation']]
        s = df[df['mutation'].str.len() == 4]
        uniq = s.drop_duplicates(subset = "b38",keep = False)
        dup = s[s["b38"].duplicated(keep = False)]
        dup = dup.sort_values(by=['b38','mutation'])
        result = dup.drop_duplicates(subset = ['b38','mutation', 'haplogroup'],keep = False)
        r = result.drop_duplicates(subset = 'b38',keep = 'first')
        r = r['b38']
        f = dup[~dup['b38'].isin(r.values)]
        f2 = f.drop_duplicates(subset = 'b38',keep = 'first')
        final = uniq.append(f2)
        filter = final[final['mutation'].str.startswith(('A','T','C','G')) & final['mutation'].str.endswith(('A','T','C','G'))]
        filter.insert(2,'ancestral',filter['mutation'].str[0])
        filter.insert(3,'descendant',filter['mutation'].str[3])
        filter = filter.drop('mutation', axis = 1)
        filter = filter.dropna()
        filter = filter[filter['ancestral'] != filter['descendant']]
        filter = filter[~filter['haplogroup'].isin(['Notes', 'SNP'])]
        filter.to_csv('filter.csv',index=False, header = True)
        return filter


    def write_vcf(self):
        seq = SeqIO.read(args.y, 'fasta')
        result = []
        position = self.filter['b38']
        for i in position:
            result.append((i, seq.seq[int(i)-1]))
        result = pd.DataFrame(result)
        result.rename(columns={0:'b38',1:'REF'}, inplace=True)
        vcf = pd.merge(result, self.filter, on = 'b38')
        vcf1 = vcf[vcf['REF'] == vcf['ancestral']]
        vcf1.columns = ['POS', 'REF', 'HAPLOGROUP', 'INFO', 'ALT']
        vcf2 = vcf[vcf['REF'] == vcf['descendant']]
        vcf2 = vcf2.drop('descendant', axis = 1)
        vcf2.insert(4,'ALT', vcf2['ancestral'])
        vcf2.columns = ['POS', 'REF', 'HAPLOGROUP', 'INFO', 'ALT']
        vcf = pd.concat([vcf1, vcf2])
        vcf.insert(0,'#CHROM','chrY')
        vcf.insert(2,'ID','.')
        vcf.insert(5,'QUAL','.')
        vcf.insert(6,'FILTER','.')
        vcf.insert(8,'FORMAT','.')
        ref_vcf = vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'HAPLOGROUP']]
        ref_vcf['POS'] = ref_vcf['POS'].astype(int)
        ref_vcf = ref_vcf.sort_values('POS')
        ref_vcf.to_csv('vcf', header = True, index = False, sep = '\t')
        
        with open('vcfhead', "w") as f:
            f.write("##fileformat=VCFv4.2")
            f.write('\n')
            f.write("##reference=hg38")
            f.write('\n')
            f.write('##contig=<ID=chrY,length=57227415,assembly=B38,species="Homo    sapiens">')
            f.write('\n')
            f.write('##INFO=<ID=A,Number=1,Type=String,Description="Ancestral        Allele">')
            f.write('\n')
            f.write('##INFO=<ID=T,Number=1,Type=String,Description="Ancestral        Allele">')
            f.write('\n')
            f.write('##INFO=<ID=C,Number=1,Type=String,Description="Ancestral        Allele">')
            f.write('\n')
            f.write('##INFO=<ID=G,Number=1,Type=String,Description="Ancestral        Allele">')
            f.write('\n')

        subprocess.call("cat vcfhead vcf |bgzip >ref_vcf.gz", shell = True)
        subprocess.call('rm vcf', shell = True)

def hp(): 
    haplogroup_dict = {}
    assert os.path.exists('ref_vcf.gz'),'the ref_vcf.gz file does not exist, please check it'
    with gzip.open('ref_vcf.gz', 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                line = i.strip().split()
                pos = line[1]
                alt = line[4]
                ref = line[3]
                ancestral = line[7]
                if alt == ancestral:
                    effect = ref
                else:
                    effect = alt
                haplogroup = line[9].split(' ')[0] 
                haplogroup_dict[pos] = [haplogroup, effect]
                
    return haplogroup_dict


def mapout(haplogroup_dict, name):
    with open(name, 'w') as f:
        json.dump(haplogroup_dict, f, ensure_ascii=False)
        f.write('\n')

if __name__ == '__main__':
    args = get_args()
    snp = args.snp
    f = snpfilter(snp)
    f.write_vcf()
    haplogroup_dict = hp()
    mapout(haplogroup_dict, 'map.json')
