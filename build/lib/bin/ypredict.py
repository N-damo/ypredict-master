#!/usr/bin/env python
import pandas as pd 
from collections import defaultdict
import json
import types
import copy_reg
import gzip
import os
from argparse import ArgumentParser
from multiprocessing import Pool

def get_args():
    parser = ArgumentParser(description = 'y haplogroup predict')
    parser.add_argument('-vcf','--vcf', help = 'input your vcffile' )
    parser.add_argument('-m','--map', default = 'map.json')
    parser.add_argument('-s','--special', default = 'hfspecial.xlsx')
    args = parser.parse_args()
    return args

            
def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)

class predict(object):

    def __init__(self, vcffile, haplogroup_dict, special):
        self.vcffile = vcffile
        self.haplogroup_dict = haplogroup_dict
        self.genotype = self.simple()
        self.name = self.genotype.keys()
        self.special = special
        self._multiprocess()


    def parse_vcf(self):
        assert os.path.exists(self.vcffile), 'the vcf file does not exist in the provided path, please check it'        
        with gzip.open(self.vcffile, 'r') as f:
            genotype = defaultdict(list)
            genotype2 = defaultdict(int)
            for line in f:
                if line.startswith(b'##'):
                    continue
                elif line.startswith(b'#CHROM'):
                    vcf = line.strip().split()
                    n = len(vcf[9:])
                else:
                    line = line.strip().split()
                    pos = line[1]
                    ref = line[3]
                    alt = line[4]
                    if pos in self.haplogroup_dict:
                        for j in range(n):
                            alleles = line[9+j].split(':')[0]   
                            haplogroup = self.haplogroup_dict[pos][0]                 
                            if alleles == '0':
                                alleles = ref
                                if alleles == self.haplogroup_dict[pos][1]:
                                    genotype[vcf[9+j]].append([pos, haplogroup, 'map'])
                                else:
                                    genotype[vcf[9+j]].append([pos, haplogroup, 'unmap'])
                            elif alleles == '1':
                                alleles = alt
                                if alleles == self.haplogroup_dict[pos][1]:
                                    genotype[vcf[9+j]].append([pos, haplogroup, 'map'])
                                else:
                                    genotype[vcf[9+j]].append([pos, haplogroup, 'unmap'])
                            else:
                                genotype2[vcf[9+j]] += 1

        for sample in genotype2:
            if genotype2[sample] > 50000:
                print ('{} like female sample and is excluding from the final result now'.format(sample))
                genotype.pop(sample)

        return genotype

    def simple(self):
        genotype = self.parse_vcf()
        genotype2 = {}
        for sample in genotype:
            d = defaultdict(list)
            for j in genotype[sample]:
                haplogroup = j[1] 
                state = j[2]
                d[haplogroup].append(state)  
            t = {}     
            for haplogroup in d:
                m = 0
                test_haplogroup_number = len(d[haplogroup])
                for states in d[haplogroup]:
                    if states == 'map':
                        m += 1
                if float(m) /test_haplogroup_number >= 0.2:
                    t[haplogroup] = 'T'
                else:
                    t[haplogroup] = 'F'
            genotype2[sample] = t

        return genotype2

    def get_parent(self, haplogroup):
        if haplogroup in list(self.special['haplogroup']):
            parent = self.special[self.special['haplogroup'] == haplogroup]['father'].values[0]
        elif haplogroup.endswith('~') and len(haplogroup) > 2:
            parent = haplogroup[:-2]
        else: 
            parent = haplogroup[:-1]
        return parent

    def get_score(self, sample):
        result = {}
        for haplogroup in self.genotype[sample]:
            ori = haplogroup      
            if self.genotype[sample][ori] == 'T':
                n_t = 1
                n_f = 0
                while haplogroup != 'Y':
                    parent = self.get_parent(haplogroup)
                    if parent in self.genotype[sample]:
                        quality = self.genotype[sample][parent]
                        if quality == 'T':
                            n_t += 1
                        else:
                            n_f += 1
                    else:
                        n_f += 1           
                    haplogroup = parent
                s1 = ((n_t**2.0)/(n_t + n_f))
                result[ori] = [s1, n_t]       
        target_haplogroup = max(result, key = result.get)
        rank = result[target_haplogroup][0]
        return target_haplogroup, rank, result

    def _multiprocess(self):
        pool = Pool(10)
        resultList = {}
        for sample in self.name:
            resultList[sample] = pool.apply_async(self.get_score, args = (sample,))
        pool.close()
        pool.join()
        results = {}
        with open('ypredict.txt', 'w') as f:
            for sample in resultList:
                target_haplogroup, rank, result = resultList[sample].get()
                results[sample] = result
                f.write('\t'.join([sample, target_haplogroup, str(rank)]) + '\n')
                print sample, target_haplogroup, rank

        df = pd.DataFrame.from_dict(results,orient='index').T
        df.to_csv('ystatistics.csv', header = True, index = True)


def mapload(name):
    with open(name, 'r') as f:
        for line in f:
            line = line.lstrip()
            haplogroup_dict = json.loads(line)
    return haplogroup_dict



if __name__ == '__main__':
    args = get_args()
    assert os.path.exists(args.special), 'the hfspecial.xlsx file does not exist, please check it'
    special = pd.read_excel(args.special, header = 0, sheet_name = 'sheet') 
    haplogroup_dict = mapload(args.map)
    predict(args.vcf, haplogroup_dict, special)


            
        
