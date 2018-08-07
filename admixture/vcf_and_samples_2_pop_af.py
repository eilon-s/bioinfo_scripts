#!/usr/bin/env python


import os,sys
import argparse
import gzip
import pandas as pd
import numpy as np

from pysam import VariantFile
from collections import defaultdict


def get_af(in_vcf_filename,in_smaples_pop_filename, output_af_filename):
  """
  Getting allele frequencies for allele frequencies in pop
  """
  
  # getting samples
  samples_df = pd.read_table(in_smaples_pop_filename)
  pop_sample_dict = defaultdict(list)
  for pop in samples_df['pop'].unique():
    pop_sample_dict[pop] = samples_df['sample'][samples_df['pop'] == pop].tolist()  

  #parsing VCF and wrting AF table
  out_header = '\t'.join(('CHROM','POS','ID','REF','ALT')) + '\t'.join(pop_sample_dict.keys())
  
  line_cnt = 0
  with open(output_af_filename, 'w') as out_fh, VariantFile(in_vcf_filename) as in_vcf:
    out_fh.write(out_header + '\n')
    
    for rec in in_vcf.fetch():
        line_cnt += 1
        if (line_cnt % 10000 == 0):
          print("parsing variant %d" % (line_cnt))
      
      
        cur_chrom = rec.chrom
        cur_pos   = rec.pos
        cur_ref   = rec.ref
        cur_alt   = rec.alts[0]
        cur_id    = rec.id
        cur_rec_info_str = '\t'.join([cur_chrom,str(cur_pos),cur_id,cur_ref,cur_alt])
        out_fh.write(cur_rec_info_str)
        
        for pop in pop_sample_dict.keys():
            ref_cnt = 0.0
            alt_cnt = 0.0
            for samp in pop_sample_dict[pop]:
                if 'GT' in rec.samples[samp].keys():
                    cur_gt = rec.samples[samp]['GT']
                    if(None not in cur_gt):
                        ref_cnt += cur_gt.count(0)
                        alt_cnt += cur_gt.count(1)
                else:
                    raise ValueError("no GT field in %s" % (cur_rec_info_str))
            
            # considering just REF and first ALT
            if ((ref_cnt + alt_cnt) > 0):
                cur_af = alt_cnt / (ref_cnt + alt_cnt)
            else:
                cur_af = np.nan # 'nan' string will be written. TODO - is ok?
            out_fh.write('\t' + str(cur_af))
        out_fh.write('\n')  



def main():
  
  parser = argparse.ArgumentParser("Get AF for all populations")  
    
  parser.add_argument("in_vcf_filename", type=str, help= "Input vcf.gz file (for example: ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz)") # ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  parser.add_argument("in_smaples_pop_filename", type=str, help= "samples pop. file (for example: integrated_call_samples_v3.20130502.ALL.panel)") # 
  parser.add_argument("output_af_filename", type=str, help= "Output file of AF by populations. NOTE - if not bi-allelic will consider only ref and first alt alleles")
  
  args = parser.parse_args()
    
  print('Starting....')
  
  
    
  get_af(args.in_vcf_filename,args.in_smaples_pop_filename, args.output_af_filename)
    
  print('Done!')

if __name__ == '__main__':
  main()
