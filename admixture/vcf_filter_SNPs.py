#!/usr/bin/env python

import os,sys
import argparse
import gzip

def write_exclude_ids(in_vcf_filename, filtered_vcf, 
                      exclude_AT_CG, exclude_multi_ids, remove_dup_id,
                      remove_no_biAllelic, remove_dup_pos):
  """
  Filtering SNPs by various criteria 
  """
  
  ids = {}
  positions = {}
  
  with gzip.open(in_vcf_filename, 'rt') as in_fh, open(filtered_vcf, 'w') as out_vcf_fh:
    for line in in_fh:
        if line.startswith('#'):
            out_vcf_fh.write(line)
        else:
            fields = line.split("\t")
            cur_pos = fields[1]
            cur_id = fields[2]
            cur_A1 = fields[3]
            cur_A2 = fields[4]
            
            
            if ((remove_no_biAllelic) and (',' in cur_A2)):
            	print("Found multi-allelic %s, %s" % (cur_id, cur_A2))
            	continue
            
            if ((remove_dup_pos) and (cur_pos in positions)):
            	print("Found dup position %s" % (cur_pos))
            	continue
            
            if ((exclude_multi_ids) and (";" in cur_id)):
              print("Found multi id %s" % (cur_id))
              continue
            
            if ( (exclude_AT_CG) and
                 ( ( (cur_A1 == 'A') and (cur_A2 == 'T') ) or
                   ( (cur_A1 == 'T') and (cur_A2 == 'A') ) or
                   ( (cur_A1 == 'G') and (cur_A2 == 'C') ) or
                   ( (cur_A1 == 'C') and (cur_A2 == 'G') ) ) ):
              print("Found A/T C/G SNP %s" % (cur_id))
              continue
            
            if ((remove_dup_id) and (cur_id in ids)):
            	print("Found dup %s" % (cur_id))
            	continue
            
            positions[cur_pos] = None
            ids[cur_id] = None
            out_vcf_fh.write(line)


def main():
  parser = argparse.ArgumentParser("""
                                   Filter VCF file for various causes of errors in downstream (admixture) analysis.
                                   TODO input is .gz output is not 
                                   """)  
  
  parser.add_argument("in_vcf_filename", type=str, help= "Input vcf.gz file")
  parser.add_argument("filtered_vcf", type=str, help= "Output filttered VCF")
  
  parser.add_argument("-a", "--exclude_AT_CG", dest='exclude_AT_CG', default=True, type=bool,
                      help='Exclude A/T C/G snps that can be strand ambiguous ')
  
  parser.add_argument("-m", "--exclude_multi_ids", dest='exclude_multi_ids', default=True, type=bool,
                      help='Exclude ids that contain ";" useful for filtering SNPs with multiple ids')
  parser.add_argument("-d", "--remove_dup_id", dest='remove_dup_id', default=True, type=bool,
                      help='remove dup ids, keeps first')
  parser.add_argument("-b", "--remove_no_biAllelic", dest='remove_no_biAllelic', default=True, type=bool,
                      help='remove not bi-allelic')
  parser.add_argument("-p", "--remove_dup_pos", dest='remove_dup_pos', default=True, type=bool,
                      help='remove duplicated position, keep first')


  args = parser.parse_args()
    
  print('Starting....')
    
  write_exclude_ids(args.in_vcf_filename,args.filtered_vcf, 
                    args.exclude_AT_CG, args.exclude_multi_ids, args.remove_dup_id, 
                    args.remove_no_biAllelic, args.remove_dup_pos)
    
  print('Done!')

if __name__ == '__main__':
  main()


