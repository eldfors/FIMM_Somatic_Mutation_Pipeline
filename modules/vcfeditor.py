# /apps/python272/bin/python

'''
Created on Apr 2, 2012

@author: Samuli Eldfors
'''

import vcf
import optparse


class VcfEditor(object):
    '''
    Class for editing VCF files.
    '''

    def __init__(self, vcf_file):
        '''
        Initialize VcfEditor with a vcf_file.
        '''
        self.vcf_file = vcf_file

    def clean_vcf(self, outfile):
        '''
        Clean the VCF file.
        '''
        vcf_reader = vcf.Reader(open(self.vcf_file, 'rb'))
        writer = vcf.Writer(open(outfile, 'w'), vcf_reader)

        for record in vcf_reader:
            record.CHROM = record.CHROM.lstrip("chr")
            if '/' in str(record.ALT[0].sequence):
                print "Skipped multiallele locus:", record
                continue
            writer.write_record(record)

        return outfile

    def clean_vcf2(self, outfile):
        '''
        Another method to clean the VCF file.
        '''
        vcf_reader = vcf.Reader(open(self.vcf_file, 'rb'))
        writer = vcf.Writer(open(outfile, 'w'), vcf_reader)

        for record in vcf_reader:
            record.CHROM = record.CHROM.lstrip("chr")

            if '/' in str(record.ALT[0].sequence):
                print "Skipped multiallele locus:", record
                continue

            alt_seq = str(record.ALT[0].sequence)
            if "+" in alt_seq:
                record.ALT[0] = record.REF + alt_seq.strip("+")
            elif "-" in alt_seq:
                deletion = alt_seq.strip("-")
                record.ALT[0] = record.REF
                record.REF = record.REF + deletion

            writer.write_record(record)

    def get_somatic(self, outfile, max_spv=1.0, min_depth=1):
        '''
        Get somatic variants from the VCF file.
        '''
        vcf_reader = vcf.Reader(open(self.vcf_file, 'rb'))
        writer = vcf.Writer(open(outfile, 'w'), vcf_reader)

        for record in vcf_reader:
            if record.INFO["SPV"] < max_spv and record.INFO["DP"] > min_depth and record.INFO["SS"] == '2':
                writer.write_record(record)

        return outfile

def main():
    '''
    Main function to run VCF editing.
    '''
    parser = optparse.OptionParser(usage="%prog infile.vcf outfile.vcf")
    parser.add_option("-s", action="store_true", dest="somatic", default=True)
    parser.add_option("-d", "--min_depth", action="store", dest="min_depth", default=1)
    parser.add_option("-p", "--max_spv", action="store", dest="max_spv", default=1.0)
    parser.add_option("-f", "--min_tumor_freq", action="store", dest="min_tumor_freq", default=0)

    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.error("No input files provided\n")

    infile = args[0]
    outfile = args[1]

