"""
VCF Parser

A python script to parse and process Variant Call Format (VCF) files. It is designed
to handle the VCF files produced by the VarsScan2 algorithm

Created on Oct 7, 2017
@author: Samuli Eldfors
"""

# Import necessary libraries
import pandas as pd
import vcf

class VcfParser:
    """
    A class used to parse VCF files.

    Methods
    -------
    parse(vcf_file):
        Parses a VCF file and returns a pandas DataFrame with processed data.

    parse_standard_vcf(vcf_file):
        Parses a standard VCF file and returns a pandas DataFrame with processed data.

    test():
        A method for testing the functionality of the class.
    """

    def __init__(self):
        '''
        Constructor
        '''
        pass

    def parse(self, vcf_file):
        """
        Parse a VCF file. This function takes as input a path to a VCF file and returns a pandas DataFrame with processed data.

        Parameters:
        vcf_file (str): The path to the VCF file.

        Returns:
        df (pd.DataFrame): A DataFrame containing the parsed data.
        """

        # Initialize required variables
        handle = open(vcf_file, "r")
        vcf_reader = vcf.Reader(handle)
        data = []

        # Process each record in the VCF file
        for record in vcf_reader:
            for i in record.ALT:
                if record.CHROM.find("chr") == -1:
                    chrom = "chr" + record.CHROM
                else:
                    chrom = record.CHROM
                var_data = "%s_%s_%s>%s" % (chrom, record.POS, record.REF, i)

                id = record.ID
                normal_var_freq = float(record.INFO['AF1'])

                depth = int(record.INFO['DP'])
                somatic_status = int(record.INFO['SS'])
                somatic_p_value = float(record.INFO['SPV'])
                germline_p_value = float(record.INFO['GPV'])
                filter = record.FILTER

                normal_dp4 = normal_sample.data.DP4.split(",")
                normal_ref_ff = int(normal_dp4[0])
                normal_ref_rev = int(normal_dp4[1])
                normal_var_ff = int(normal_dp4[2])
                normal_var_rev = int(normal_dp4[3])

                try:
                    eff = record.INFO['EFF']
                    data_row = {"chrom": chrom, "pos": record.POS, "ref": record.REF, "var": record.ALT[0].sequence,
                                "var_data": var_data, "id": id, "normal_reads1": normal_reads1,
                                "normal_reads2": normal_reads2, "normal_var_freq": normal_var_freq, "depth": depth,
                                "somatic_status": somatic_status, "somatic_p_value": somatic_p_value,
                                "germline_p_value": germline_p_value, "filter": filter, 'eff': eff,
                                "normal_ref_ff": normal_ref_ff, "normal_ref_rev": normal_ref_rev,
                                "normal_var_ff": normal_var_ff, "normal_var_rev": normal_var_rev}

                    data.append(data_row)
                    snpeff_file = True
                except KeyError:

                    data_row = {"chrom": chrom, "pos": record.POS, "ref": record.REF, "var": record.ALT[0].sequence,
                                "var_data": var_data, "id": id, "normal_reads1": normal_reads1,
                                "normal_reads2": normal_reads2, "normal_var_freq": normal_var_freq, "depth": depth,
                                "somatic_status": somatic_status, "somatic_p_value": somatic_p_value,
                                "germline_p_value": germline_p_value, "filter": filter, "normal_ref_ff": normal_ref_ff,
                                "normal_ref_rev": normal_ref_rev, "normal_var_ff": normal_var_ff,
                                "normal_var_rev": normal_var_rev}

                    data.append(data_row)
                    snpeff_file = False

        # After processing all records, convert data list to pandas DataFrame and return
        df = pd.DataFrame(data)
        df.index = df.var_data

        # Additional processing
        if snpeff_file == True:
            df = df[['chrom', 'pos', 'ref', 'var', 'id', 'filter', 'eff',  'depth', 'germline_p_value', 'normal_reads1', 'normal_reads2', 'normal_var_freq', 'somatic_p_value', 'somatic_status', 'normal_ref_ff', 'normal_ref_rev', 'normal_var_ff', 'normal_var_rev']]
        else:
            df = df[['chrom', 'pos', 'ref', 'var', 'id', 'filter', 'depth', 'germline_p_value', 'normal_reads1', 'normal_reads2', 'normal_var_freq', 'somatic_p_value', 'somatic_status', 'normal_ref_ff', 'normal_ref_rev', 'normal_var_ff', 'normal_var_rev']]
        return df

    def parse_standard_vcf(self, vcf_file):
        """
        Parse a standard VCF file. This function takes as input a path to a VCF file and returns a pandas DataFrame with processed data.

        Parameters:
        vcf_file (str): The path to the VCF file.

        Returns:
        df (pd.DataFrame): A DataFrame containing the parsed data.
        """

        # Initialize required variables
        handle = open(vcf_file, "r")
        vcf_reader = vcf.Reader(handle)
        data = []

        # Process each record in the VCF file
        for record in vcf_reader:
            for i in record.ALT:
                if record.CHROM.find("chr") == -1:
                    chrom = "chr" + record.CHROM
                else:
                    chrom = record.CHROM
                var_data = "%s_%s_%s>%s" % (chrom, record.POS, record.REF, i)

                id = record.ID
                normal_var_freq = float(record.INFO['AF1'])

                depth = int(record.INFO['DP'])
                somatic_status = int(record.INFO['SS'])
                somatic_p_value = float(record.INFO['SPV'])
                germline_p_value = float(record.INFO['GPV'])
                filter = record.FILTER

                normal_dp4 = normal_sample.data.DP4.split(",")
                normal_ref_ff = int(normal_dp4[0])
                normal_ref_rev = int(normal_dp4[1])
                normal_var_ff = int(normal_dp4[2])
                normal_var_rev = int(normal_dp4[3])

                try:
                    eff = record.INFO['EFF']
                    data_row = {"chrom": chrom, "pos": record.POS, "ref": record.REF, "var": record.ALT[0].sequence,
                                "var_data": var_data, "id": id, "normal_reads1": normal_reads1,
                                "normal_reads2": normal_reads2, "normal_var_freq": normal_var_freq, "depth": depth,
                                "somatic_status": somatic_status, "somatic_p_value": somatic_p_value,
                                "germline_p_value": germline_p_value, "filter": filter, 'eff': eff,
                                "normal_ref_ff": normal_ref_ff, "normal_ref_rev": normal_ref_rev,
                                "normal_var_ff": normal_var_ff, "normal_var_rev": normal_var_rev}

                    data.append(data_row)
                    snpeff_file = True
                except KeyError:

                    data_row = {"chrom": chrom, "pos": record.POS, "ref": record.REF, "var": record.ALT[0].sequence,
                                "var_data": var_data, "id": id, "normal_reads1": normal_reads1,
                                "normal_reads2": normal_reads2, "normal_var_freq": normal_var_freq, "depth": depth,
                                "somatic_status": somatic_status, "somatic_p_value": somatic_p_value,
                                "germline_p_value": germline_p_value, "filter": filter, "normal_ref_ff": normal_ref_ff,
                                "normal_ref_rev": normal_ref_rev, "normal_var_ff": normal_var_ff,
                                "normal_var_rev": normal_var_rev}

                    data.append(data_row)
                    snpeff_file = False

        # After processing all records, convert data list to pandas DataFrame and return
        df = pd.DataFrame(data)
        df.index = df.var_data

        # Additional processing
        if snpeff_file == True:
            df = df[['chrom', 'pos', 'ref', 'var', 'id', 'filter', 'eff', 'depth', 'germline_p_value', 'normal_reads1',
                     'normal_reads2', 'normal_var_freq', 'somatic_p_value', 'somatic_status', 'normal_ref_ff',
                     'normal_ref_rev', 'normal_var_ff', 'normal_var_rev']]
        else:
            df = df[['chrom', 'pos', 'ref', 'var', 'id', 'filter', 'depth', 'germline_p_value', 'normal_reads1',
                     'normal_reads2', 'normal_var_freq', 'somatic_p_value', 'somatic_status', 'normal_ref_ff',
                     'normal_ref_rev', 'normal_var_ff', 'normal_var_rev']]
        return df

