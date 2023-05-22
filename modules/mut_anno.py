'''
Created on May 24, 2013
@author: Samuli Eldfors
'''

import pandas as pd
import eff_string_parser
import vcfParser
import vcf

class MutAnno(object):

    def __init__(self):
        self.version = "v1.5"
        self.vcfparser = vcfParser.VcfParser()
        self.eff_parser = eff_string_parser.EffStringParser()

    def annotate_snpeff_vcf2(self, snpeff_vcf, sample_code, results_dir, gene_expression, cgc_file,
                             ensembl_gene_descriptons, gene_exclusion_list_file, pan_cancer_file):

        out_dir = results_dir
        sample = sample_code
        vcf_file = snpeff_vcf
        gexp_file = gene_expression

        # Optional input
        if gexp_file != "":
            ge = pd.read_table(gexp_file)

        # Output files
        out_xls = out_dir + "/" + sample + "_ns_somatic_mutations_" + self.version + ".xls"
        cgc_out_xls = out_dir + "/" + sample + "_cgc_ns_somatic_mutations_" + self.version + ".xls"
        hc_snps_file = out_dir + "/" + sample + "_hc_snp_" + self.version + ".txt"
        filemaker_export_file = out_dir + "/" + sample + "_filemaker_" + version + ".tsv"
        unfiltered_mutations_file = out_dir + "/" + sample + "_unfiltered_somatic_calls_" + self.version + ".tsv"

        mut = self.vcfparser.parse(vcf_file)

        ens_gd = pd.read_table(ensembl_gene_descriptons)

        gef = gene_exclusion_list_file
        ge_list = pd.read_table(gef)

        mut["Effect"], mut["Effect_Impact"], mut["Functional_Class"], mut["Codon_Change"], mut["Amino_Acid_Change"], \
            mut[
                "Amino_Acid_Length"], mut["Gene_Name"], mut["Gene_Biotype"], mut["Coding"], mut["Transcript"], mut[
            "Exon"] = zip(
            *mut['eff'].map(self.eff_parser.get_first_eff))
        mut_anno = mut[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Effect', 'Effect_Impact', 'Amino_Acid_Change',
             'Transcript',
             'Exon', 'Codon_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq', 'tumor_reads1',
             'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]

        # Cancer_Gene_Census
        xls = pd.ExcelFile(cgc_file)
        cgc_dat = xls.parse('List', index_col=None, na_values=['NA'])
        cgc = cgc_dat[['Symbol', 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics']]

        mut_anno3 = mut_anno[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Exon', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Transcript', 'Amino_Acid_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]

        # Add CGC annotation
        mut_anno3 = mut_anno3.merge(cgc, how='left', left_on='Gene_Name', right_on="Symbol", suffixes=("", "_CGC"))

        # Add pan_cancer annotation

        pan_cancer_data = pd.read_table(pan_cancer_file, index_col=0)
        mut_anno3 = mut_anno3.merge(pan_cancer_data, how='left', left_on='Gene_Name', right_index=True)

        # Remove genes on exlcusion list
        mut_anno3 = mut_anno3.merge(ge_list, how='left', left_on='Gene_Name', right_on='Excluded_Gene')
        mut_anno3['Excluded_Gene'].fillna("Null", inplace='True')
        mut_anno3 = mut_anno3[mut_anno3.Excluded_Gene == "Null"]
        mut_anno3.Excluded_Gene.value_counts()

        mut_anno4 = mut_anno3[(mut_anno3.Effect_Impact == 'HIGH') | (mut_anno3.Effect_Impact == 'MODERATE')]

        # Add ensembl description
        mut_anno4 = mut_anno4.merge(ens_gd, how="left", left_on="Gene_Name", right_on="Associated Gene Name")

        # Add gene expression

        if gexp_file != "":
            ge_exists = 1
            mut_anno4 = mut_anno4.merge(ge, how="left", left_on="Gene_Name", right_on="gene_short_name")
        else:
            ge_exists = 0

        mut_anno4 = mut_anno4.sort('somatic_p_value', ascending=True)

        if ge_exists == 1:
            mut_anno4 = mut_anno4[
                ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Description', 'Effect', 'Effect_Impact',
                 'Codon_Change', 'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2',
                 'normal_var_freq', 'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value', 'FPKM',
                 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics', 'OncodriveFM-qval',
                 'Drivers', 'MuSIC', 'ActiveDriver', 'MutSig', 'OncodriveFM', 'OncodriveCLUST-qval', 'OncodriveCLUST']]
        elif ge_exists == 0:
            mut_anno4 = mut_anno4[
                ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Description', 'Effect', 'Effect_Impact',
                 'Codon_Change', 'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2',
                 'normal_var_freq', 'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value',
                 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics', 'OncodriveFM-qval',
                 'Drivers', 'MuSIC', 'ActiveDriver', 'MutSig', 'OncodriveFM', 'OncodriveCLUST-qval', 'OncodriveCLUST']]

        mut_anno4 = mut_anno4.drop_duplicates()
        mut_anno5 = mut_anno4[mut_anno4['id'].isnull()]

        # Write output files
        # Nonsynonymous mutations

        mut_anno5.to_excel(out_xls, index=False)
        mut_anno_cgc = mut_anno4[mut_anno4['Cancer Somatic Mut'] == 'yes']
        mut_anno_cgc.to_excel(cgc_out_xls)

        # Filemaker export format

        mut_anno5['sample'] = sample
        filemaker_data = mut_anno5[
            ['sample', 'chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]
        filemaker_data.to_csv(filemaker_export_file, sep="\t", index=False)

        # Unfiltered mutation calls

        mut_anno3_2 = mut_anno3[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Exon', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Transcript', 'Amino_Acid_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]
        mut_anno3_2.to_csv(unfiltered_mutations_file, sep="\t", index=False)

        # High confidence somatic SNPs

        snps = mut_anno3[mut_anno3['var'].map(lambda x: x in ['A', 'C', 'G', 'T']) & mut_anno3['ref'].map(
            lambda x: x in ['A', 'C', 'G', 'T'])]
        snps = snps[snps['id'].isnull()]
        snps.drop_duplicates(inplace=True)
        hc_snp = snps[snps['id'].isnull() & (snps['somatic_p_value'] < 0.01)]
        hc_snp.to_csv(hc_snps_file, index=False, sep="\t")

    def annotate_snpeff_vcf(self, snpeff_vcf, sample_code, results_dir, gene_expression, cgc_file,
                            ensembl_gene_descriptons, gene_exclusion_list_file, pan_cancer_file):

        out_dir = results_dir
        sample = sample_code
        vcf_file = snpeff_vcf
        gexp_file = gene_expression

        # Optional input
        if gexp_file != "":
            ge = pd.read_table(gexp_file)

        # Output files
        out_xls = out_dir + "/" + sample + "_ns_somatic_mutations_" + self.version + ".xls"
        cgc_out_xls = out_dir + "/" + sample + "_cgc_ns_somatic_mutations_" + self.version + ".xls"
        hc_snps_file = out_dir + "/" + sample + "_hc_snp_" + self.version + ".txt"
        filemaker_export_file = out_dir + "/" + sample + "_filemaker_" + self.version + ".tsv"
        unfiltered_mutations_file = out_dir + "/" + sample + "_unfiltered_somatic_calls_" + self.version + ".tsv"

        mut = self.vcfparser.parse(vcf_file)

        ens_gd = pd.read_table(ensembl_gene_descriptons)

        gef = gene_exclusion_list_file
        ge_list = pd.read_table(gef)

        mut["Effect"], mut["Effect_Impact"], mut["Functional_Class"], mut["Codon_Change"], mut["Amino_Acid_Change"], \
        mut[
            "Amino_Acid_Length"], mut["Gene_Name"], mut["Gene_Biotype"], mut["Coding"], mut["Transcript"], mut[
            "Exon"] = zip(
            *mut['eff'].map(self.eff_parser.get_high_eff))
        mut_anno = mut[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Effect', 'Effect_Impact', 'Amino_Acid_Change',
             'Transcript',
             'Exon', 'Codon_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq', 'tumor_reads1',
             'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]

        # Cancer_Gene_Census
        xls = pd.ExcelFile(cgc_file)
        cgc_dat = xls.parse('List', index_col=None, na_values=['NA'])
        cgc = cgc_dat[['Symbol', 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics']]

        mut_anno3 = mut_anno[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Exon', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Transcript', 'Amino_Acid_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]

        # Add CGC annotation
        mut_anno3 = mut_anno3.merge(cgc, how='left', left_on='Gene_Name', right_on="Symbol", suffixes=("", "_CGC"))

        # Add pan_cancer annotation

        pan_cancer_data = pd.read_table(pan_cancer_file, index_col=0)
        mut_anno3 = mut_anno3.merge(pan_cancer_data, how='left', left_on='Gene_Name', right_index=True)

        # Remove genes on exlcusion list
        mut_anno3 = mut_anno3.merge(ge_list, how='left', left_on='Gene_Name', right_on='Excluded_Gene')
        mut_anno3['Excluded_Gene'].fillna("Null", inplace='True')
        mut_anno3 = mut_anno3[mut_anno3.Excluded_Gene == "Null"]
        mut_anno3.Excluded_Gene.value_counts()

        mut_anno4 = mut_anno3[(mut_anno3.Effect_Impact == 'HIGH') | (mut_anno3.Effect_Impact == 'MODERATE')]

        # Add ensembl description
        mut_anno4 = mut_anno4.merge(ens_gd, how="left", left_on="Gene_Name", right_on="Associated Gene Name")

        # Add gene expression

        if gexp_file != "":
            ge_exists = 1
            mut_anno4 = mut_anno4.merge(ge, how="left", left_on="Gene_Name", right_on="gene_short_name")
        else:
            ge_exists = 0

        mut_anno4 = mut_anno4.sort('somatic_p_value', ascending=True)

        if ge_exists == 1:
            mut_anno4 = mut_anno4[
                ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Description', 'Effect', 'Effect_Impact',
                 'Codon_Change', 'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2',
                 'normal_var_freq', 'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value', 'FPKM',
                 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics', 'OncodriveFM-qval',
                 'Drivers', 'MuSIC', 'ActiveDriver', 'MutSig', 'OncodriveFM', 'OncodriveCLUST-qval', 'OncodriveCLUST']]
        elif ge_exists == 0:
            mut_anno4 = mut_anno4[
                ['chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Description', 'Effect', 'Effect_Impact',
                 'Codon_Change', 'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2',
                 'normal_var_freq', 'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value',
                 'Cancer Somatic Mut', 'Cancer Germline Mut', 'Cancer Molecular Genetics', 'OncodriveFM-qval',
                 'Drivers', 'MuSIC', 'ActiveDriver', 'MutSig', 'OncodriveFM', 'OncodriveCLUST-qval', 'OncodriveCLUST']]

        mut_anno4 = mut_anno4.drop_duplicates()
        mut_anno5 = mut_anno4[mut_anno4['id'].isnull()]

        # Write output files
        # Nonsynonymous mutations

        mut_anno5.to_excel(out_xls, index=False)
        mut_anno_cgc = mut_anno4[mut_anno4['Cancer Somatic Mut'] == 'yes']
        mut_anno_cgc.to_excel(cgc_out_xls)

        # Filemaker export format

        mut_anno5['sample'] = sample
        filemaker_data = mut_anno5[
            ['sample', 'chrom', 'pos', 'ref', 'var', 'id', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Exon', 'Transcript', 'Amino_Acid_Change', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]
        filemaker_data.to_csv(filemaker_export_file, sep="\t", index=False)

        # Unfiltered mutation calls

        mut_anno3_2 = mut_anno3[
            ['chrom', 'pos', 'ref', 'var', 'id', 'Exon', 'Gene_Name', 'Effect', 'Effect_Impact', 'Codon_Change',
             'Transcript', 'Amino_Acid_Change', 'depth', 'normal_reads1', 'normal_reads2', 'normal_var_freq',
             'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 'somatic_p_value']]
        mut_anno3_2.to_csv(unfiltered_mutations_file, sep="\t", index=False)

        # High confidence somatic SNPs

        snps = mut_anno3[mut_anno3['var'].map(lambda x: x in ['A', 'C', 'G', 'T']) & mut_anno3['ref'].map(
            lambda x: x in ['A', 'C', 'G', 'T'])]
        snps = snps[snps['id'].isnull()]
        snps.drop_duplicates(inplace=True)
        hc_snp = snps[snps['id'].isnull() & (snps['somatic_p_value'] < 0.01)]
        hc_snp.to_csv(hc_snps_file, index=False, sep="\t")

