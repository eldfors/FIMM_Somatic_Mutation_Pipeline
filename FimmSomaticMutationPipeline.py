#!/Users/samuli/anaconda/bin/python

"""
Created on Jan 27, 2017
@author: Samuli Eldfors

This script generates and manipulates pileup files and uses several
external tools to analyze and annotate somatic mutations in tumor-normal sequencing data.
"""

import os
import os.path

import subprocess
import ConfigParser
import optparse

from modules import mut_anno, vcfeditor

version = "1.5"

parser = optparse.OptionParser(usage="%prog config.ini [options]", version="%prog " + version)

(options, args) = parser.parse_args()

if len(args) < 1:
    print("FimmSomaticMutationPipeline.py version: " + version + "\n")
    parser.error("No config file provided\n")

# Configuration
config = ConfigParser.ConfigParser()
config_file = args[0]
config.read(config_file)

# snpeff
snpeff_exe = config.get('snpeff', 'snpeff_exe')
snpeff_config = config.get('snpeff', 'snpeff_config')
snpeff_database = config.get('snpeff', 'snpeff_database')
snpsift_exe = config.get('snpeff', 'snpsift_exe')
dbsnp_vcf = config.get('snpeff', 'dbsnp_vcf')

# mpileup
samtools_exe = config.get('mpileup', 'samtools')
reference_genome = config.get('mpileup', 'reference_genome')
target_file = config.get('mpileup', 'target_file')
normal_bam = config.get('mpileup', 'normal_bam')
tumor_bam = config.get('mpileup', 'tumor_bam')

# varscan
varscan_exe = config.get('varscan', 'varscan_exe')
normal_pileup = config.get('varscan', 'normal_pileup')
tumor_pileup = config.get('varscan', 'tumor_pileup')
normal_purity = config.get('varscan', 'normal_purity')
tumor_purity = config.get('varscan', 'tumor_purity')
min_var_freq = config.get('varscan', 'min_var_freq')
strand_filter = config.get('varscan', 'strand_filter')
min_reads2 = config.get('varscan', 'min_reads2')
somatic_p_value = config.get('varscan', 'somatic_p_value')
min_coverage_normal = config.get('varscan', 'min_coverage_normal')
min_coverage_tumor = config.get('varscan', 'min_coverage_tumor')

# mut_anno
cgc_file = config.get('mut_anno', 'cgc_file')
ensembl_gene_descriptons = config.get('mut_anno', 'ensembl_gene_descriptons')
gene_exclusion_list = config.get('mut_anno', 'gene_exclusion_list')
pan_cancer_file = config.get('mut_anno', 'pan_cancer_file')
gene_expression = config.get('mut_anno', 'gene_expression_file')

# vcf2maf
vcf2maf = config.get('vcf2maf', 'vcf2maf')
vep_path = config.get('vcf2maf', 'vep_path')
vep_data = config.get('vcf2maf', 'vep_data')

# oncotator
oncotator_data_dir = config.get('oncotator', 'oncotator_data_dir')

# data
patient_id = config.get('mutannot', 'patient_id')
sample_id = config.get('mutannot', 'sample_id')

# Initialize the parser with the usage and version information
parser = optparse.OptionParser(usage = "%prog file.config [options]", version= "%prog " + version)

(options, args) = parser.parse_args()

if len(args) < 1:
    print ("varscanpipe.py version: " + version + "\n")
    parser.error("No config file provided\n")

print "Patient: " + patient_id
print "Sample: " + sample_id

analysis_root = os.getcwd()
log_file = os.path.join(analysis_root, sample_id + "_varscanpipe_log.txt")

# Functions
def create_results_dir():
    path = analysis_root + "/results"
    if not os.path.exists(path):
        os.mkdir(path)
    os.chdir(path)
    return path

def generate_pileup_files():
    if normal_pileup:
        print "Normal pileup found: " + normal_pileup
        npileup = normal_pileup

    else:
        print "Normal pileup not found. Generating pileup:"
        normal_pileup_out = os.getcwd() + os.sep + patient_id + "_normal.pileup"
        normal_mpileup_cmd = "%s mpileup -B -q1 -l %s -f %s %s > %s" % (
        samtools_exe, target_file, reference_genome, normal_bam, normal_pileup_out)
        print normal_mpileup_cmd
        process01 = subprocess.Popen(normal_mpileup_cmd, shell=True, executable="/bin/bash")
        process01.wait()
        npileup = normal_pileup_out

    if tumor_pileup:
        print "Tumor pileup found: " + tumor_pileup
        tpileup = tumor_pileup

    else:
        print "Tumor pileup not found. Generating pileup:"
        tumor_pileup_out = os.getcwd() + os.sep + patient_id + "_tumor.pileup"
        tumor_mpileup_cmd = "%s mpileup -B -q1 -l %s -f %s %s > %s" % (
        samtools_exe, target_file, reference_genome, tumor_bam, tumor_pileup_out)
        print tumor_mpileup_cmd
        process02 = subprocess.Popen(tumor_mpileup_cmd, shell=True, executable="/bin/bash")
        process02.wait()
        tpileup = tumor_pileup_out

    return npileup, tpileup

def normalize_vcf(in_vcf):
    out_vcf = os.path.realpath("%s.indel.norm.vcf" % (sample_id))

    cmd = "bcftools norm -D --check-ref w -f %s -o %s %s" % (reference_genome, out_vcf, in_vcf)
    print cmd

    process = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    process.wait()

    return out_vcf

def run_varscan(npileup, tpileup):
    varscan_cmd = "java -jar %s somatic %s %s %s --dream3-settings 1 --fpfilter --copynumber --tumor-purity %f --normal-purity %f --min-var-freq %f --strand-filter %i --output-vcf 1 --min-reads2 %i --somatic-p-value %f --min-coverage-normal %i --min-coverage-tumor %i 2> %s" % \
    (varscan_exe, npileup, tpileup, sample_id, float(tumor_purity), float(normal_purity), float(min_var_freq), int(strand_filter), int(min_reads2), float(somatic_p_value), int(min_coverage_normal), int(min_coverage_tumor), log_file)

    # Varscan somatic
    print varscan_cmd
    process1 = subprocess.Popen(varscan_cmd, shell=True, executable="/bin/bash")
    process1.wait()

    snp_file = os.path.realpath("%s.snp.vcf" % (sample_id))
    indel_file = os.path.realpath("%s.indel.vcf" % (sample_id))

    normalized_indel_file = normalize_vcf(indel_file)

    return (snp_file, normalized_indel_file)

def clean_vcf(vcf_file):

    outfile = analysis_root + os.sep + os.path.basename(vcf_file).replace(".vcf", ".ref.vcf")
    vcf_editor = vcfeditor.VcfEditor(vcf_file)
    vcf_editor.clean_vcf2(outfile)
    print "clean_vcf output written to: %s" % outfile

    os.chmod(outfile, 0777)
    return outfile

def somatic_vcf(vcf_file):

    outfile = analysis_root + os.sep + os.path.basename(vcf_file).replace(".ref.vcf", ".raw_somatic_calls.vcf")
    vcf_editor = vcfeditor.VcfEditor(vcf_file)
    vcf_editor.get_somatic(outfile)
    print "somatic_vcf output written to: %s" % outfile

    os.chmod(outfile, 0777)
    return outfile

def run_snpeff(vcf_file):

    path = analysis_root + "/snpeff"
    if not os.path.exists(path):
        os.mkdir(path)
    os.chdir(path)

    snpeff_file = os.path.realpath(path + os.sep + os.path.basename(vcf_file).replace(".vcf", ".annot.vcf"))
    stats_file = os.path.realpath(path + os.sep + os.path.basename(vcf_file) + ".snpeff_summary.html")

    cmd = "java -jar %s eff -ud 0 -formatEff %s -c %s %s -canon -o vcf -s %s -v | java -jar %s annotate %s - > %s 2>> %s" % (
    snpeff_exe, snpeff_database, snpeff_config, vcf_file, stats_file, snpsift_exe, dbsnp_vcf, snpeff_file, log_file)
    print cmd
    process = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    process.wait()

    os.chdir(analysis_root)
    return snpeff_file

def cat_vcf_files(vcf_files):

    vcf_snp, vcf_indel = (vcf_files)
    outfile = analysis_root + os.sep + os.path.basename(vcf_snp).replace(".snp.raw_somatic_calls.vcf",
                                                                         ".raw_somatic_calls.vcf")
    cmd = "vcf-concat %s %s | vcf-sort > %s" % (vcf_snp, vcf_indel, outfile)
    print cmd
    process = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    process.wait()

def main():
    # Generate pileup files
    npileup, tpileup = generate_pileup_files()

    # Run VarScan on pileup files
    varscan_vcf_files = run_varscan(npileup, tpileup)

    # Initialize an empty list for somatic VCF files
    somatic_vcf_files = []

    # Process each VCF file generated by VarScan
    for vcf_file in varscan_vcf_files:

        # Generate a somatic VCF file
        somatic_vcf_file = somatic_vcf(vcf_file)

        # Append the somatic VCF file to the list
        somatic_vcf_files.append(somatic_vcf_file)

    # Combine all somatic VCF files into one
    combined_somatic_vcf = cat_vcf_files(somatic_vcf_files)

    # Annotate the combined somatic VCF file with snpEff
    snpeff_file = run_snpeff(combined_somatic_vcf)

    # Create results directory
    results_dir = create_results_dir()

    # Initialize a mutation annotator
    annotator = mut_anno.MutAnno()

    # Annotate the snpEff VCF file
    annotator.annotate_snpeff_vcf(snpeff_file, sample_id, results_dir, gene_expression, cgc_file, ensembl_gene_descriptons, gene_exclusion_list, pan_cancer_file)

if __name__ == '__main__':
    main()

