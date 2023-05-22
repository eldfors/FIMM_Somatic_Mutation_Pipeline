'''
Created on Jan 21, 2013
@author: Samuli Eldfors
'''

class Eff(object):
    '''
    Variant Effect Information class
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.Effect = ""
        self.Effect_Impact = ""
        self.Functional_Class = ""
        self.Codon_Change = ""
        self.Amino_Acid_Change = ""
        self.Amino_Acid_Length = 0
        self.Gene_Name = ""
        self.Gene_Biotype = ""
        self.Coding = ""
        self.Transcript = ""
        self.Exon = ""

    def update_severity(self):
        eff_dict = {'SPLICE_SITE_ACCEPTOR' : 1, 'SPLICE_SITE_DONOR' : 2, 'START_LOST' : 3, 'EXON_DELETED' : 4, 'FRAME_SHIFT' : 5, 'STOP_GAINED' : 6, 'STOP_LOST' : 7, 'RARE_AMINO_ACID' : 8, 'NON_SYNONYMOUS_CODING' : 9, 'CODON_CHANGE' : 10, 'CODON_INSERTION' : 11, 'CODON_CHANGE_PLUS_CODON_INSERTION' : 12, 'CODON_DELETION' : 13, 'CODON_CHANGE_PLUS_CODON_DELETION' : 14, 'UTR_5_DELETED' : 15, 'UTR_3_DELETED' : 16, 'SYNONYMOUS_START' : 17, 'NON_SYNONYMOUS_START' : 18, 'START_GAINED' : 19, 'SYNONYMOUS_CODING' : 20, 'SYNONYMOUS_STOP' : 21, 'UTR_5_PRIME' : 22, 'UTR_3_PRIME' : 23, 'REGULATION' : 24, 'UPSTREAM' : 25, 'DOWNSTREAM' : 26, 'GENE' : 27,  'TRANSCRIPT' : 28, 'EXON' : 29, 'INTRON_CONSERVED' : 30, 'INTRON' : 31, 'INTRAGENIC' : 32, 'INTERGENIC' : 33, 'INTERGENIC_CONSERVED' : 34, 'NONE' : 35, 'CHROMOSOME' : 36, 'CUSTOM' : 37, 'CDS' : 38}
        self.severity = eff_dict.get(self.Effect, None)

    def eff_data(self):
        eff_data = (self.Effect, self.Effect_Impact, self.Functional_Class, self.Codon_Change, self.Amino_Acid_Change, self.Amino_Acid_Length, self.Gene_Name, self.Gene_Biotype, self.Coding, self.Transcript, self.Exon)
        print eff_data
        return eff_data

    def print_eff(self):
        result_line = "Effect: {} Effect Impact: {}\n".format(self.Effect, self.Effect_Impact)
        print result_line
