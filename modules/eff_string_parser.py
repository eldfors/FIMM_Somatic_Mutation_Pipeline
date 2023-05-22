'''
Created on Jan 21, 2013
@author: Samuli Eldfors
'''

import eff

class EffStringParser(object):

    def parse_effs(self, eff_list):
        return eff_list

    def parse_eff(self, eff_string):
        eff_obj = eff.Eff()
        eff_obj.Effect, eff_data_str = eff_string.lstrip().split("(")
        eff_data = eff_data_str.rstrip(")").split("|")

        eff_obj.Effect_Impact = eff_data[0]
        eff_obj.Functional_Class = eff_data[1]
        eff_obj.Codon_Change = eff_data[2]
        eff_obj.Amino_Acid_Change = eff_data[3]
        eff_obj.Amino_Acid_Length = eff_data[4]
        eff_obj.Gene_Name = eff_data[5]
        eff_obj.Gene_Biotype = eff_data[6]
        eff_obj.Coding = eff_data[7]
        eff_obj.Transcript = eff_data[8]
        eff_obj.Exon = eff_data[9]

        return eff_obj

    def get_first_eff(self, eff_list):
        high_eff = self.parse_eff(eff_list[0])
        high_eff_info = high_eff.eff_data()
        return high_eff_info

    def get_high_eff(self, effs_string):
        eff_objs = {}
        effs_list = self.parse_effs(effs_string)

        for i in effs_list:
            eff_obj = self.parse_eff(i)
            eff_objs[eff_obj.severity] = eff_obj

        key_list = sorted(eff_objs.keys())
        high_eff = eff_objs[key_list[0]]
        high_eff_info = high_eff.eff_data()

        return high_eff_info

    def test(self):
        print "Multi EFF string parsing test\n"
        test_effs_string2 = 'missense_variant(MODERATE|MISSENSE|gAg/gTg|p.Glu1118Val/c.3353A>T|1553|TRPM2|protein_coding|CODING|ENST00000397932|22|T), intragenic_variant(MODIFIER|||n.44418447T>A||TRPM2-AS||NON_CODING|||T)'
        eff_data = self.get_first_eff(test_effs_string2)
        print eff_data

if __name__ == '__main__':
    a = EffStringParser()
    a.test()
