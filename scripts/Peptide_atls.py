import pandas as pd
from selenium.webdriver.chrome.options import Options
import requests, sys
from dataclasses import dataclass
from pyteomics import parser


import numpy as np




#prefs = {"download.default_directory": 'C:\\Users\\bairamkulov_dd\\PycharmProjects\\MRM_SRM_PARSER/data'}
import proteomeScoutAPI

chrome_options = Options()
chrome_options.add_argument("--disable-extensions")
chrome_options.add_argument('--start-maximized')
prefs = {
    'profile.default_content_settings.popups': 0,
    'download.default_directory': r'C:\Users\bairamkulov_dd\PycharmProjects\MRM_SRM_PARSER\csv\\',
    "directory_upgrade": True
}
chrome_options.add_experimental_option("prefs", prefs)

@dataclass(frozen = True)
class ParserBuilding:
    file: str

    @property
    def proteomeScout(self):
        """Return list of proteins modifications"""

    @property
    def uniprot(self, iD):
        """Return table of modification from uniprot"""

    @property
    def selectingByRules(self):
        """Return peptides selecting by Nisvidsky rules"""

    @property
    def getTableOfModification(self):
        """Return table of peptides modification"""

    @property
    def inSilicoDigestion(self):
        """in silico digestions fasta sequencing"""

    @property
    def selectTrypticPeptides(self):
        """Select tryptic peptides"""

    @property
    def selectingByModifications(self):
        """Selecting peptides by PTM"""


class SelectingPeptides(ParserBuilding):


    def uniprot(self, id):

        res_ptms = []
        res_type = []
        res_begin = []
        res_end = []
        res_peptide = []
        res_position = []
        print("######################################### Get proteins information from uniprot ###########################################")
        try:
            requestURL = f"https://www.ebi.ac.uk/proteins/api/proteomics-ptm/{id}"

            r = requests.get(requestURL, headers={"Accept": "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            responseBody = r.text

            df = pd.read_json(responseBody)
            for i in df['features']:
                res_type.append(i['type'])
                res_begin.append(i['begin'])
                res_end.append(i['end'])
            # res_xrefs.append(i['xref'])
                res_peptide.append(i['peptide'])
                res_peptide.append(i['unique'])

            df['type'] = pd.DataFrame(res_type)
            df['begin'] = pd.DataFrame(res_begin)
            df['end'] = pd.DataFrame(res_end)
        #df['xref'] = pd.DataFrame(res_xrefs)
            df['peptide'] = pd.DataFrame(res_peptide)
        #df['unique'] = pd.DataFrame(res_unique)
            for i in df['features']:
                for k in i['ptms']:
                    res_ptms.append(k['name'])
                    res_position.append(k['position'])

            df['modifications'] = pd.DataFrame(res_ptms)
            df['position_modification'] = pd.DataFrame(res_position)
            df = df.drop('features', axis = 1)
            print(df)
            return df
        except Exception as e:
            print(e)


    def getTableOfModification(self):
        res_key = []
        data = pd.read_excel(self.file)
        print("##################################### Get PTM #############################")
        for key in data['Unnamed: 0']:
            res_key.append(self.uniprot(key))
        modification_table =pd.concat(res_key)
        modification_table.to_csv('../data/modification_table.txt')
        return modification_table


    def inSilicoDigestion(self):
        table_list = []
        print("##################### In silico digestions stage #########################")
        modification_table = self.getTableOfModification()
        accession = pd.unique(modification_table['accession'])
        seq = pd.unique(modification_table['sequence'])
        for s in seq:
            for i in accession:
                peptide_df = pd.DataFrame(parser.xcleave(s, parser.expasy_rules['trypsin'], 0), columns = ['Position', 'Sequence'])
                peptide_df['accession'] = i
                table_list.append(peptide_df)

        table_sequence = pd.concat(table_list)
        table_sequence.to_csv('../data/table_of_digestion.txt')
        return  table_sequence


    def selectTrypticPeptides(self):
        drop_peptides = []
        peptides_table = self.inSilicoDigestion()
        print("######################## Selecting unique peptides stage ##########################")
        for i in peptides_table['Sequence']:
            if len(i) < 6:
                drop_peptides.append(i)

        peptides_table = peptides_table[peptides_table.Sequence.isin(drop_peptides) == False]
        peptides_table['Next value'] = peptides_table['Position'].shift(-1)
        peptides_table['Previous'] = peptides_table['Position'].shift(1)
        peptides_table['difference_next_and_position'] = peptides_table['Next value'] - peptides_table['Position']
        peptides_table['difference_previous_and_position'] = peptides_table['Position'] - peptides_table['Previous']
        peptides_table = peptides_table[(peptides_table['difference_previous_and_position'] != 3) & (peptides_table['difference_next_and_position'] != 3)]
        return peptides_table


    def selectingByRules(self):
        bad_sequencing = []
        peptideTable = self.selectTrypticPeptides()
        print("############################# Selecting peptides by rule #####################################")
        print(f"Peptides number before add rule:{len(peptideTable)}")
        print(peptideTable)
        for i in peptideTable['Sequence']:
            if 'KK' in i  or 'AA' in i:
                bad_sequencing.append(i)

        print(len(bad_sequencing))

        peptideTable = peptideTable[peptideTable.Sequence.isin(bad_sequencing) == False]
        print(f"Peptides number after add rule:{len(peptideTable)}")
        print(peptideTable)
        return peptideTable


    def selectingByModifications(self):
        bad_positions_start = []
        bad_position_end = []
        peptideTable = self.selectingByRules()
        ptm = pd.read_csv('../data/modification_table.txt')
        merge_table  = pd.merge(peptideTable,ptm, on = 'accession')
        try:
            for pos, nextv, begin, end in zip(merge_table['Position'], merge_table['Next value'], merge_table['begin'], merge_table['end']):
                if pos != pos  or nextv != nextv or begin != begin or end  != end:
                    print("The end of columns")
                for p in  range(int(begin), int(end)):
                    if pos <= p <= end:
                        bad_positions_start.append(pos)
                        bad_position_end.append(pos)
            print(list(set(bad_positions_start)))
            print(list(set(bad_position_end)))
            peptideTable = peptideTable[peptideTable.Position.isin(list(set(bad_positions_start))) == False]
            peptideTable = peptideTable[peptideTable.Position.isin(list(set(bad_position_end))) == False]
            peptideTable.to_csv("../data/peptide.tsv")
            print(peptideTable)
            return peptideTable

        except Exception as e:
            print(e)


    def proteomeScout(self):
        PTM_API = ProteomeScoutAPI('../data/all_moodifications.tsv')
        PTM_API.get_PTMs('P05155')







def main():
    print("Hello")
    cls = SelectingPeptides('../data/proteins.xlsx')
    cls.proteomeScout()





if __name__ == '__main__':
    main()
