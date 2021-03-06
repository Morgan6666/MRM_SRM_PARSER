# ProteomeScoutAPI
#
# ===============================================================================
# ABOUT
# ===============================================================================
# Version 1.5
# 
# January 2016 
#
# By Alex Holehouse, Washington University in St. Louis 
# Contact alex.holehouse@gmail.com or https://github.com/alexholehouse
#
#
#
# ===============================================================================
# OVERVIEW
# ===============================================================================
# ProteomeScoutAPI is a Python module which can be used to connect to and parse
# ProteomeScout flatfiles. 
#
# Specifically, the goal of this module is to allow anyone to interact with 
# ProteomeScout data without the need to
#
# 1) Repeatedly query the ProteomeScout sever
#
# 2) Have any knowledge of SQL, or use an SQL-Python ORM
#
# 3) Facilitate rapid exploration of the ProteomeScout dataset
#
# ===============================================================================
# Usage
# ===============================================================================
# The general approach to usage is as follows
#
# 0) Import the proteomeScoutAPI module
#  
# 1) Create the API object by loading a flat file. ProteomeScout flat files can 
#    be obtained from https://proteomescout.wustl.edu/compendia
#
#    PTM_API = ProteomeScoutAPI('<flat filename goes here>')
#
# 2) Query the API object using the functions. Available functions are
#
#    PTM_API.get_phosphosites(ID)
#    PTM_API.get_PTMs(ID)
#    PTM_API.get_mutations(ID)
#
#    For more information on these functions I suggest reading the rest of the
#    source code, or once you've loaded an API object type
#
#    help(PTM_API.get_mutations)
#
# 3) The PTM_APU.uniqueKeys is a list of the unique accession numbers to provide
#    an easy way to loop over all the unique entries. NOTE THAT IDS HAVE REDUNDANCY
#    which is deliberate (and makes interfacing easy) but cycling over all the IDs
#    in the API object would be incorrect and lead to double counting
#
# ===============================================================================
# EXAMPLES
# ===============================================================================
# from proteomeScoutAPI import ProteomeScoutAPI 
#
# 
# ID = "O46631"
#  
# PTM_API = ProteomeScoutAPI("proteomescout_mammalia_20140831.tsv")
#
# PTM_API.get_mutations(ID)
#
# PTM_API.get_PTMs(ID)
#
# # the following loop prints all the phosphosites in all the proteins
# for ID in PTM_API.uniqueKeys():
#    print PTM_API.get_phosphosites(ID)
#
# ===============================================================================


# Exception if file is bad
from pbkdf2 import xrange


class BadProteomeScoutFile(Exception):
    pass

class ProteomeScoutAPIBug(Exception):
    pass


class ProteomeScoutAPI:


    def __init__(self,filename):
        """ 
        filename should be a ProteomeScout flatfile
        """        

        self.database={}
        self.uniqueKeys=[]

        # this will throw and exception if there's a problem with the file
        self.__checkFile(filename)
        
        self.__buildAPI(filename)


    def __checkFile(self, filename):
        """
        Internal function to check we (apparently) have a valid
        proteomeScout file
        """
        
        try:
            with open(filename, 'r') as f:
                first_line = f.readline()
                
            if not len(first_line.split("\t")) == 19:
                raise BadProteomeScoutFile("N/A")
            
                
        except:
            BadProteomeScoutFile("Invalid ProteomeScout flat file %s.\nFile is invalid or corrupted" % str(filename))
                
        

    def __buildAPI(self, datafile):

        # read in the file line by line
        with open(datafile) as f:
            content = f.readlines()
            

        
        # read the flatfile headers and parse
        # the contents
        headers_raw = content[0].split('\t')

        headers=[]
        for i in headers_raw:
            headers.append(i.strip())

        # remove the header line from the read in data
        content.pop(0)

        
        # for each line in the data
        for line in content:

            # split the record and get the canonical ID
            record          = line.split('\t')
            IDlist_raw      = record[1].split(";")
            ProteomeScoutID = record[0]
            IDlist = []
            for ID in IDlist_raw:
                IDlist.append(ID.strip())
                
            # add the first ID to the list of unique keys
            self.uniqueKeys.append(ProteomeScoutID)

            # now construct the object dictionary. Note we
            # have hardcoded the number of header columns here
            # though if in future versions of the ProteomeScout
            # dataset more columns are added you coulded extend this
            # here
            OBJ={}
            for i in xrange(1,19):
                OBJ[headers[i]] = record[i]


            # ALWAYS add the record via the ProteomeScout uniquekey
            self.database[ProteomeScoutID] = OBJ
        
            # for each other ID associated with the record, assign that accession-record
            # mapping assuming the record is not overwriting a more complete record which
            # already exists. 
            #
            # This means users can use whatever accession type they want to interact 
            # with the ProteomeScoutAPI, while recognizing that many different accessions
            # will be referring to the same protein
            for ID in IDlist:

                # A possible issue is the possibility that the same accession
                # points to several ProteomeScout records. Rather than merging records,
                # which can introduce problems, we've chosen to deal with this by defaulting
                # to the record with the largest number of modifications associated with it
                #
                # In most examples of record duplicaion there is clearly a 'major' and 'minor'
                # record. This occurs through insufficient cross-referencing in databases
                # outside of ProteomeScout. The 'minor' record is typically a subset of
                # the 'major' record
                if ID in self.database:

                    # if the newly found record has more PTMs associated with it than the 
                    # original record then overwrite 
                    if len(OBJ['modifications'].split(";")) > len(self.database[ID]['modifications'].split(";")):
                        self.database[ID] = OBJ

                else:                                        
                    self.database[ID] = OBJ

    def get_PTMs(self, ID):
        """
        Return all PTMs associated with the ID in question.

        POSTCONDITIONS:

        Returns a list of tuples of modifications
        [(position, residue, modification-type),...,]
        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications        

        """
                
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        mods = record["modifications"]
        
        mods_raw=mods.split(";")
        mods_clean =[]
        for i in mods_raw:
            tmp = i.strip()
            tmp = tmp.split("-")
            
            # append a tuple of (position, residue, type)
            mods_clean.append((tmp[0][1:], tmp[0][0], "-".join(tmp[1:])))
        return mods_clean
    
    def get_domains(self, ID, domain_type):
        """
        Return all domains associated with the ID in question.
        For pfam domains domain_type is 'pfam'
        For UniProt domains domain_type is 'uniprot'

        POSTCONDITIONS:

        Returns a list of tuples of domains 
        if there is a problem with the start and end position, these will be
        returned as -1
        [(domain_name, start_position, end_position),...,]
        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications        

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
        domain_type = domain_type.lower()
        if domain_type == 'pfam':
            doms = record["pfam_domains"]
        elif domain_type == 'uniprot':
            doms = record["uniprot_domains"]
        else:
            print("%s is an unrecognized domain type. Use 'pfam' or 'uniprot'"%(domain_type))
            return -2
        
        doms_raw=doms.split(";")
        doms_clean =[]
        for i in doms_raw:
            if i:
                tmp = i.strip()
                tmp = tmp.split(":")
                if len(tmp)>=2:
                    name, sites = tmp
                    tmp = sites.split("-")
                    doms_clean.append((name, tmp[0], tmp[1]))
                else:
                    print("ERROR: the domain did not match expected %s"%(i))
                    #doms_clean.append((tmp, -1, -1))
        return doms_clean
          

    def get_nearbyPTMs(self,ID,pos, window):
        """
        Return all PTMs associated with the ID in question that reside within
        +/-window, relative to the designated position (pos)

        POSTCONDITIONS:

        Returns a list of tuples of modifications
        [(position, residue, modification-type),...,]

        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications
        
        """
        mods = self.get_PTMs(ID)
        modsInWin = []
        if mods == -1:
            return -1
        elif len(mods)==0:
            return []
        else:
            for m in mods:
                site = int(m[0])
                if site >= pos-window and site <= pos+window:
                    modsInWin.append(m)

        return modsInWin



    def get_species(self,ID):
        """
        Return the species associated with the ID in question.
        POSTCONDITIONS:

        Returns a string of the species name

        Returns '-1' if unable to find the ID

        Returns '' (empty list) if no species


        """
        try:
            record = self.database[ID]
        except KeyError:
            return '-1'
        species = record["species"]
        return species
    
    def get_sequence(self, ID):
        """
        Return the sequence associated with the ID in question.
        POSTCONDITIONS:

        Returns a string of the sequence

        Returns '-1' if unable to find the ID

        Returns '' (empty list) if no sequence

        """
        try: 
            record = self.database[ID]
        except KeyError:
            return '-1'
        sequence = record["sequence"]
        return sequence

    
    def get_phosphosites(self,ID):
        """
        Return all phosphosites associated with the ID in question.

        POSTCONDITIONS:

        Returns a list of tuples of phosphosites
        [(position, residue, phosphosite-type),...,]
        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications        

        """

        mods = self.get_PTMs(ID)
        
        if mods == -1:
            return -1
        
        phospho=[]
        for mod in mods:
            if mod[2].find("Phospho") >= 0:
                phospho.append(mod)

        return phospho

    def get_mutations(self, ID):
        
        """
        Return all mutations associated with the ID in question.
        mutations = PTM_API.get_mutations(ID)
        
        POSTCONDITIONS:

        Returns a list of tuples of mutations 
        [(original residue, position, new residue, annotation),...,]
        
        Returns -1 if unable to find the ID
        Returns -2 if the number of mutations and annotations do not match

        Returns [] (empty list) if no mutations        

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        mutations = record["mutations"]
        mutations_ann = record["mutation_annotations"]
        if len(mutations) == 0:
            return []
        
        mutations_raw=mutations.split(";")
        mutations_clean=[]
        mutations_ann_raw = mutations_ann.split("|")
        if len(mutations_raw) != len(mutations_ann_raw):
            print("Error: Not the same number of annotations (%d) and mutations (%d)\n"%(len(mutations_ann_raw), len(mutations_raw)))
            return -2

        for idx, i in enumerate(mutations_raw):
            tmp = i.strip()            
            tmpArr = tmp.split(":")
            tmp = tmpArr[0];
            if len(tmpArr) == 2:
                label = tmpArr[1]
            else:
                label = ''
                
            # append a tuple of (position, residue, type)
            mutations_clean.append((tmp[1:-1], tmp[0], tmp[-1], label,
                mutations_ann_raw[idx].strip()))

        return mutations_clean


    def get_GO(self, ID):
        """
        Return all GO terms associated with the ID in question

        POSTCONDITIONS:
        
        Returns a list of GO Terms

        Returns a -1 if unable to find the ID

        Returns a [] (empty list) if no GO terms

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        GO_terms = record["GO_terms"]
        if len(GO_terms)==0:
            return []

        GO_termsArr = GO_terms.split(";")
        GO_terms_clean = []
        for i in GO_termsArr:
            GO_terms_clean.append(i.strip())

        return GO_terms_clean

    def get_accessions(self,ID):
        """
        Return a list of accessions associated with the protein

        POSTCONDITIONS:
        
        Returns a list of the protein's accessions

        Returns a -1 if unable to find the ID

        """

        try:
            record = self.database[ID]
        except KeyError:
            return -1

        raw_accessions   = record['accessions']
        split_accessions = raw_accessions.split(';')

        acc_list = []
        for acc in split_accessions:
            acc_list.append(acc.strip())

        return(acc_list)


    def get_evidence_per_site(self, ID):
        """
        Returns a list of PTMs with the number of different evidence sources
        included on a per-PTM basis

        """
        
        try:
            record = self.database[ID]
        except KeyError:
            return -1
    

        evidence_per_PTM = [len(x.split(',')) for x in record['evidence'].split(';')]
        PTMs = self.get_PTMs(ID)

        num_evidence = len(evidence_per_PTM)
        num_PTMs     = len(PTMs)
        
        if num_evidence != num_PTMs:
            raise ProteomeScoutAPIBug("While trying to determine evidence-per-PTM we obtained a mismatch in the data. Please report this along with the ID [%s] to alex.holehouse@wustl.edu" % (ID))

        PTMs_with_evidence = []
        for PTMIDX in xrange(0, num_PTMs):
            PTMs_with_evidence.append([PTMs[PTMIDX][0], PTMs[PTMIDX][1], PTMs[PTMIDX][2], evidence_per_PTM[PTMIDX]])

        return PTMs_with_evidence

    def get_protein_name(self, ID, first_name_only=True):
        """
        Returns a list of PTMs with the number of different evidence sources
        included on a per-PTM basis

        """
        try:
            record = self.database[ID]
        except KeyError:
            return -1

        if first_name_only:
            return record['protein_name'].split(';')[0]
        else:
            return record['protein_name']
        
            
                                    
            
            

            
            

        




    
                


        
        
