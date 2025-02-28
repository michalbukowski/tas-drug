# https://github.com/michalbukowski/tas-drug
# (c) 2023-2025 Michał Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# Department of Analytical Biochemistry, Faculty of Biochemistry, Biophysics and Biotechnology,
# Jagiellonian University, Krakow, Poland
# Distributed under GPL-3.0 license
# This notebook is a part of the TAS-Drug repository and is associated with the following study:
# Bukowski M, Banasik M, Chlebicka K, Bednarczyk K, Bonar E, Sokołowska D, Żądło T,
# Dubin G, Władyka B. (2025) Analysis of co-occurrence of type II toxin–antitoxin systems
# and antibiotic resistance determinants in Staphylococcus aureus. mSystems 0:e00957-24.
# https://doi.org/10.1128/msystems.00957-24

# 3 helper classes that organises information on coding sequences into a hierarchy
# of genes in loci and loci in assemblies (genomes)

class Gene:
    def __init__(self, name: str, start: int, end: int, qcovs: float, ppos: float):
        '''Requires 5 arguments describing a sequnece that
           is a translated BLAST hit:
           name: str    -- gene name
           start: int   -- start coordinate in the source DNA sequence
           end: int     -- end coordinate in the source DNA sequence
           qcovs: float -- protein query coverage
           ppos: float  -- protein sequence similarity
        '''
        self.name   = name
        self.start  = start
        self.end    = end
        self.qcovs  = qcovs
        self.ppos   = ppos
        
    def __repr__(self):
        return f'{self.name}({self.ppos})'

class Locus:
    def __init__(self, name: str, acc: str, strand: int, genes: list[Gene]):
        '''Requires 4 arguments describing a locus:
           name: str          -- locus name
           acc: str           -- accession number of the source DNA sequence
           strand: int        -- source DNA sequence strand
           genes: list[Gene]  -- a list of genes belonging to the locus
                                 (object of the Gene type)
        '''
        self.name   = name
        self.acc    = acc
        self.strand = strand
        self.genes  = genes
        
    def __repr__(self):
        return '/'.join(repr(gene) for gene in self.genes)

class Assembly:
    def __init__(self, asmid: str):
        '''Contains two empty lists that are to be populated by references
           to Locus type objects and store information on TA systems and
           drug resistance determinants found in a given assembly (genome).
           Requires 1 argument describing a blank assembly (genome):
           asmid: str -- assembly accession number
        '''
        self.asmid = asmid
        self.tas   = []
        self.drugs = []
    
    def __repr__(self):
        return ' '.join( repr(item) for item in self.tas + self.drugs)

