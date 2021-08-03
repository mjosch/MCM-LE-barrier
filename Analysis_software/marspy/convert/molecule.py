import re

import numpy as np
import scyjava as sc
import seaborn as sns
from scyjava.convert._pandas import table_to_pandas


class Molecule:

    def __init__(self, uid, archive):
        self.uid = uid
        self.archive = archive
        self.meta_uid = self.archive.get(self.uid).getMetadataUID()
        self.params = dict(sc.to_python(self.archive.get(self.uid).getParameters()))
        self.tags = list(sc.to_python(self.archive.get(self.uid).getTags()))
        self.regions = list(sc.to_python(self.archive.get(self.uid).getRegionNames()))
        # don't keep entire df in memory for future versions
        self.df = table_to_pandas(self.archive.get(self.uid).getTable())
        self.seg_dfs = None

    def __str__(self):
        return f'Greetings from Molecule {self.uid}.'

    def __len__(self):
        return len(self.df)

    def __del__(self):
        pass
        print(f'Molecule {self.uid} deleted.')


class SingleMolecule(Molecule):

    def __init__(self, uid, protein, archive):
        Molecule.__init__(self, uid, archive)
        self.protein = protein


class DnaMolecule(Molecule):

    def __init__(self, uid, proteins, archive):
        Molecule.__init__(self, uid, archive)

        # DnaMolecule specific attributes
        self.proteins = {protein: 0 for protein in proteins}

        # grab keys from dict
        for protein in self.proteins:
            # protein specific prefixes with nomenclature protein_prefixes:
            (exec(
                f"self.{protein + '_prefixes'} = set(re.findall(f'{protein}_\d+_',' '.join(self.df.columns)))"))

            # Store number of molecules based off of actual dataTable headers
            self.proteins[protein] = len(set(re.findall(f'{protein}_\d+_', ' '.join(self.df.columns))))

        # generate prefixes based union of protein_prefixes
        self.prefixes = set()
        for protein in self.proteins.keys():
            for prefix in getattr(self, f'{protein}_prefixes'):
                self.prefixes.add(prefix)

        # region objects for DnaMolecules
        self.regions = list()
        # all region names
        for region_name in sc.to_python(self.archive.get(self.uid).getRegionNames()):
            _region = self.archive.get(self.uid).getRegion(region_name)
            match_prefix = None
            # separate prefix from column name
            for prefix in self.prefixes:
                if re.match(prefix, _region.getColumn()):
                    # correct prefix found
                    match_prefix = prefix
                    break

            # append Region object to molecule regions
            self.regions.append(Region(uid=self.uid,
                                       name=region_name,
                                       start=_region.getStart(),
                                       end=_region.getEnd(),
                                       prefix=match_prefix,
                                       column=_region.getColumn().split(match_prefix)[-1]))

    def calc_length_dna(self):
        """
        Calculates the Molecule's DNA length in px.
        """
        return np.sqrt((self.params['Dna_Bottom_X2'] - self.params['Dna_Top_X1']) ** 2 +
                       (self.params['Dna_Bottom_Y2'] - self.params['Dna_Top_Y1']) ** 2)

    def plot(self):
        for prefix in self.prefixes:
            try:
                sns.lineplot(x=prefix + 'Time_(s)', y=prefix + 'Position_on_DNA', data=self.df)
            except ValueError:
                continue


class Region:
    """
    Region object holding following attributes: name, start, end and molecule column.
    """

    def __init__(self, uid, name, start, end, prefix, column):
        # uid for debugging
        self.uid = uid
        self.name = name
        self.start = start
        self.end = end
        self.prefix = prefix
        # to which column was the region set (without the prefix)
        self.column = column


# a custom error if something goes wrong
class MarsPyException(Exception):
    def __init___(self):
        Exception.__init__(self)


class MarsPyWarning(UserWarning):
    def __init__(self):
        UserWarning.__init__(self)
