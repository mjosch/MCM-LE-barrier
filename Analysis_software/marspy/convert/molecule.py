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


class SegmentsTable:
    """
    SegmentsTable object holding the actual df, with additional information in attributes.
    Also contains specific methods for filtering, bleaching steps, pause detection.
    """

    def __init__(self, molecule, prefix, col_x, col_y, region):
        # uid for debugging
        self.uid = molecule.uid
        # which prefix does it belong to
        self.prefix = prefix
        # remove prefix for x and y columns => general naming
        self.col_x = col_x.split(self.prefix)[-1]
        self.col_y = col_y.split(self.prefix)[-1]
        self.region = region
        # actual SegmentsTable()
        self.df = sc.to_python(molecule.archive.get(molecule.uid).getSegmentsTable(col_x, col_y, region))
        # type of SegmentTable (default None)
        self.type = None
        # keep track if seg_dfs were already filtered before data is interpreted
        self.filtered = False
        # assign type of SegmentTable (one needs to be True)
        if self.col_y == 'Intensity':
            self.type = 'bleaching'
        elif self.col_y == 'Position_on_DNA':
            self.type = 'rate'
        else:
            err_message = f"Conflict in molecule {self.uid}!\nSegmentTable type could not be assigned!"
            raise MarsPyException(err_message)

    def filter_segments(self, sigma_b_max=0):
        """
        Mode 2: SegmentsTable type: 'rate' -
        Reject all segments with B value (velocity) < b_min and sigma_B > sigma_b_max (poor fits)
        """

        # Mode 1 - type 'rate'
        if self.type == 'rate':

            # need to remove rows at the end
            remove_rows = set()
            # loop through rows and find segments which match exclusion criteria
            for i in range(len(self.df)):
                if self.df.loc[i, 'Sigma_B'] >= sigma_b_max:
                    remove_rows.add(i)

            # remove all rows & update indices
            self.df.drop(list(remove_rows), axis=0, inplace=True)
            self.df.reset_index(drop=True, inplace=True)

            # successfully filtered
            self.filtered = True

    def detect_pauses(self, thresh=200, global_thresh=True, length=1, col='B'):
        """
        Detection pauses in SegmentTable (only for type = 'rate')
        global_thresh: Set to True if a fixed threshold for all molecules should be used
        thresh: threshold to detect pauses.
        length: minimal pause duration (s)
        If global_thresh is False, a molecule-specific threshold is calculated with thresh^-1 * np.mean(col)
        col: column evaluated for pauses
        """
        # only SegmentsTables with type 'rate'
        if self.type == 'rate':
            self.df['pause_' + col] = False
            cutoff = thresh
            # iterate once for each entry in seg_df
            for i in range(len(self.df)):
                if not global_thresh:
                    # redefine cutoff
                    cutoff = self.df[~self.df['pause_' + col]][col].mean() / thresh

                for row in self.df.index:
                    self.df.loc[row, 'pause_' + col] = (abs(self.df.loc[row, col]) < cutoff) and \
                                                       (self.df.loc[row, 'X2'] - self.df.loc[row, 'X1'] >= length)
            # if two subsequent segments are pauses merge them
            remove_rows = set()
            for i in range(1, len(self.df)):
                # both pauses
                # additional requirements: time values match and y values within 1 kb
                if (self.df.loc[i - 1, 'pause_B'] and self.df.loc[i, 'pause_B'] and
                        self.df.loc[i - 1, 'X2'] == self.df.loc[i, 'X1'] and
                        abs(self.df.loc[i - 1, 'Y2'] - self.df.loc[i, 'Y1'] < 1000)):
                    # recalculate values (x2 and y2 values in row i stay the same)
                    x1 = self.df.loc[i - 1, 'X1']
                    y1 = self.df.loc[i - 1, 'Y1']
                    a = np.average(self.df[i - 1:i + 1]['A'],
                                   weights=self.df[i - 1:i + 1]['X2'] - self.df[i - 1:i + 1]['X1'])
                    sigma_a = np.average(self.df[i - 1:i + 1]['Sigma_A'],
                                         weights=self.df[i - 1:i + 1]['X2'] - self.df[i - 1:i + 1]['X1'])
                    b = np.average(self.df[i - 1:i + 1]['B'],
                                   weights=self.df[i - 1:i + 1]['X2'] - self.df[i - 1:i + 1]['X1'])
                    sigma_b = np.average(self.df[i - 1:i + 1]['Sigma_B'],
                                         weights=self.df[i - 1:i + 1]['X2'] - self.df[i - 1:i + 1]['X1'])

                    # add row i-1 for removal and reassign calculated values for row i
                    remove_rows.add(i - 1)
                    self.df.loc[i, 'X1'] = x1
                    self.df.loc[i, 'Y1'] = y1
                    self.df.loc[i, 'A'] = a
                    self.df.loc[i, 'Sigma_A'] = sigma_a
                    self.df.loc[i, 'B'] = b
                    self.df.loc[i, 'Sigma_B'] = sigma_b

            # remove all rows & update indices
            self.df.drop(list(remove_rows), axis=0, inplace=True)
            self.df.reset_index(drop=True, inplace=True)


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
