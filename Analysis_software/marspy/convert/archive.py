import pandas as pd
from jnius import autoclass

from awesome_data import DataSet
from marspy.convert.molecule import *


class Archive:

    def __init__(self, filepath):
        self.filepath = filepath
        self.name = self.filepath.split('/')[-1]
        self.File = autoclass('java.io.File')
        self.yamaFile = self.File(self.filepath)

    def get_molecule_by_uid(self, uid):
        raise NotImplementedError

    def get_molecules_by_tags(self, tags):
        raise NotImplementedError

    def validate_params(self):
        pass


class SingleMoleculeArchive(Archive):
    instances = []

    def __init__(self, filepath, accept_tag, label=dict()):
        Archive.__init__(self, filepath)
        self.instances.append(self)
        self.Archive = autoclass('de.mpg.biochem.mars.molecule.SingleMoleculeArchive')
        self.archive_link = self.Archive(self.yamaFile)
        self.metadata_uids = tuple(sc.to_python(self.archive_link.getMetadataUIDs()))
        self.label = label

        # nucleotide
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nucleotide')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nucleotide = 'n/a'
            print(f'nucleotide not found. Setting default to {self.nucleotide}')
        # parameter properly set
        else:
            self.nucleotide = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')

        # NaCl concentration
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nacl')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nacl = 'n/a'
            print(f'NaCl concentration not found. Setting default to {self.nacl}')
        # parameter properly set
        else:
            self.nacl = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nacl')

        # MCM variant
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('mcm')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('mcm')) == 0:
            # default n/a
            self.mcm = 'n/a'
            print(f'MCM variant not found. Setting default to {self.mcm}')
        # parameter properly set
        else:
            self.mcm = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('mcm')

        self.protein = list(self.label.keys())[0]

        # instantiate a new SingleMolecule for each uid and store instances as list
        self.molecules = [SingleMolecule(uid, self.protein, archive=self.archive_link) for uid in
                          sc.to_python(self.archive_link.getMoleculeUIDs()) if
                          self.archive_link.get(uid).hasTag(accept_tag)]
        self.tags = set()
        for molecule in self.molecules:
            self.tags.update(molecule.tags)

    def get_molecule_by_uid(self, uid):
        """
        Returns molecule object with provided UID.
        """
        return list(filter(lambda molecule: molecule.uid == uid, self.molecules))[0]

    def get_molecules_by_tags(self, tags):
        """
        Provide tags as list.
        Returns filter of all molecules which have all the specified tags
        """
        return filter(lambda molecule: set(tags).issubset(set(molecule.tags)), self.molecules)

    def __len__(self):
        return len(self.molecules)


class DnaMoleculeArchive(Archive):
    instances = []

    def __init__(self, filepath, accept_tag, labels=dict()):
        Archive.__init__(self, filepath)
        self.instances.append(self)
        self.Archive = autoclass('de.mpg.biochem.mars.molecule.DnaMoleculeArchive')
        self.archive_link = self.Archive(self.yamaFile)
        self.metadata_uids = tuple(sc.to_python(self.archive_link.getMetadataUIDs()))
        self.dna_molecule_count = 0
        for metadata in self.metadata_uids:
            self.dna_molecule_count += dict(sc.to_python(self.archive_link.getMetadata(metadata).getParameters()))[
                'DnaMoleculeCount']
        # subtract # of reject_dna tags
        self.dna_molecule_count -= len(list(filter(lambda uid:
                                                   self.archive_link.get(uid).hasTag('reject_dna'),
                                                   sc.to_python(self.archive_link.moleculeUIDs))))
        self.labels = labels

        # nucleotide
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nucleotide')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nucleotide = 'n/a'
            print(f'nucleotide not found. Setting default to {self.nucleotide}')
        # parameter properly set
        else:
            self.nucleotide = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')

        # NaCl concentration
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('nacl')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nucleotide')) == 0:
            # default n/a
            self.nacl = 'n/a'
            print(f'NaCl concentration not found. Setting default to {self.nacl}')
        # parameter properly set
        else:
            self.nacl = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('nacl')

        # MCM variant
        # check if all metadata parameters match & raise warning if conditions to in one archive are not identical
        if len({self.archive_link.getMetadata(metadata_uid).getStringParameter('mcm')
                for metadata_uid in self.metadata_uids}) > 1:
            raise MarsPyWarning()
        # if StringParameter is not set, getStringParameter returns empty string ''
        if len(self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('mcm')) == 0:
            # default n/a
            self.mcm = 'n/a'
            print(f'MCM variant not found. Setting default to {self.mcm}')
        # parameter properly set
        else:
            self.mcm = self.archive_link.getMetadata(self.metadata_uids[0]).getStringParameter('mcm')

        self.proteins = set()

        # will get all columns in DataTable with 'Protein_n_Position_on_Dna'
        for match in re.findall('\w+_Position_on_DNA', '$'.join(set(sc.to_python(
                self.archive_link.properties().getColumnSet())))):
            self.proteins.add(match.split('_')[0])

        # instantiate a new DnaMolecule for each uid and store instances as list
        self.molecules = [DnaMolecule(uid, self.proteins, archive=self.archive_link) for uid in
                          sc.to_python(self.archive_link.getMoleculeUIDs())
                          if self.archive_link.get(uid).hasTag(accept_tag)]

        # define archive tags union of all molecule tags
        # define archive prefixes as union of all molecule prefixes (will be used for top level columns in big df later)
        self.tags = set()
        self.prefixes = set()
        for molecule in self.molecules:
            self.tags.update(molecule.tags)
            self.prefixes.update(molecule.prefixes)

    def validate_params(self):
        """
        Integrity check of passed Archive.
        """
        # compare number protein in params vs actual one (retrieved from metadata)
        for molecule in self.molecules:
            # take global protein to confirm dict was pasted correctly
            for protein in self.proteins:
                if not (molecule.proteins[protein] == molecule.params['Number_' + protein]):
                    err_message = f"Conflict in molecule {molecule.uid}!\n\
                    Number of {protein} retrieved from metadata: {molecule.proteins[protein]}\n\
                    Number of {protein} based on Parameter: {molecule.params['Number_' + protein]}"
                    raise MarsPyException(err_message)

        return 'passed'

    def add_df_noidle(self, prefix, specifier):
        """
        Generates a copy of molecule.df (df_noidle) with all rows removed falling in specified regions
        """
        for molecule in self.molecules:

            molecule.df_noidle = molecule.df.copy()
            # list of rows marked for removal
            drop_rows = []

            for region in molecule.regions:
                if specifier in region.name and region.prefix == prefix:
                    for row, index in molecule.df.iterrows():
                        if molecule.df.loc[row, prefix + region.column] in range(int(region.start),
                                                                                 int(region.end) + 1):
                            drop_rows.append(row)

            molecule.df_noidle.drop(drop_rows, inplace=True)

    def get_molecule_by_uid(self, uid):
        """
        Returns molecule object with provided UID.
        """
        return list(filter(lambda molecule: molecule.uid == uid, self.molecules))[0]

    def get_molecules_by_tags(self, tags):
        """
        Provide tags as list.
        Returns filter of all molecules which have all the specified tags
        """
        return filter(lambda molecule: set(tags).issubset(set(molecule.tags)), self.molecules)

    def __len__(self):
        return len(self.molecules)


def instantiate_archive(name, datasets):
    """
    Instantiates passed archive from underlying dataset
    """
    # check if we have the right data type
    for data in datasets:
        if not isinstance(data, DataSet):
            raise MarsPyException('Dataset contains non-compatible data type.')
    data = list(filter(lambda dataset: dataset.name == name, datasets))[0]
    if data.archive_type == 'DnaMoleculeArchive':
        DnaMoleculeArchive(filepath=data.filepath + data.name, accept_tag=data.accept_tag, labels=data.labels)
    elif data.archive_type == 'SingleMoleculeArchive':
        SingleMoleculeArchive(filepath=data.filepath + data.name, accept_tag=data.accept_tag, label=data.labels)
    else:
        raise MarsPyException(f'Failed to instantiate Archive {data.name}.')


def describe_archives(archives):
    """
    Describes passed archives by returning a pandsa DataFrame. Pass archives as iterable object
    """
    df = pd.DataFrame(
        columns=['# of datasets', '# of molecules', 'labeled proteins', 'nucleotide', 'NaCl concentration',
                 'MCM variant', 'archive validation'])
    for archive in archives:
        _temp_df = pd.DataFrame(index=[archive.name.split('.')[0]],
                                data=[[len(archive.metadata_uids), len(archive),
                                       '; '.join([label + '-' + protein for protein, label in archive.labels.items()]),
                                       archive.nucleotide, archive.nacl, archive.mcm, archive.validate_params()]],
                                columns=['# of datasets', '# of molecules', 'labeled proteins', 'nucleotide',
                                         'NaCl concentration', 'MCM variant', 'archive validation'])
        df = pd.concat([df, _temp_df])
    df = df.infer_objects()
    return df
