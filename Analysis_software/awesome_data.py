from collections import namedtuple

# before you start: change filepath to archives location
filepath = '/Volumes/pool-duderstadt/Matthias/git/PhD/data/MCM-cohesin barrier/Final_archives/'

# do not change the code below
DataSet = namedtuple('DataSet', ['filepath', 'name', 'archive_type', 'accept_tag', 'labels'])

datasets = [DataSet(filepath=filepath,
                    name='Cohesin-MCM_wt_highsalt.yama', archive_type='DnaMoleculeArchive',
                    accept_tag='accept', labels=dict(Cohesin='Halo-JF549', MCM='ybbR-LD655')),
            DataSet(filepath=filepath,
                    name='Cohesin-MCM_wt_lowsalt.yama', archive_type='DnaMoleculeArchive',
                    accept_tag='accept', labels=dict(Cohesin='Halo-JF549', MCM='ybbR-LD655')),
            DataSet(filepath=filepath,
                    name='Cohesin-only_lowsalt.yama', archive_type='DnaMoleculeArchive',
                    accept_tag='accept', labels=dict(Cohesin='Halo-JF549'))
            ]
