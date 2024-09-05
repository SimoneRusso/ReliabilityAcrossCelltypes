from allensdk.core.cell_types_cache import CellTypesCache
import os

os.chdir(r'C:\Users\simon\PycharmProjects\pythonProject1\allensdk\CellType\Data')
os.chdir(r'D:\AllenInstitute\CellType')

ctc = CellTypesCache(manifest_file='cell_types/manifest.json')

# a list of cell metadata for cells with reconstructions, download if necessary
cells = ctc.get_cells(require_reconstruction=False)

# open the electrophysiology data of one cell, download if necessary
for counterCell in range(len(cells)):
    data_set = ctc.get_ephys_data(cells[counterCell]['id'])

# read the reconstruction, download if necessary
#reconstruction = ctc.get_reconstruction(cells[0]['id'])