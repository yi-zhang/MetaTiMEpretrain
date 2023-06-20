import sys
import pandas as pd
mypyPath = '../scICA/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)
import loompy
import os
import config

loomfiles = [t for t in os.listdir( config.SCDIR )  if t[-5:]== '.loom' ]

import json

# assume you have the following dictionary


obsdict={}
for loomfile in loomfiles:
    name = loomfile.replace('.loom','')
    with loompy.connect( os.path.join( config.SCDIR, loomfile) ) as ds:
        obsdict.update({ name: ds.ca.keys()})
        
with open( os.path.join(config.SCMETADIR, 'scRNA_obs_columns.json'), "w") as write_file:
    json.dump(obsdict, write_file, indent=4)
print('Metadata extracted in ')
print(os.path.join(config.SCMETADIR, 'scRNA_obs_columns.json'))


# layers

obsdict={}
for loomfile in loomfiles:
    name = loomfile.replace('.loom','')
    with loompy.connect( os.path.join( config.SCDIR, loomfile) ) as ds:
        obsdict.update({ name: ds.layers.keys()})
        
with open( os.path.join(config.SCMETADIR, 'scRNA_layers.json'), "w") as write_file:
    json.dump(obsdict, write_file, indent=4)
print('Metadata extracted in ')
print(os.path.join(config.SCMETADIR, 'scRNA_layers.json'))
