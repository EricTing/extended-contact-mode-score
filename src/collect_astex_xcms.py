#!/usr/bin/env python

import pandas as pd
from astex_xcms import AstexXcmsBackCompatibility

complexes = [_.rstrip() for _ in file("../dat/astex_sdf.lst")]

all_dset = pd.DataFrame()
for sdf_id in complexes:
    path = AstexXcmsBackCompatibility(sdf_id).output().path
    dset = pd.read_csv(path, sep=' ')
    dset['complex'] = sdf_id
    all_dset = pd.concat([all_dset, dset], axis=0)

all_dset.to_csv("../dat/axtex_xcms.csv")
