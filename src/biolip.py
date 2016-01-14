#!/usr/bin/env python


import os
import pybel
import shlex
import subprocess32
import tempfile


class FastSearch:
    def __init__(self, lig_path,
                 index="/ddnB/work/jaydy/dat/BioLip/ligand_nr.fs",
                 tani=0.5, minimum_size=6):
        self.lig_path = lig_path
        self.index = index
        self.tani = tani
        self.minimum_size = minimum_size

    def search(self, verbose=True):
        try:
            ofn = tempfile.mkstemp(suffix='.sdf')[1]
            cmd = "babel %s %s -s %s -at%f" % (self.index,
                                               ofn,
                                               self.lig_path,
                                               self.tani)
            if verbose:
                print cmd
            args = shlex.split(cmd)
            subprocess32.call(args)
            mols = [mol for mol in pybel.readfile(format="sdf", filename=ofn)]
            return mols
        except Exception as detail:
            print "FAIL Fastsearch for %s" % (self.lig_path)
            print detail
        finally:
            os.remove(ofn)


class BioLipQuery(FastSearch):

    def correspondingPrtPdb(self, lig):
        self.prt_dir = "/work/jaydy/dat/BioLip/prt/"
        mid_two, code = lig.title.split('/')
        pdb = code.split('_')[0] + code.split('_')[2] + ".pdb"
        path = os.path.join(self.prt_dir, mid_two, pdb)
        assert os.path.exists(path)
        return path


def test():
    query = BioLipQuery("/work/jaydy/dat/BioLip/sml_ligand_nr/03/103m_NBN_A_1.pdb")
    mols = query.search()
    pdb_path = query.correspondingPrtPdb(mols[0])
    print pdb_path


if __name__ == '__main__':
    pass
