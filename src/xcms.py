#!/usr/bin/env python

from sklearn.metrics import matthews_corrcoef  # matthews_corrcoef has to be imported before everything to avoid segfault
import luigi
import os
import pybel
from scipy.spatial.distance import euclidean
from Bio.PDB import PDBParser

from urls import APOC_WORKING_DIR
from run_control_apoc import ApocResultParer, LpcPocketPathTask
from run_control_apoc import LpcKcombuResult, LpcApocResultTask


def buildArrayOfContact(residue_coords, atom_coords):
    return [1 if euclidean(res, atom) < 4.2 else 0
            for res in residue_coords
            for atom in atom_coords]


class LpcApocXcms(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def _convert(self, name):
        return name.split('_')

    def requires(self):
        return [LpcKcombuResult(self.tname,
                                self.qname,
                                self.subset),
                LpcKcombuResult(self.qname,
                                self.tname,
                                self.subset),
                LpcApocResultTask(self.tname,
                                  self.qname,
                                  self.subset),
                LpcPocketPathTask(self.tname),
                LpcPocketPathTask(self.qname)]

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
                                self.my_midtwo())
        return dir_path

    def output(self):
        path = os.path.join(self.mydir(),
                            "%s__%s.xcms" % (self.tname, self.qname))
        return luigi.LocalTarget(path)

    def _select_residues_coords(self, name, chain_id, selected):
        pdb_path = LpcPocketPathTask(name).output().path
        pdb_parser = PDBParser(QUIET=True)
        pdb_structure = pdb_parser.get_structure(name, pdb_path)
        chain = pdb_structure[0][chain_id]
        coords = []
        for res_id in selected:
            for atom in chain[res_id].get_unpacked_list():
                coords.append(atom.get_coord())

        return coords

    def _select_ligand_atom_coords(self, tname, qname):
        sdf_path = LpcKcombuResult(tname, qname, self.subset).output().path
        ligand = pybel.readfile("sdf", sdf_path).next()
        ligand.removeh()
        coords = [atom.coords for atom in ligand.atoms]
        return coords

    def run(self):
        with LpcApocResultTask(self.tname,
                               self.qname,
                               self.subset).output().open('r') as f:
            apoc_parser = ApocResultParer(f.read())

        f = self.output().open('w')
        f.write("Inputs\n")
        f.write("================================================================================\n")
        f.write("template:\t\t%s\nquery:\t\t%s\n" % (self.tname, self.qname))
        f.write("Pocket for template:\t\t%s\n" % LpcPocketPathTask(self.tname).output().path)
        f.write("Pocket for query:\t\t%s\n" % LpcPocketPathTask(self.qname).output().path)
        f.write("Apoc result:\t\t%s\n" % LpcApocResultTask(self.tname, self.qname, self.subset).output().path)
        f.write("MCS for template in sdf format:\t\t%s\n" % LpcKcombuResult(self.tname, self.qname, self.subset).output().path)
        f.write("MCS for query in sdf format:\t\t%s\n" % LpcKcombuResult(self.qname, self.tname, self.subset).output().path)

        pocket_alignment = apoc_parser.queryPocket(self.tname, self.qname)
        if pocket_alignment.has_pocket_alignment:
            f.write("Calculating CMS\n")
            f.write("================================================================================\n")
            t_res = self._select_residues_coords(self.tname,
                                                 pocket_alignment.template_chainid,
                                                 pocket_alignment.template_res)
            t_lig = self._select_ligand_atom_coords(self.tname,
                                                    self.qname)

            t_contact = buildArrayOfContact(t_res, t_lig)
            f.write("Contact in template's complex\n")
            f.write(" ".join(map(str, t_contact)) + "\n")

            # the second complex
            q_res = self._select_residues_coords(self.qname,
                                                 pocket_alignment.query_chainid,
                                                 pocket_alignment.query_res)

            q_lig = self._select_ligand_atom_coords(self.qname,
                                                    self.tname)

            q_contact = buildArrayOfContact(q_res, q_lig)
            f.write("Contact in query's complex\n")
            f.write(" ".join(map(str, q_contact)) + "\n")

            cms = matthews_corrcoef(t_contact, q_contact)
            f.write("Contact Mode Score:\t\t%f\n" % cms)

        f.close()
        print "xcms output %s" % (self.output().path)


def main():
    pass

if __name__ == '__main__':
    main()
