#!/usr/bin/env python

# matthews_corrcoef has to be imported before everything to avoid segfault
from sklearn.metrics import matthews_corrcoef
import luigi
import os
from scipy.spatial.distance import euclidean
from Bio.PDB import PDBParser

from urls import APOC_WORKING_DIR
from apoc_inputs import DecompressedPdb
from run_control_apoc import ApocResultParer, LpcPocketPathTask
from run_control_apoc import LpcKcombuResult, LpcApocResultTask
from run_control_apoc import PkcombuAtomMatchParser, getPdbAtomsBySerialNum


def buildArrayOfContact(residue_coords, atom_coords):
    return [1 if euclidean(res, atom) < 4.5 else 0
            for res in residue_coords
            for atom in atom_coords]


class LpcApocXcms(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def _convert(self, name):
        return name.split('_')

    def requires(self):
        t_pdb = self.tname.split('_')[0]
        q_pdb = self.qname.split('_')[0]
        return [LpcKcombuResult(self.tname,
                                self.qname,
                                self.subset),
                LpcApocResultTask(self.tname,
                                  self.qname,
                                  self.subset),
                DecompressedPdb(t_pdb),
                DecompressedPdb(q_pdb),
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
                            "%s__%s_atomic.xcms" % (self.tname, self.qname))
        return luigi.LocalTarget(path)

    def _kcombu_results(self):
        path = LpcKcombuResult(self.tname,
                               self.qname,
                               self.subset).output().path
        return PkcombuAtomMatchParser(path)

    def _select_residues(self, name, chain_id, selected):
        pdb_id = name.split('_')[0]
        pdb_path = DecompressedPdb(pdb_id).output().path
        pdb_parser = PDBParser(QUIET=True)
        pdb_structure = pdb_parser.get_structure(name, pdb_path)
        chain = pdb_structure[0][chain_id]
        coords = []
        names = []
        for res_id in selected:
            for atom in chain[res_id].get_unpacked_list():
                names.append(atom.get_name())
                coords.append(atom.get_coord())
        return coords, names

    def _select_ligand_atoms(self, tname, qname):
        kcombu_task = LpcKcombuResult(tname, qname)

        t_pdb = kcombu_task.requires()[0].output().path
        q_pdb = kcombu_task.requires()[1].output().path

        kcombu_parser = PkcombuAtomMatchParser(kcombu_task.output().path)
        list_t, list_q = kcombu_parser.getMatchingSerialNums()

        t_atoms = getPdbAtomsBySerialNum(t_pdb, list_t)
        q_atoms = getPdbAtomsBySerialNum(q_pdb, list_q)

        t_eles = [atom.element for atom in t_atoms]
        q_eles = [atom.element for atom in q_atoms]

        #  to guarantee atoms are of the same type
        assert(t_eles == q_eles)
        return t_atoms, q_atoms

    def _select_ligand_atom_coords(self):
        t_atoms, q_atoms = self._select_ligand_atoms(self.tname, self.qname)
        t_coords = [atom.coord for atom in t_atoms]
        q_coords = [atom.coord for atom in q_atoms]
        return t_coords, q_coords

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
        f.write("Kcombu matching atom result:\t\t%s\n" % LpcKcombuResult(self.tname, self.qname, self.subset).output().path)

        kcombu_data = self._kcombu_results().data
        f.write("Kcombu tanimoto:\t\t%f\n" % kcombu_data.tanimoto)

        t_coords, q_coords = self._select_ligand_atom_coords()

        pocket_alignment = apoc_parser.queryPocket(self.tname, self.qname)
        if pocket_alignment.has_pocket_alignment:
            f.write("Total number of matching residues:\t\t%d\n"
                    % (len(pocket_alignment.template_res)))
            f.write("Total number of matching ligand coordinates:\t\t%d\n"
                    % len(t_coords))
            f.write("PS-score:\t\t%f\n" % pocket_alignment.ps_score)
            f.write("p-value:\t\t%f\n" % pocket_alignment.p_value)
            f.write("z-score:\t\t%f\n" % pocket_alignment.z_score)

            t_prt_coords, t_prt_names = self._select_residues(self.tname,
                                                              pocket_alignment.template_chainid,
                                                              pocket_alignment.template_res)

            q_prt_coords, q_prt_names = self._select_residues(self.qname,
                                                              pocket_alignment.query_chainid,
                                                              pocket_alignment.query_res)

            assert(t_prt_names == q_prt_names)

            f.write("Total number of matching protein coordinates:\t\t%d\n"
                    % len(t_prt_coords))

            f.write("Calculating CMS\n")
            f.write("================================================================================\n")

            t_contact = buildArrayOfContact(t_prt_coords, t_coords)
            f.write("Contact in template's complex\n")
            f.write(" ".join(map(str, t_contact)) + "\n")


            q_contact = buildArrayOfContact(q_prt_coords, q_coords)
            f.write("Contact in query's complex\n")
            f.write(" ".join(map(str, q_contact)) + "\n")

            cms = matthews_corrcoef(t_contact, q_contact)
            f.write("Contact Mode Score:\t\t%f\n" % cms)

        f.close()
        print "xcms output %s" % (self.output().path)


def main(tname, qname):
    if tname != "tname":
        luigi.build([LpcApocXcms(tname, qname, "subject")],
                    local_scheduler=True)


if __name__ == '__main__':
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    main(tname, qname)
