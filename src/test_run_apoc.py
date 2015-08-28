import luigi
import unittest
from run_control_apoc import ApocResultParer, LpcApocResultTask
from run_control_apoc import LpcKcombuResult, LpcPocketPathTask
from apoc_inputs import DecompressedPdb
from my_pdb import buildSelect
from xcms import LpcApocXcms
from Bio.PDB import PDBParser, PDBIO


class Test(unittest.TestCase):

    def setUp(self):
        self.task = LpcApocResultTask("3sis_MN0_A_6535", "10gs_VWW_A_210")
        luigi.build([self.task],
                    local_scheduler=True)

    def test_a_ApocResultPaser(self):
        with self.task.output().open("r") as f:
            parser = ApocResultParer(f.read())
            self.assertEqual(0.25506, parser.global_property.tm_score)
            self.assertEqual(5.39, parser.global_property.rmsd)
            self.assertEqual(0.034, parser.global_property.seq_identity)

            self.assertEqual(1.87, parser.pocket_property.rmsd)
            self.assertEqual(0.125, parser.pocket_property.seq_identity)
            self.assertEqual(0.29410, parser.pocket_property.ps_score)
            self.assertSequenceEqual([["A", 187], ["A", 50]],
                                     parser.matching_list[0])

    def test_b_buildSelect(self):
        pdb_parser = PDBParser(QUIET=True)

        t_ifn = "../dat/pdb3sis.ent"
        t_structure = pdb_parser.get_structure('t', t_ifn)

        io = PDBIO()
        io.set_structure(t_structure)

        f = self.task.output().open('r')
        parser = ApocResultParer(f.read())
        selected_residues = [_[0] for _ in parser.matching_list]
        selected = buildSelect(selected_residues)
        io.save("test.pdb", selected())

        f.close()

    def test_c_kcombu(self):
        luigi.build([LpcKcombuResult('104m_NBN_A_156', '1m8e_7NI_A_906')],
                    local_scheduler=True)

    def test_d_pdb(self):
        luigi.build([DecompressedPdb('104m')],
                    local_scheduler=True)

    def test_e_pocket(self):
        luigi.build([LpcPocketPathTask('3vn9_ANK_A_401')],
                    local_scheduler=True)
        luigi.build([LpcPocketPathTask('4ej7_ATP_C_401')],
                    local_scheduler=True)

    def test_f_pocket(self):
        luigi.build([LpcApocXcms('3vn9_ANK_A_401', '4ej7_ATP_C_401')],
                    local_scheduler=True)

if __name__ == "__main__":
    unittest.main()
