import luigi
import os
import unittest
from xcms import LpcApocXcms
from run_control_apoc import PkcombuAtomMatchParser, getPdbAtomsBySerialNum
from run_control_apoc import ApocResultParer, LpcApocResultTask
from run_control_apoc import LpcKcombuResult, LpcPocketPathTask
from apoc_inputs import DecompressedPdb


class Test(unittest.TestCase):

    def setUp(self):
        self.task = LpcApocResultTask("3v76_FDA_A_547", "3zxs_FAD_A_1509")
        luigi.build([self.task],
                    local_scheduler=True)

    def test_a_ApocResultPaser(self):
        with self.task.output().open("r") as f:
            parser = ApocResultParer(f.read())
            global_alignment = parser.queryGlobal("3v76A", "3zxsA")
            self.assertEqual(0.26646, global_alignment.tm_score)
            self.assertEqual(6.32, global_alignment.rmsd)
            self.assertEqual(0.070, global_alignment.seq_identity)

            pocket_alignment = parser.queryPocket("3v76_FDA_A_547", "3zxs_FAD_A_1509")
            self.assertEqual(0.34803, pocket_alignment.ps_score)
            self.assertEqual(3.25, pocket_alignment.rmsd)
            self.assertEqual('A', pocket_alignment.template_chainid)
            self.assertListEqual([402, 44], pocket_alignment.template_res)
            self.assertListEqual([269, 303], pocket_alignment.query_res)

            pocket_alignment = parser.queryPocket("3v76_FDA_A_547", "3zxs_SF4_A_1510")
            self.assertEqual(0.36987, pocket_alignment.ps_score)
            self.assertEqual(2.55, pocket_alignment.rmsd)
            if pocket_alignment.has_pocket_alignment:
                self.assertEqual('A', pocket_alignment.template_chainid)
                self.assertListEqual([402, 44], pocket_alignment.template_res)
                self.assertListEqual([269, 303], pocket_alignment.query_res)

    def test_c_kcombu(self):
        to_build = LpcKcombuResult('3vn9_ANK_A_401', '4ej7_ATP_C_401')
        try:
            os.remove(to_build.output().path)
        except:
            pass
        luigi.build([to_build],
                    local_scheduler=True)

        t_pdb = to_build.requires()[0].output().path
        q_pdb = to_build.requires()[1].output().path

        kcombu_parser = PkcombuAtomMatchParser(to_build.output().path)
        list_t, list_q = kcombu_parser.getMatchingSerialNums()

        t_atoms = getPdbAtomsBySerialNum(t_pdb, list_t)
        q_atoms = getPdbAtomsBySerialNum(q_pdb, list_q)

        t_eles = [atom.element for atom in t_atoms]
        q_eles = [atom.element for atom in q_atoms]

        self.assertListEqual(t_eles, q_eles)

    def test_d_pdb(self):
        luigi.build([DecompressedPdb('104m')],
                    local_scheduler=True)

    def test_e_pocket(self):
        luigi.build([LpcPocketPathTask('3vn9_ANK_A_401')],
                    local_scheduler=True)
        luigi.build([LpcPocketPathTask('4ej7_ATP_C_401')],
                    local_scheduler=True)


if __name__ == "__main__":
    unittest.main()
