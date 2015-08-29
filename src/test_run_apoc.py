import luigi
import unittest
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

    # def test_f_xcms(self):
    #     luigi.build([LpcApocXcms('3vn9_ANK_A_401', '4ej7_ATP_C_401'),
    #                  LpcApocXcms('3v76_FDA_A_547', '3zxs_FAD_A_1509'),
    #                  LpcApocXcms('4a2a_ATP_B_1391', '4a5a_ANP_A_700')],
    #                 local_scheduler=True)

if __name__ == "__main__":
    unittest.main()
