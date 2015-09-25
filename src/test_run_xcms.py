import luigi
import subprocess32
import unittest
import os
from run_control_apoc import LpcKcombuResult, LigandExpStructureInMol2
from run_control_apoc import PkcombuAtomMatchParser
from run_control_apoc import LpcApocResultTask, ApocResultParer
from run_control_apoc import LigandMatchingList, ProteinMatchingList
from apoc_inputs import DecompressedPdb


class Test(unittest.TestCase):

    def test_a_mol2(self):

        to_build = [LigandExpStructureInMol2('3vn9_ANK_A_401'),
                    LigandExpStructureInMol2('3v76_FDA_A_547'),
                    LigandExpStructureInMol2('4a2a_ATP_B_1391'),
                    LigandExpStructureInMol2('4ej7_ATP_C_401')]

        for task in to_build:
            try:
                os.remove(task.output().path)
            except:
                pass

        luigi.build(to_build, local_scheduler=True)

    def test_b_kcombu(self):
        to_build = [LpcKcombuResult('3vn9_ANK_A_401', '4ej7_ATP_C_401'),
                    LpcKcombuResult('3v76_FDA_A_547', '3zxs_FAD_A_1509'),
                    LpcKcombuResult('4a2a_ATP_B_1391', '4a5a_ANP_A_700')]

        for task in to_build:
            try:
                os.remove(task.output().path)
            except:
                pass

        luigi.build(to_build, local_scheduler=True)

    def test_c_kcombu_parser(self):
        oam_fn = LpcKcombuResult('3vn9_ANK_A_401', '4ej7_ATP_C_401').output().path
        parser = PkcombuAtomMatchParser(oam_fn)
        list_a, list_b = parser.getMatchingSerialNums()
        print len(list_a), len(list_b)
        self.assertEqual(list_a[1], 2)
        self.assertEqual(list_b[5], 26)
        self.assertEqual(list_b[25], 13)
        parser.writeMatchingSerialNums("./3vn9_ANK_A_401_4ej7_ATP_C_401.lml")

    def test_d_apoc_parser(self):
        t_name = "3vn9_ANK_A_401"
        q_name = "4ej7_ATP_C_401"

        to_build = [LpcApocResultTask(t_name, q_name)]
        luigi.build(to_build, local_scheduler=True)

        apoc_out = LpcApocResultTask(t_name, q_name).output().path
        with open(apoc_out, 'r') as ifs:
            parser = ApocResultParer(ifs.read())
        pocket_alignment = parser.queryPocket(t_name, q_name)
        self.assertEqual(pocket_alignment.template_res[5], 65)
        self.assertEqual(pocket_alignment.query_res[20], 203)

        parser.writeProteinMatchingList(t_name, q_name,
                                        "./%s_%s.pml" % (t_name, q_name))

    def test_e_matchinglist(self):
        t_name = "3vn9_ANK_A_401"
        q_name = "4ej7_ATP_C_401"
        subset = "subject"

        to_build = [LigandMatchingList(t_name, q_name, subset),
                    ProteinMatchingList(t_name, q_name, subset)]

        luigi.build(to_build, local_scheduler=True)

    def test_f_xcms(self):
        t_name = "3vn9_ANK_A_401"
        q_name = "4ej7_ATP_C_401"
        subset = "subject"

        t_lig_path = LigandExpStructureInMol2(t_name).output().path
        q_lig_path = LigandExpStructureInMol2(q_name).output().path

        t_pdb = t_name.split('_')[0]
        q_pdb = q_name.split('_')[0]

        luigi.build([DecompressedPdb(t_pdb),
                     DecompressedPdb(q_pdb)],
                    local_scheduler=True)

        t_prt_path = DecompressedPdb(t_pdb).output().path
        q_prt_path = DecompressedPdb(q_pdb).output().path

        lml_path = "./%s_%s.lml" % (t_name, q_name)
        pml_path = "./%s_%s.pml" % (t_name, q_name)

        cmds = ["run_xcms",
                "--la", t_lig_path,
                "--lb", q_lig_path,
                "--pa", t_prt_path,
                "--pb", q_prt_path,
                "--lml", lml_path,
                "--pml", pml_path]

        try:
            subprocess32.call(cmds)
        except:
            pass


if __name__ == "__main__":
    unittest.main()
