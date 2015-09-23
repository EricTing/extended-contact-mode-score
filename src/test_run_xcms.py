import luigi
import unittest
import os
from run_control_apoc import LpcKcombuResult, LigandExpStructureInMol2


class Test(unittest.TestCase):

    def test_a_mol2(self):

        to_build = [LigandExpStructureInMol2('3vn9_ANK_A_401'),
                    LigandExpStructureInMol2('3v76_FDA_A_547'),
                    LigandExpStructureInMol2('4a2a_ATP_B_1391')]

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


if __name__ == "__main__":
    unittest.main()
