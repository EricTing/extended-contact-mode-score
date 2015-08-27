import luigi
import os
import unittest
from run_apoc import ApocResultParer, LpcApocResultTask


class Test(unittest.TestCase):

    def setUp(self):
        self.task = LpcApocResultTask("3sis_MN0_A_6535", "10gs_VWW_A_210")
        print self.task.output().path
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
            self.assertSequenceEqual(("A 187 G", "A 50 G"),
                                     parser.matching_list[0])


if __name__ == "__main__":
    unittest.main()
