import luigi
import unittest
import os
from xcms import LpcApocXcms


class Test(unittest.TestCase):

    def test_a_xcms(self):
        to_build = [LpcApocXcms('3vn9_ANK_A_401', '4ej7_ATP_C_401'),
                    LpcApocXcms('3v76_FDA_A_547', '3zxs_FAD_A_1509'),
                    LpcApocXcms('4a2a_ATP_B_1391', '4a5a_ANP_A_700')]

        for task in to_build:
            try:
                os.remove(task.output().path)
            except:
                pass

        luigi.build(to_build, local_scheduler=True)


if __name__ == "__main__":
    unittest.main()
