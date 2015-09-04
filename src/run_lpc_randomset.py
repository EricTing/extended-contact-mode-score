#!/usr/bin/env python

import luigi
from xcms import LpcApocXcms


def main(tname, qname):
    luigi.build([LpcApocXcms(tname, qname, "rs2")],
                local_scheduler=True)

if __name__ == '__main__':
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    main(tname, qname)
