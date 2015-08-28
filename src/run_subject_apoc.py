#!/usr/bin/env python

import luigi
from run_control_apoc import LpcApocResultTask, LpcKcombuResult

if __name__ == '__main__':
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    if tname != "tname":
        luigi.build([LpcApocResultTask(tname, qname, "subject"),
                     LpcKcombuResult(qname, tname, "subject"),
                     LpcKcombuResult(tname, qname, "subject")],
                    local_scheduler=True)
