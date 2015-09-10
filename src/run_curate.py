#!/usr/bin/env python

import luigi
from curate import Curate


def main(tname, subset):
    luigi.build([Curate.PairWisePsScore(tname=tname,
                                        subset=subset),
                 Curate.PairWiseTanimoto(tname=tname,
                                         subset=subset)],
                local_scheduler=True)
    pass

if __name__ == '__main__':
    import sys
    main(sys.argv[1], "subject")
