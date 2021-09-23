import numpy as np
import pandas as pd
import os
import sys
import ascatparser


args = ascatparser.parse_args(sys.argv[1:])

allelic_capseg, gistic = ascatparser.helper(args)

allelic_capseg.to_csv(args.allelic_output, sep='\t', index=False)

gistic.to_csv(args.gistic_output, sep='\t', index=False)

# This is the command I want to use to be able to run this script as intended:
# python /src/ascatparser.py ${caveman_tumor} ${copynumber_tumor} ${copynumber_normal} ${ID}.ascat_allelic_capseg.tsv ${ID}.ascat.seg 0
# python ./ascatparser.py test_files/input_files/CTSP-ACYS-TTP1-E-1-1-D-A791-36.copynumber.caveman.csv test_files/input_files/CTSP-ACYS-TTP1-E-1-1-D-A791-36.copynumber.txt.gz test_files/input_files/CTSP-ACYS-NB1-A-1-0-D-A791-36.count test_files/output_files/test_ACYS_ascat_allelic_capseg.tsv test_files/output_files/test_ACYS.ascat.seg

