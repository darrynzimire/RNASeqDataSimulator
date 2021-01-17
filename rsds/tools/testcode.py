# encoding=utf-8

import pickle as pickle
import gzip
import pyfaidx
import numpy as np
from rsds.process_inputFiles import proc_tx_expmodel as ptx

readtot = 500000
f = 'test_subset_p.gz'
model = ptx(f)
id = model[2]
# print(id)
counts_p = np.rint(np.multiply(id, readtot)).astype(int)
print(sum(counts_p))
