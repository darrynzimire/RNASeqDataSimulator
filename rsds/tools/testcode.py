# encoding=utf-8

import pickle as pickle
import gzip
import pyfaidx


f = gzip.open('chr1FLmodel_p.gzip', 'rb')
model = pickle.load(f)
distparam= model[0]

aic = model[1]
bic = model[2]

print(model)