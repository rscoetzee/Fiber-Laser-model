# -*- coding: utf-8 -*-
"""
Created on Tue May 24 11:51:27 2016

@author: Riaan
"""

import numpy as np

z = np.random.rand(1000,2000)


from datetime import datetime
startTime = datetime.now()

for i in range(0,len(z)):
    X = np.sum(z,axis=0)
    
print(datetime.now() - startTime)

