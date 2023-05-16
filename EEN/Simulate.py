from __future__ import division

import random
import math


iterations = 1000
random.seed(10815657)
count = 0
ExonDict={}
for introns in range(1,20):
    for i in range(0, iterations):
        random_draws =[]
        for foo in range(0, introns):
            random_draws.append(random.uniform(0, 1))
        random_draws.append(0)
        random_draws.append(1)
        random_draws.sort()
        sumlens=0
        for index in range(0,len(random_draws)-1):
            exonlen=random_draws[index+1]-random_draws[index]
            sumlens=sumlens+exonlen**2
        itlen=1/sumlens
        if ExonDict.has_key(introns):
            ExonDict[introns].append(itlen)
        else:
            ExonDict[introns]=[itlen]

for item in ExonDict:
    mylens=ExonDict[item]
    print item, sum(mylens)/len(mylens)
    
                        
