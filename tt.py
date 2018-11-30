import pandas as pd
import numpy as np

df=pd.DataFrame(np.random.rand(8,4),columns=['1s','1d','2s','2d'])

g=df.columns.map(lambda x:x[0])
df.groupby(g,axis=1).sum()


