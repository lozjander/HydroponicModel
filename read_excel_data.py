import pandas as pd

fpath = "D:\PhD Backup - 04.12.2018\Python Model\Tomato Model\Data\heuvelink (1995) Detailed.xlsx"
#df = pd.read_excel(fpath, header=[0,1], index_col=[0,1])
df = pd.read_excel(fpath)
d = df['1 - Crosses'].values[1:]
