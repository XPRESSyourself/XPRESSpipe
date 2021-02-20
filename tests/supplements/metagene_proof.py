import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline

"""
- Meta-coordinate is taken as only the mid-point of the read to avoid any biases from longer reads counting more to meta-profile than shorter reads
- Normalize each transcript so that super-expressors contribute as much to the meta-profile as lower expressors
    - If this doesn't happen, a transcript with 1000 reads will be completely swamped out by one that has 500000 reads
"""

# Example where super-expressor would normally swamp out other reads where the over trend is really a flat line
fig = plt.figure()
ax = fig.add_subplot(111)

df = pd.DataFrame()
df['transcript'] = ['t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t2','t2','t2','t2','t3','t3','t3']
df['meta'] = [1,1,1,1,1,1,2,2,2,3,3,3,1,2,3,3,1,2,3]

df_t = pd.DataFrame()
df_t['t1'] = [6,3,3]
df_t['t2'] = [1,1,2]
df_t['t3'] = [1,1,1]
df_t.index = [1,2,3]
df_t.plot.line(ax=ax)

# No norm
profile = []
for x in [1,2,3]:
    profile.append(df[df['meta'] == x].shape[0])

df_tt = pd.DataFrame()
df_tt['meta_no_norm'] = profile
df_tt.index = [1,2,3]
df_tt.plot.line(ax=ax)

# Norm
profile = []
for x in [1,2,3]:
    profile.append(df[df['meta'] == x].groupby('transcript').count().mean()[0])

df_tt = pd.DataFrame()
df_tt['meta_norm'] = profile
df_tt.index = [1,2,3]
df_tt.plot.line(ax=ax)


plt.close()

# Check that this method doesn't just collapse everything to one naturally and you will pick up end biases
fig = plt.figure()
ax = fig.add_subplot(111)

df = pd.DataFrame()
df['transcript'] = ['t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t1','t2','t2','t2','t2','t3','t3','t3']
df['meta'] = [3,3,3,3,3,3,3,1,1,1,2,2,2,3,3,3,3,3,3,1,2,3,3,3,3,3]

df_t = pd.DataFrame()
df_t['t1'] = [3,3,13]
df_t['t2'] = [1,1,2]
df_t['t3'] = [0,0,3]
df_t.index = [1,2,3]
df_t.plot.line(ax=ax)

# No norm
profile = []
for x in [1,2,3]:
    profile.append(df[df['meta'] == x].shape[0])

df_tt = pd.DataFrame()
df_tt['meta_no_norm'] = profile
df_tt.index = [1,2,3]
df_tt.plot.line(ax=ax)

# Norm
profile = []
for x in [1,2,3]:
    profile.append(df[df['meta'] == x].groupby('transcript').count().mean()[0])


df_tt = pd.DataFrame()
df_tt['meta_norm_ned'] = profile
df_tt.index = [1,2,3]
df_tt.plot.line(ax=ax)
