import os

def str(myfloat):
    return '%.5g'%myfloat

for fname in os.listdir('feb10/'):
    if 'hist_evolve_7.5' in fname:
        splits = fname.split('_')
        newname = '%.5f'%float(fname.split('_')[3])
        oldfloat = float(fname.split('_')[3])
        print(str(oldfloat))
        splits[3] = newname
        newfname = '_'.join(splits)
#        os.rename('feb10/'+fname,'feb10/'+newfname)
