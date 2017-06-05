from pyteomics import mzxml
import sys
from itertools import groupby
import collections
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns

scans= []
sample=sys.argv[1]

r = mzxml.read(sample)
while True:
	try:
		scans.append(r.next())
	except:
		break
freq_list=[]
print "making frequency list"
for scan in scans: 
	if scan['msLevel'] == '2':
		base_mz = int(scan['precursorMz'][0]['precursorMz']*10000)
		precursor_intensity = scan['precursorMz'][0]['precursorIntensity']
		freq_list.append(base_mz)
freq_list=sorted(freq_list)
d = {x:freq_list.count(x) for x in freq_list}
od =sorted(d.items(), key=lambda x: x[1])
print freq_list

#sns.distplot(d.keys())
plt.hist([item/float(10000) for item in freq_list],bins=len(freq_list))
#plt.bar(list(d.keys()), d.values(), color='r')
plt.show()
#plt.savefig('test')

#print freq_list
#freq_list=[len(list(group)) for key, group in groupby(freq_list)]
#print sorted(freq_list)