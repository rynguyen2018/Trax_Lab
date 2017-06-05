from __future__ import division
from pyteomics import mzxml
import sys
import math
from operator import itemgetter
import numpy as np
import csv
from itertools import izip
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt


input_file=sys.argv[1]





def PCA(data, dims_rescaled_data=2):
    """
    returns: data transformed in 2 dims/columns + regenerated original data
    pass in: data as 2D NumPy array
    """
    from scipy import linalg as LA
    m, n = data.shape
    # mean center the data
    # data -= data.mean(axis=0)
    # calculate the covariance matrix
    R = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since R is symmetric, 
    # the performance gain is substantial
    evals, evecs = LA.eigh(R)
    # sort eigenvalue in decreasing order
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    # sort eigenvectors according to same index
    evals = evals[idx]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
    evecs = evecs[:, :dims_rescaled_data]
    # carry out the transformation on the data using eigenvectors
    # and return the re-scaled data, eigenvalues, and eigenvectors
    return np.dot(evecs.T, data.T).T, evals, evecs


def plot_pca(data):
    from matplotlib import pyplot as MPL
    print "Creating plot for data..."
    clr1 =  '#2026B2'
    fig = MPL.figure()
    ax1 = fig.add_subplot(111)
    data_resc, data_orig = PCA(data)
    ax1.plot(data_resc[:, 0], data_resc[:, 1], '.', mfc=clr1, mec=clr1)
    MPL.show()

'''Retrieves all MS2 spectra from a sample file. log-normalize all MS2 base peak intensities.'''
def preprocess_sample(sample):
	print "Reading " + sample
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break

	print str(len(scans)) + " scans found in " + sample
	base_peaks = {}
	all_peaks = []
	for scan in scans:
		if scan['msLevel'] == '2':
			ms1_scan_num= int(scan['precursorMz'][0]['precursorScanNum'])
			base_mz = scan['precursorMz'][0]['precursorMz']
			precursor_intensity = scan['precursorMz'][0]['precursorIntensity']
			intensity_array=scan['intensity array'].tolist()
			msI_TIC= scans[ms1_scan_num]['totIonCurrent']
			ms2_TIC= scan['totIonCurrent']
			percentComposition= ms2_TIC/msI_TIC
			#Normalize and log transform peak intensities in each scan
			#intensities = normalize(np.log(1+np.asarray(scan['intensity array']).reshape(1,-1)), norm='l1')[0]
			for x in range(0,len(intensity_array)): 
			 	intensity_array[x]= intensity_array[x]/ms2_TIC	
			mzs = scan['m/z array']
			num = int(scan['num'])
			percentComposition= ms2_TIC/msI_TIC
			base_peaks[num] = {"num":num, "base_mz":base_mz, "intensities":intensity_array, "mzs":mzs, "precursor_intensity": precursor_intensity,"Percent of Sample": percentComposition, "MS2 TIC": ms2_TIC}

			all_peaks = all_peaks + mzs.tolist()

	peak_min = int(math.floor(min(all_peaks)))
	peak_max = int(math.ceil(max(all_peaks)))

	#Get rid of really big variables...
	all_peaks = None
	scans = None
	r = None
	
	#Returns a list of MS2 spectra organized by scan #, and the largest and smallest precursor peaks across all MS2 in this file.
	return peak_min, peak_max, base_peaks

""" 
	Orders Peak data by percent composition of identified compounds that have MSII 
	After finding clustering of molecules, maybe it could look at its clustering and see which one has the highest abundance of such a molecule,
	or rather, if it is abundant at all? 

"""
def orderedPercentComp(new_peak_data): 
	percentComp_list=[]
	for key in new_peak_data.keys(): 
		percentComp_list.append((new_peak_data[key]['num'],new_peak_data[key]['Percent of Sample'],new_peak_data[key]['MS2 TIC']))
	lis= sorted(percentComp_list,key=lambda x: x[1], reverse=True)
	return lis


def vectorize_peak(peak_min, peak_max, sample_data, sample_name):	

	#initializes peak vector (intervals of 0.1 so multiple all values by 10)
	vector_length = ( peak_max - peak_min ) * 10
	peak_vectors_list = []
	print "Creating peak vectors for " + sample_name

	peak_vectors = {}
	new_peaks_list = sorted(sample_data.values(), key=itemgetter('base_mz'))
	#Creates peak vectors for each scan in this sample
	for scan in new_peaks_list:
		peak_vector = np.zeros(vector_length, dtype=np.float32)
	 	i = 0
		for p in scan['mzs']:
			pos = int((math.floor(p*10)/10 - peak_min) * 10)
			peak_vector[pos] = scan['intensities'][i]
			i += 1

		peak_vectors[scan['num']] = peak_vector

	  	peak_vectors_list.append(scan['num'])
	return peak_vectors


peak_min = 999999999
peak_max = 0
	#peak_data = {}
new_min, new_max, new_peak_data = preprocess_sample(input_file)
if new_min < peak_min:
	peak_min = new_min
if new_max > peak_max:
	peak_max = new_max
ordered_Peaks= orderedPercentComp(new_peak_data)

peak_data = vectorize_peak(peak_min, peak_max, new_peak_data, input_file)
value_list=[]
for key in peak_data.keys():
	 value_list.append(peak_data[key])
peak_data_2d= np.array(value_list)
plot_pca(peak_data_2d)
