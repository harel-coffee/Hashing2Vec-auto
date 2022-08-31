#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
from sklearn.decomposition import TruncatedSVD
import random
# import seaborn as sns
import os.path as path
import os
# import matplotlib
# import matplotlib.font_manager
# import matplotlib.pyplot as plt # graphs plotting
# import Bio
from Bio import SeqIO # some BioPython that will come in handy
#matplotlib inline

from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score


# from matplotlib import rc
# # for Arial typefont
# matplotlib.rcParams['font.family'] = 'Arial'

from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
from sklearn import svm

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

import pandas as pd
import sklearn
from sklearn import preprocessing
from sklearn.model_selection import train_test_split 
from sklearn.preprocessing import StandardScaler  
from sklearn.neural_network import MLPClassifier 
from sklearn.metrics import classification_report, confusion_matrix 

from sklearn.neighbors import KNeighborsClassifier

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

from pandas import DataFrame

from sklearn.model_selection import KFold 
from sklearn.model_selection import RepeatedKFold

from sklearn.metrics import confusion_matrix

from numpy import mean

import seaborn as sns

import itertools
from itertools import product

import csv

from sklearn.model_selection import ShuffleSplit # or StratifiedShuffleSplit

from sklearn.decomposition import KernelPCA

import timeit
from fnvhash import fnv1a_32


print("done")


# In[2]:


seq_data = np.load("E:/RA/IJCAI/Dataset/Original/seq_data_7000.npy")
attribute_data = np.load("E:/RA/IJCAI/Dataset/Original/seq_data_variant_names_7000.npy")

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

def hash(item,seed,m_val):
    return fnv1a_32(item.encode(),seed) % m_val


# In[36]:





# unique_seq_kmers_final_list = [''.join(c) for c in product('ACDEFGHIKLMNPQRSTVWXY*', repeat=kmer_length)] 

# # n = 2  # no of items to add
# # N = 4  # size of the each counter
# m = len(unique_seq_kmers_final_list)  # total number of the buckets
# k = 1  # number of hash functions




# unique_value_cnt_holder = 0

# optimal_bin_size = m

# while unique_value_cnt_holder<len(unique_seq_kmers_final_list):
#     hash_int_value = []
#     for i in range(len(unique_seq_kmers_final_list)):
#         index = hash(unique_seq_kmers_final_list[i],k-1,optimal_bin_size)
#         hash_int_value.append(index)
#     #     print("Value: ",word_present[i],", Index value: ",index)
#     unique_value_cnt_holder = len(np.unique(hash_int_value))
#     print("Unique Hash Values: ",unique_value_cnt_holder, ", Optimal Bin Size:",optimal_bin_size,", Total Unique k-mers:",len(unique_seq_kmers_final_list))
#     optimal_bin_size = optimal_bin_size+100


# In[38]:


# optimal_bin_size = 404048
# for i in range(len(unique_seq_kmers_final_list)):
#     index = hash(unique_seq_kmers_final_list[i],k-1,optimal_bin_size)
#     hash_int_value.append(index)
    
# hash_int_value[0]


# In[ ]:


num_hash_fun = 1
kmer_length = 3 #
frequency_vector = []

# unique_seq_kmers_final_list = [''.join(c) for c in product('ACDEFGHIKLMNPQRSTVWXY*', repeat=kmer_length)] 

for seq_ind in range(len(seq_data)):
    print("index: ",seq_ind,"/",len(seq_data))
    se_temp = seq_data[seq_ind]
    kmers_list = build_kmers(se_temp,kmer_length)



    #create dictionary
    idx = pd.Index(kmers_list) # creates an index which allows counting the entries easily
    # print('Here are all of the viral species in the dataset: \n', len(idx),"entries in total")
    aq = idx.value_counts()
    counter_tmp = aq.values
    kmers_tmp = aq.index
    # counter_tmp,gmers_tmp


    #create frequency vector
    #cnt_check2 = 0
    optimal_bin_size = 404048 # for k=3
    listofzeros = [0] * optimal_bin_size
    ##################################################
    for ii in range(len(kmers_tmp)):
        index = hash(kmers_tmp[ii],num_hash_fun-1,optimal_bin_size)
        listofzeros[index] = counter_tmp[ii]
#         hash_int_value.append(index)
    ##################################################
    
#     for ii in range(len(gmers_tmp)):
#         seq_tmp = gmers_tmp[ii]
#     #     listofzeros = [0] * len(unique_seq_kmers_final_list)
#     #     for j in range(len(seq_tmp)):
#         ind_tmp = unique_seq_kmers_final_list.index(seq_tmp)
#         listofzeros[ind_tmp] = counter_tmp[ii]
    frequency_vector.append(listofzeros)


# In[7]:


np.save("E:/RA/IJCAI/Dataset/Original/Hasing_kmers_7000_Freq_vec_k_" + str(kmer_length)  +".npy",frequency_vector)


# In[ ]:





# In[ ]:




