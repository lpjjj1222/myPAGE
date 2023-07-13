#!/usr/bin/env python
# coding: utf-8

# In[271]:


import utils
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import patsy
import seaborn as sns


# In[272]:


gene_pair_path = '/Users/lpjjj/Downloads/scPAGE-main-supplement/bulkPAGE/filtered_pairs_p_7.txt'
gene_pairs_res = pd.read_table(gene_pair_path,sep='\t')
num_pair = len(gene_pairs_res)

gene_pairs = []
for row in range(num_pair):
    gene_pairs.append(list(gene_pairs_res.iloc[row,]))
print(gene_pairs)


# In[273]:


#print(gene_pairs)
test_data_path='/Users/lpjjj/Downloads/scPAGE-main-supplement/model2/for_training/185263.txt'
test_data = pd.read_table(test_data_path,sep='\t', index_col=[0])
#print(test_data)

test_data = utils.profile_preprocessing(test_data)  # samples by genes matrix
samples = test_data.index
print("Number of genes in the test data: %s" % test_data.shape[1])
print("Number of samples in the test data: %s" % test_data.shape[0])


# In[274]:


def get_gene_list(gene_pairs):
    gene_list = []
    for pair in gene_pairs:
        gene_list.append(pair[0])
        gene_list.append(pair[1])
    gene_list = list(set(gene_list))
    return gene_list

def get_test_drop_gene(test_data,gene_pairs):
    gene_list = get_gene_list(gene_pairs)
    test_genes = test_data.columns
    drop_gene = [gene for gene in gene_list if gene not in test_genes]
    if len(drop_gene)!=0:
        print("%d genes in scGPS not detected in the test dataset" % len(drop_gene))
        print("They are: %s" % str(drop_gene))
    else:
        print("All genes in scGPS detected in the test dataset")
    return drop_gene

def pair_2_pair_index(drop_gene,pair_list):
    dele=[]
    pair_index = np.array([[pair[0] for pair in pair_list],[pair[1] for pair in pair_list]])  # reshape to 2 by x
    if len(drop_gene) > 0:
        for gene in drop_gene:
            for i in range(pair_index.shape[1]):
                if gene in pair_index[:,i]:
                    dele.append(i)
        dele = list(set(dele))
        new_index = np.delete(pair_index,dele,axis=1)
    else:
        new_index = pair_index

    return new_index


# In[275]:


drop_gene = get_test_drop_gene(test_data,gene_pairs)
pair_index=pair_2_pair_index(drop_gene,gene_pairs)


# In[276]:


def pairconvert_test(data,gene_pair_index):
    # 取基因对的差值，在取差值前将每一个数+1后取log2
    sub = np.array(np.log2(test_data.iloc[:][pair_index[1,:]]+1)) - np.array(np.log2(test_data.iloc[:][pair_index[0,:]]+1))  # gj - gi > 0 (gi, gj in gene_pair_index)
    return sub

rankdata = pairconvert_test(test_data,pair_index)
print(rankdata)

num_pair = rankdata.shape[1]
num_pair

print(pair_index)



# In[277]:


test_label_path='/Users/lpjjj/Downloads/scPAGE-main-supplement/model2/label/label_185263.txt'
test_label_df = pd.read_table(test_label_path)
test_label = np.array(list(test_label_df['Label']))
test_label


# In[278]:



X=pd.DataFrame(data=rankdata)
#print(X)
#X.columns=gene_pairs
# print(pair_index[0])
# print(pair_index[1])
pair_together_index=[]
pair=np.array(pair_index)[0]

pair_con_name=[]
for i in range(len(pair)):
    con=pair_index[0][i]+'_'+pair_index[1][i]
    pair_con_name.append(con)
tolist_con=list(pair_con_name)

print(tolist_con)
# print(type(pair_con_name))
# print(X)
X.columns=tolist_con
X    


# In[109]:


X['status']=test_label
print(X)
import openpyxl
X.to_excel('/Users/lpjjj/Downloads/scPAGE-main-supplement/bulkPAGE/185263_p7.xlsx')

