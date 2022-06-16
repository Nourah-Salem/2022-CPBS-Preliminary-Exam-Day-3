import numpy as np
import pandas as pd

df1 = pd.read_csv('D:/PhD at Anschutz/semester2/7712/prelim/data1.csv')
df2 = pd.read_csv('D:/PhD at Anschutz/semester2/7712/prelim/data2.csv')

def count_kmers(read, k):
    """Count kmer occurrences in a given read.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    Examples
    --------
    >>> count_kmers("GATGAT", 3)
    {'ATG': 1, 'GAT': 2, 'TGA': 1}
    """
    # Start with an empty dictionary
    counts = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    return counts


df_0 = df2[:1529]
df_1 = df2
df_2 = df1[:500]


allkmer_0 = []
for i in df_0['read']:
    allkmer_0.append(i)
allkmer_0 =str(allkmer_0)
allkmer_0 =allkmer_0.replace(",","")
allkmer_0 =allkmer_0.replace(" ","")
allkmer_0 =allkmer_0.replace("'","")
allkmer_0 =allkmer_0.replace("[","")
allkmer_0 =allkmer_0.replace("]","")
allkmer_0  = allkmer_0[2:] +  allkmer_0[:-2] 

kmer_0= count_kmers(allkmer_0,8)




allkmer_1 = []
for i in df_1['read']:
    allkmer_1.append(i)
allkmer_1 =str(allkmer_1)
allkmer_1 =allkmer_1.replace(",","")
allkmer_1 =allkmer_1.replace(" ","")
allkmer_1 =allkmer_1.replace("'","")
allkmer_1 =allkmer_1.replace("[","")
allkmer_1 =allkmer_1.replace("]","")
allkmer_1  = allkmer_1[2:] +  allkmer_1[:-2] 

kmer_1= count_kmers(allkmer_1,3)



allkmer_2 = []
for i in df_2['read']:
    allkmer_2.append(i)
allkmer_2 =str(allkmer_2)
allkmer_2 =allkmer_2.replace(",","")
allkmer_2 =allkmer_2.replace(" ","")
allkmer_2 =allkmer_2.replace("'","")
allkmer_2 =allkmer_2.replace("[","")
allkmer_2 =allkmer_2.replace("]","")
allkmer_2  = allkmer_2[2:] +  allkmer_2[:-2] 

kmer_2= count_kmers(allkmer_2,8)


'''
Normalised freqeuncy of kᵢ
= Number of occurrences of kᵢ / total number of k-mers
(where kᵢ is the iᵗʰ k-mer)'''

sum_kmers0 = 0
for index, key in enumerate (kmer_0):
    sum_kmers0 = sum_kmers0+kmer_0[key]
    
    
Norm_f_kmer0 = []
for index, key in enumerate (kmer_0):
    Norm_f_kmer0.append(kmer_0[key]/sum_kmers0)


sum_kmers1 = 0
for index, key in enumerate (kmer_1):
    sum_kmers1 = sum_kmers1+kmer_1[key]
    
    
Norm_f_kmer1 = []
for index, key in enumerate (kmer_1):
    Norm_f_kmer1.append(kmer_1[key]/sum_kmers1)
    

sum_kmers2 = 0
for index, key in enumerate (kmer_2):
    sum_kmers2 = sum_kmers2+kmer_2[key]
    
    
Norm_f_kmer2 = []
for index, key in enumerate (kmer_2):
    Norm_f_kmer2.append(kmer_2[key]/sum_kmers2)

    

Norm_f_kmer0.sort(reverse=True)
Norm_f_kmer1.sort(reverse=True)
Norm_f_kmer2.sort(reverse=True)

sample100_0 =  Norm_f_kmer0[:50]
sample100_1 =  Norm_f_kmer1[:50]
sample100_2 =  Norm_f_kmer2[:50]

    
#Generate a dummy dataset.
X = np.column_stack([sample100_0, sample100_1,sample100_2])


import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import pandas as pd
 
def PCA(X , num_components):
     
    #Step-1
    X_meaned = X - np.mean(X , axis = 0)
     
    #Step-2
    cov_mat = np.cov(X_meaned , rowvar = False)
     
    #Step-3
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)
     
    #Step-4
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]
     
    #Step-5
    eigenvector_subset = sorted_eigenvectors[:,0:num_components]
     
    #Step-6
    X_reduced = np.dot(eigenvector_subset.transpose() , X_meaned.transpose() ).transpose()
     
    return X_reduced


X_reduced = PCA(X , 3)


#Creating a Pandas DataFrame of reduced Dataset
principal_df = pd.DataFrame(X_reduced , columns = ['PC1','PC2','PC3'])
 
#Concat it with target variable to create a complete Dataset
# principal_df = pd.concat([principal_df , pd.DataFrame(target)] , axis = 1)


 
plt.figure(figsize = (6,6))
sb.scatterplot(data = principal_df , x = 'PC2',y = 'PC3' , s = 60 , palette= ['blue','red', 'green'])
