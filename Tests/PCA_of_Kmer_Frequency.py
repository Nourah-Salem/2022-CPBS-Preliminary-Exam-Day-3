import numpy as np
import pandas as pd

'''theis script takes each DNA seq, cut it into Kmers with the size that the user choses
after cutting the Kmers, the frequency of each kmer across the genome is calculated and normalized.
Next, we calculated the PCA for 2 genomes through calculating the covariance to get the eigen_values and eigen_vectors and then calculate the PCA 
'''

df1 = pd.read_csv('./Processed_data/SARS_Cov2_data.csv')
df2 = pd.read_csv('./Processed_data/bacterial_data.csv')

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


# we took a 500 sequences from each community to make equal length columns
df_0 = df1[:500]
df_1 = df2[:500]

# create the label for each community to assign them in the PCA
label_vec = ["SARS Cov2"]*25
label_vec = label_vec + ["Escherichia coli"]*25


# additional processing of the sequences before cutting the kmers:
# Kmer for the viral community
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

kmer_0= count_kmers(allkmer_0,3)

#Kmer for the baterial comuity 
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


'''
Normalised freqeuncy of kᵢ
= Number of occurrences of kᵢ / total number of k-mers
(where kᵢ is the iᵗʰ k-mer)'''

# applied of the SARS-CoV2
sum_kmers0 = 0
for index, key in enumerate (kmer_0):
    sum_kmers0 = sum_kmers0+kmer_0[key]
    
    
Norm_f_kmer0 = []
for index, key in enumerate (kmer_0):
    Norm_f_kmer0.append(kmer_0[key]/sum_kmers0)
    

# applied on the bacterial community
sum_kmers1 = 0
for index, key in enumerate (kmer_1):
    sum_kmers1 = sum_kmers1+kmer_1[key]
    
    
Norm_f_kmer1 = []
for index, key in enumerate (kmer_1):
    Norm_f_kmer1.append(kmer_1[key]/sum_kmers1)
    

# selecting the kmers that has the highest frequency to represent each community
Norm_f_kmer0.sort(reverse=True)
Norm_f_kmer1.sort(reverse=True)
sample_0 =  Norm_f_kmer0[:50]
sample_1 =  Norm_f_kmer1[:50]

    
#combine the 2 comunities into one array.
X = np.column_stack([sample_0, sample_1])


import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import pandas as pd
 
''' we wrote steps of the PCA as follows'''
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


X_reduced = PCA(X , 2)


#Creating a Pandas DataFrame of reduced Dataset
principal_df = pd.DataFrame(X_reduced , columns = ['PC1','PC2'])
  
plt.figure(figsize = (10,5))
sb.scatterplot(data = principal_df , x = 'PC1',y = 'PC2' , s = 60 , hue = label_vec , palette= ['blue','red'])
plt.savefig('./Output/pca_sklearn.png', format="png", dpi = 300, bbox_inches='tight')


