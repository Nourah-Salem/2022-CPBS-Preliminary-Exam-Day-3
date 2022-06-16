from sklearn.decomposition import PCA
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

'''this script is confirm the results we get from writing the PCA model from scratch.
Here, use the funtion PCA imported from sklearn.decomposition'''

df1 = pd.read_csv('./Processed_data/SARS_Cov2_data.csv')
df2 = pd.read_csv('./Processed_data/bacterial_data.csv')

def count_kmers(read, k):
    
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

    
#Generate a dummy dataset.
pca = PCA(n_components=2)
embedded = pca.fit_transform(X)

fig = plt.figure(figsize=(10, 5))
sns.scatterplot(embedded[:,0], embedded[:,1],s = 60 , hue = label_vec , palette= ['blue','red']).plot()   
plt.savefig('./Output/pca_sklearn.png', format="png", dpi = 300, bbox_inches='tight')
print ("PCA done and image saved in the Output file")
