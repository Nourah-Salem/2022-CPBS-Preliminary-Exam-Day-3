''' Coverting the list of reads for the bacterial genome from the FASTA format to the FASTQ format'''

import pandas as pd

def FASTA2DF (file):
    # The function takes the FASTA file, convert it to 3-column dataframe (one for the ID, one for the read and the length of the read)
    
    # open file and iterate through the lines, composing each single line as we go
    id_lines = []
    seq_lines = []
    with open(file) as fp:
         for line in fp:
             if line.startswith('>'):
                 id_lines.append(line)
             else:
                 line.replace('\n', "")
                 seq_lines.append(line)
    
    # Calling DataFrame constructor after zipping
    # both lists, with columns specified
    df = pd.DataFrame(list(zip(id_lines, seq_lines)),
                   columns =['id', 'read']).replace(r'\n','', regex=True)
    
    #Calculates the length of each read, append the length
    df1= df.read.str.len().tolist()
    df1 = pd.Series(df1 , name = 'length')
    
    
    #Sort the reads descendingly
    df2 = pd.concat([df, df1], axis=1, join="inner").sort_values("length",ascending=False)
    df2.reset_index(drop=True, inplace=True)
    
        
    return df2


        
if __name__ == "__main__":
    
    df = FASTA2DF("./Raw_Data/data.fasta")
    df = df.sort_values(by='id')
    df.reset_index(drop=True, inplace=True)    
    df.to_csv("./Processed_data/data_bacteria.csv")
