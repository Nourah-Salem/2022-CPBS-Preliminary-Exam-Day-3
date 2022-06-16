

from FASTA_CSV import FASTA2DF
from PCA_of_Kmer_Frequency import count_kmers
from PCA_of_Kmer_Frequency import PCA 

import unittest
import pandas as pd
import pandas.api.types as ptypes
import re
import numpy as np

class project_tests(unittest.TestCase):
    
    def test_FASTA2DF_return(self):
        
        '''this test function tests the correctness of the dataframe generated from FASTA (and the same goes for the FASTQ function), it tests the type of the columns that were generated from FASTA file'''
        
        # expected
         
        ''' the shape of the output dataframe is similar to the example below:
            
        df = pd.DataFrame({
            
            'id': ['text'],
            
            'read': ['text'],
            
            'length': [int]}) '''
        
        #actual
        df = FASTA2DF("D:/PhD at Anschutz/semester2/7712/prelim/bacterial_data.fasta")
        
        # assert
        assert ptypes.is_string_dtype(df['id'])
        assert ptypes.is_string_dtype(df['read'])
        assert ptypes.is_numeric_dtype(df['length'])
        
    
    
    def read_char_check(self):
        
        ''' this test function tests if the reads have any character other than A, T, C, G'''
        
        non_DNA = re.compile('[QWERYUIOP{}SDFHJKL:"?><MNBVXZ@_!#$%^&*()<>?/\|}{~:]')
        df = FASTA2DF("D:/PhD at Anschutz/semester2/7712/prelim/bacterial_data.fasta")
        
        # chech if the reads have any og the characters in the non_DNA list
        if(non_DNA.search(df['read'][1]) == None):
            # generate True variable of a read has any special char
            X = True
         
        else:
            # generate False variable of a read doesn't have any special char
            X = False
            
        self.assertTrue(X)
        
        return None

    
    def num_kmer_generated(self):
        
        '''this test function checks the number of the kmers generated from out Kmer builder function the correct number of the kmer depend on have only the DNA sequence, if there is any additional special character (e.g "{"), the test will fail'''
        
        read = 'AATTCCGGAA'
        k= 8
        kmer= count_kmers(read , k)
        
        
        # expected
        ''' the number of kmer follows the equation: (count_kmer = R - k + 1), R is 
        the length of the unput sequence, k is the window size of the kemr.'''
        
        count_kmer = len(read) - k +1
        
        # assert
        self.assertEqual(kmer,count_kmer)
        
        return None
        
    
    def PCA_input_output_same_shape_check(self):
        
        '''this test function checks if the shape of the input to the PCA is correct or not.
        Here, we compare the shape of the PCA input to an idea shpae that we generated from random data'''
        
        # expected 
        input_shape = np.random.randint(10,50,60).reshape(50,2) 
        
        output_shape = PCA(input_shape)
        
        # assert
        self.assertAlmostEqual(input_shape.shape, output_shape.shape)
        
        return None
    

if __name__ == '__main__':
    unittest.main()
