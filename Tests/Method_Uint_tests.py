

from EDA_sample_prep import FASTA2DF
from Kmer_builder import count_kmers
from Kmer_builder import PCA 

import unittest
import pandas as pd
import pandas.api.types as ptypes
import re
import numpy as np

class project_tests(unittest.TestCase):
    
    def test_FASTA2DF_return(self):
        
        # expected
         
        ''' the shape of the output dataframe is similar to the example below:
            
        df = pd.DataFrame({
            
            'id': ['text'],
            
            'read': ['text'],
            
            'length': [int]}) '''
        
        #actual
        df = FASTA2DF("D:/PhD at Anschutz/semester2/7712/prelim/data2.fasta")
        
        # assert
        assert ptypes.is_string_dtype(df['id'])
        assert ptypes.is_string_dtype(df['read'])
        assert ptypes.is_numeric_dtype(df['length'])
        
    
    
    def read_char_check(self):
        
        non_DNA = re.compile('[QWERYUIOP{}SDFHJKL:"?><MNBVXZ@_!#$%^&*()<>?/\|}{~:]')
        df = FASTA2DF("D:/PhD at Anschutz/semester2/7712/prelim/data2.fasta")
        
        if(non_DNA.search(df['read'][1]) == None):
            X = True
         
        else:
            X = False
            
        self.assertTrue(X)
        
        return None

    
    def num_kmer_generated(self):
        
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
        
        # expected 
        input_shape = np.random.randint(10,50,60).reshape(20,3) 
        
        output_shape = PCA(input_shape)
        
        # assert
        self.assertAlmostEqual(input_shape.shape, output_shape.shape)
        
        return None
    

if __name__ == '__main__':
    unittest.main()
