import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def stats(foldernamelist):
    for foldername in foldernamelist:
        filename = os.path.join(foldername, 'statistics')
        outfoldername = os.path.join(foldername, 'statistics_plots')
        
        # Create output folder if it does not exist
        if not os.path.exists(outfoldername):
            os.makedirs(outfoldername)
        
        # Extract variable names from the header
        with open(filename) as f:
            header = f.readlines()  
        nonheaderID = [x[0]!='#' for x in header] 
        for index, linecontent in reversed(list(enumerate(header))):
            if linecontent[0]!='#':
                del header[index]
        header = [x.strip() for x in header] 
        
        # Extract data
        df = pd.read_csv(filename, comment='#', header=None, delim_whitespace=True)
        step = np.array(df.iloc[:,0])
        time = np.array(df.iloc[:,1])
        
        
        for i in range(len(header)):
            # Skip non-numeric columns
            if isinstance(df.iloc[0, i], str):
                continue
            