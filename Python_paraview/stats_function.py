# plot all columns of aspect statistics file and save figures
# v2 loops over multiple folders

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


### INPUT ###

foldernamelist = glob.glob('/scratch/projects/bbp00064/bbkkaili/0001_com_extmodels/111_007_0001_0.25stronglowercrust/001_25_15/003_nosurface/002_30com/0005_density2896/001_totalerosion/004_timestepasnormal/output/')


print ("foldernamelist")
for foldername in foldernamelist:
    
    # assign aspect folder containing statistics file
    #foldername = foldernamelist[0]
    filename = foldername + 'statistics'
    outfoldername = foldername +'statistics_plots'
    
    
    ### CREATE OUTPUT FOLDER (IF IT DOES NOT EXIST YET) ###
    import os, errno
    try:
        os.makedirs(outfoldername)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    
    
    ### EXTRACT VARIABLE NAMES ###
    # read entire file and store each line
    with open(filename) as f:
        header = f.readlines()  
    # remove all lines that do not start with "#" (starting from the back)
    nonheaderID = [x[0]!='#' for x in header] 
    for index,linecontent in reversed(list(enumerate(header))):
        if linecontent[0]!='#':
            del header[index]
    # remove whitespace characters like `\n` at the end of each line
    header = [x.strip() for x in header] 
    
    
    ### EXTRACT DATA
    df = pd.read_csv(filename, comment='#', header=None, delim_whitespace=True)
    step = np.array(df.iloc[:,0])
    time = np.array(df.iloc[:,1])
    
    
    ### PLOT DATA ###
    # prepare second axis
    NumberOfTicks = 6
    dtick = int(round( float(len(time)) / NumberOfTicks ))
    xticks2 = []
    for i,t in enumerate(time):
        if i%dtick == 0:
            xticks2.append(t)
    
    for i in range(len(header)):
    #for i in range(1):
        # check if column data is plottable (not a string)
        if type(df.iloc[0,i]) is str:
            continue
        # plot wrt steps
        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(111)
        ax1.plot(step,df.iloc[:,i], color='k')
        ax1.set_xlabel(header[0])
        ax1.set_ylabel(header[i])
        ax1.set_xlim((0,step[-1]))
        # plot wrt time as second x-axis
        ax2 = ax1.twiny()
        ax2.plot(time, df.iloc[:,i], color='royalblue')
        ax2.set_xlabel(header[1], color='royalblue')
        # x1ticks = ax1.get_xticks
        # x1ticklabels = ax1.get_xticklabels
        ax2.set_xticks(xticks2)
        ax2.xaxis.label.set_color('royalblue')
        ax2.spines['top'].set_color('royalblue')
        ax2.tick_params(axis='x', colors='royalblue')
        ax2.set_xlim((0,time[-1]))
        # save figures with leading zeros
        plt.savefig(outfoldername + '/' + str(i+1).zfill(3) + '.png', dpi=200)
        plt.show()
   
        # extract maximum velocity and dt
        if header[i] == '# 3: Time step size (years)':
            dt=df.iloc[:,i]
        if header[i] == '# 19: Max. velocity (m/year)':
            maxvel=df.iloc[:,i]

## PLot CFL number through time
#dx=5000 # m
#CFL = maxvel * dt /dx
## plot wrt steps
#fig = plt.figure(figsize=(8,4))
#ax1 = fig.add_subplot(111)
#ax1.plot(step,CFL, color='k')
#ax1.set_xlabel(header[0])
#ax1.set_ylabel('CFL number (-)')
#ax1.set_xlim((0,step[-1]))
## plot wrt time as second x-axis
#ax2 = ax1.twiny()
#ax2.plot(time, CFL, color='royalblue')
#ax2.set_xlabel(header[1], color='royalblue')
## x1ticks = ax1.get_xticks
## x1ticklabels = ax1.get_xticklabels
#ax2.set_xticks(xticks2)
#ax2.xaxis.label.set_color('royalblue')
#ax2.spines['top'].set_color('royalblue')
#ax2.tick_params(axis='x', colors='royalblue')
#ax2.set_xlim((0,time[-1]))
## save figures with leading zeros
#plt.savefig(outfoldername + '/' + str(i+2).zfill(3) + '.png', dpi=200)
#plt.show()
#
## PLot dx/dt number through time
#dx=312 # m
#CFL = maxvel * dt /dx
## plot wrt steps
#fig = plt.figure(figsize=(8,4))
#ax1 = fig.add_subplot(111)
#ax1.plot(step,CFL, color='k')
#ax1.set_xlabel(header[0])
#ax1.set_ylabel('CFL number (-), Elem size ' + str(dx) + ' m')
#ax1.set_xlim((0,step[-1]))
## plot wrt time as second x-axis
#ax2 = ax1.twiny()
#ax2.plot(time, CFL, color='royalblue')
#ax2.set_xlabel(header[1], color='royalblue')
## x1ticks = ax1.get_xticks
## x1ticklabels = ax1.get_xticklabels
#ax2.set_xticks(xticks2)
#ax2.xaxis.label.set_color('royalblue')
#ax2.spines['top'].set_color('royalblue')
#ax2.tick_params(axis='x', colors='royalblue')
#ax2.set_xlim((0,time[-1]))
## save figures with leading zeros
#plt.savefig(outfoldername + '/' + str(i+2).zfill(3) + '.png', dpi=200)
#plt.show()
