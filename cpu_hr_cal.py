'''
 This script is meant to calculate the number of
cpu hours run on super computers from the output
of the sacct command. e.g:

sacct --user=x_sepki -X --format=JobID,USER,Submit,Start,End,Elapsed,Nodelist,State,AllocCPUS%9,CPUTimeRaw --allocations --state=ca,cd,f,nf,to,r --starttime=2019-10-06 --endtime=2019-10-24 --state=completed

adding (> filename.txt) at the end of the command will store the information in the specified filename.txt
'''

import numpy as np

filename = '/home/seki2695/test.txt'

# open file in read-only format
file = open(filename, 'r')
time = []
node = []
i = 0
for line in file:
    test = line.split()
    if (test[7]=='COMPLETED'):
        dd_hh = (test[5].split(':',3)[0]).split('-',1)
        if (len(dd_hh))==2:
            dd = float(dd_hh[0]) # days
            hh = float(dd_hh[1]) # hours
        else:
            dd = 0
            hh =  float(dd_hh[0])
        mm = float(test[5].split(':',3)[1]) # minutes
        ss = float(test[5].split(':',3)[2]) # seconds
        node.append(int(test[8]))
        #print(dd, hh, mm, ss)
        time.append(dd*24 + hh + mm/60. + ss/3600.) # time array in hours

# converting to numpy array
time = np.array(time, dtype = float) 
node = np.array(node, dtype = float)

print ('Hours of processing on each node     = ' + str(np.round(np.sum(time))) +
       '\nAverage inversion time / sb-region   = ' + str(np.round(np.mean(time), decimals = 2)) +
       '\nAverage number of cpu-cores per node = ' +str(int(np.round(np.mean(node)))) +
       '\nTotal number of regions              = ' + str(len(time)) + 
       '\nTotal number of CPU hours            = ' + str(np.round(np.sum(time*node)))
)
