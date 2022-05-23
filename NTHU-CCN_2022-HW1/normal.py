# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 08:23:56 2022
@author: Paul
@file:   normal.py
@dependencies: 
    conda env:  pytest
    python:     3.7.11
    numpy:      1.21.2
    matplotlib: 3.5.0
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    # hyper-parameters
    SAMPLE_SIZE = 10000 # total number of smaples, which is the output shape
    BIN_SIZE = 100      # the number of equal-width bins in the range


    # numpy.random.normal(loc=0.0, scale=1.0, size=None)
    # params:  
    #   loc:     float(s)          -Mean, centre
    #   scale:   float(s)          -Standard deviation, spread
    #   size:    int, (ints..)     -Output shape
    # returns: 
    #   samples: ndarray or scalar -Drawn samples

    # mean = 0 and standard deviation = 1, which means variance also equal to 1
    mu, sigma = 0.0, 1.0
    # store all the samples of a normal distribution to samples
    samples = np.random.normal(loc=mu, scale=sigma, size=SAMPLE_SIZE) 
    print(type(samples)) # <class 'numpy.ndarray'>

    # verify the mean and the variance (standard deviation) 
    print('mean difference:', abs(mu - np.mean(samples)))
    # delta degrees of freedom, but I can't really tell the difference lol
    # set ddof=0 return population standard deviation, biased estimation
    # set ddof=1 return sample standard deviation, unbiased estimation
    print('standard deviation difference:', abs(sigma - np.std(samples, ddof=1)))


    # matplotlib.pyplot.hist(x, bins=None, density=False, stacked=False, ...)
    # params:     
    #   x:       (n,) array(s)     -Input values, can take single array or sequence of arrays
    #   bins:    int, seq. or str  -If int, it defines the number of equal-width bins in range
    #   density: bool              -If True, area under the histogram integrates to 1
    #   stacked: bool              -If True, multiple data are stacked on top of each other
    # returns:    
    #   counts:  array(s)          -The values of the histogram bins
    #   bins:    array             -The edges of the bins
    #   patches: BarContainer      -Just ignored, when there's only one input dataset

    plt.figure(1)
    # plot the pmf of the samples
    count, bins, ignored = plt.hist(x=samples, bins=BIN_SIZE, density=True, stacked=True)
    # plot the actual normal distribution
    plt.plot(bins, 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-(bins - mu)**2 / (2*(sigma**2))), 
             linewidth=2, color='r')

    plt.xlabel('random values')
    plt.ylabel('probability')
    # format string f'...{SAMPLE_SIZE}...' will substitute variable in runtime
    plt.title(f'The PMF over {SAMPLE_SIZE} random samples')
    plt.text(-5, 0.37, r'$(\mu = 0,\ \sigma = 1)$') # text(x-coord, y-coord, '{text display}')
    # be very careful on setting axis, cause it will dramatically affect how we perceive data
    # plt.axis([-5.5, 5.5, 0, 0.42])
    plt.grid(True)

    # first check whether the specified path is an existing file
    path = './normal_distribution/normal.png'
    # if no file exist, then save current file
    if os.path.isfile(path) is False:
        plt.savefig('./normal_distribution/normal.png')
    plt.show()

# main()
if __name__ == '__main__':
    main()


    
    
