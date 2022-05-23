from __future__ import division
# import
#from matplotlib import rc
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
import sys
import random
from numpy import *
from scipy.special import *
#import math
#import scipy.stats

# change font size in legend
#params = {'legend.fontsize': 20}
#plt.rcParams.update(params)

##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font', family='serif')
#rc('text', usetex=True)

### glbal variables
class glb:
    snrdBmin = 0                                    
    snrdBmax = 32
    snrdBinterval = 2
    R = 1 # tx rate
    T = pow(2,2*R)-1 # threshold
    tol_e = 0.00001 # tolerance of error
    MaxSamples = pow(10,5)

    FILE = True
    

### convert decibels to absolute values ###
def db2val(db):
    return pow(10,db/10.0)

class Cchannel:
    def get_chgain(lambd):        
        h = random.exponential(lambd)
        return h
    get_chgain = staticmethod(get_chgain)


def main():

    if glb.FILE:
        filename = "snrdB.txt"
        FILE0 = open(filename,"w")
        
        filename = "pout_vaf_sim.txt"
        FILE1 = open(filename,"w")

        filename = "pout_vaf_ana.txt"
        FILE2 = open(filename,"w")

        filename = "gain_vaf.txt"
        FILE3 = open(filename,"w")

    
    plist_pout_sim = []
    plist_pout_ana = []
    gainlist = []
    snrdBlist = range(glb.snrdBmin,glb.snrdBmax,glb.snrdBinterval)

    for snrdB in snrdBlist:
        
        snr = db2val(snrdB)
        
        # Theoretical
        pout_ana = 1-2*glb.T/np.sqrt(snr*snr)*exp(-glb.T*(1/snr+1/snr))*kn( 1, 2*glb.T/np.sqrt(snr*snr) )
        print(pout_ana)
        Stop = False
        n_err = 0
        #n_e = pow(1.96/glb.tol_e,2) # stopping criteria
        n_s1 = pout_ana*(1-pout_ana)*pow(1.96/ pow(10,-4), 2) # stopping criteria
        print(n_s1)
        n_s = 0 # num of samples

        G2_sum = 0.0
        while (not Stop):
            n_1 = 1/math.sqrt(2)*(np.random.normal(0,1)+np.random.normal(0,1)*1j) # whilte gaussian of direct hop
            h_1 = 1/math.sqrt(2)*(np.random.normal(0,1)+np.random.normal(0,1)*1j) # rayleigh channel of direct hop
            n_2 = 1/math.sqrt(2)*(np.random.normal(0,1)+np.random.normal(0,1)*1j) # whilte gaussian of direct hop
            h_2 = 1/math.sqrt(2)*(np.random.normal(0,1)+np.random.normal(0,1)*1j) # rayleigh channel of direct hop

            # compute SNR
            a_1 = pow(abs(h_1),2)
            a_2 = pow(abs(h_2),2)
            n_1 = pow( abs(pow(10,-snrdB/20.0)*n_1), 2)
            n_2 = pow( abs(pow(10,-snrdB/20.0)*n_2), 2)
            G2 = 1/(n_1+a_1)
            #snr_eq = a_1/n_1*a_2/n_2/(a_1/n_1+a_2/n_2+1)
            snr_eq = a_1/n_1*a_2/n_2/(a_2/n_2+1/G2/n_2)
            G2_sum += G2
            

            #snr_0 = Cchannel.get_chgain(1)*snr

            # count the errors
            #if np.log2(1+snr_0) < glb.R:
            if snr_eq < glb.T:
                n_err += 1
            

            n_s += 1
            #Stop = (n_s > glb.MaxSamples) or (n_err > n_e)
            Stop = (n_s > glb.MaxSamples) #(n_s > n_s1) or (n_s > glb.MaxSamples)
            #print("n_err=", n_err, "Stop=", Stop)
               
        #print(snr_sum/n_s)
        #print("snrdB=", snrdB, "Required errors=", n_e, "Actual runs=", n_s)        
        
        pout_sim = n_err/float(n_s)
        G2_avg = G2_sum/float(n_s)

        print("snrdB=", snrdB, "Required runs=", n_s1, "Actual runs=", n_s, "pout_sim=", pout_sim, "G2=", G2_avg)
        
        plist_pout_sim.append(pout_sim)
        plist_pout_ana.append(pout_ana)
        gainlist.append(G2_avg)

        

        if glb.FILE:
            FILE0.write('%s \n' % snrdB)
            FILE1.write('%s \n' % pout_sim)            
            FILE2.write('%s \n' % pout_ana)
            FILE3.write('%s \n' % G2_avg)

    if glb.FILE:
        FILE0.close()
        FILE1.close()
        FILE2.close()
        FILE3.close()

 
    # diversity order
    do = 1
    a1 = [glb.snrdBmin,1];
    a2 = [glb.snrdBmax,pow(10,(-do*glb.snrdBmax/10))];

    plt.figure(1)    
    line1 = plt.semilogy(snrdBlist,plist_pout_sim,label='Sim')
    plt.setp(line1, 'color', 'r', 'marker', 'o', 'markersize', 10, 'lw', 0)
    #plt.hold(True)
    line2 = plt.semilogy(snrdBlist,plist_pout_ana, label='Ana')
    plt.setp(line2, 'color', 'b', 'linestyle', '-', 'linewidth', 2)
    plt.legend(loc='upper right', shadow=True, fontsize=20)
    plt.xlabel('SNR [dB]',fontsize=18)
    plt.ylabel('Outage probability', fontsize=18)
    plt.grid(which='both')
    # set some Legend properties
    #lh = plt.gca().get_legend()
    #ltext = lh.get_texts()
    #plt.setp(ltext, fontsize='small')
    plt.ylim([10**-3,1])
    plt.xlim([0,30])
##    plt.hold(False)    
    plt.savefig("pout_vaf.eps", edgecolor='none', format='eps', bbox_inches='tight')
    plt.show()
    plt.close()
    #sys.exit(0)

    plt.figure(2)
    line3 = plt.semilogy(snrdBlist,gainlist)
    plt.xlabel('SNR (dB)', fontsize=18)
    plt.ylabel('Gain power', fontsize=18)
    plt.grid(which='both')
    plt.savefig("gain_vaf.eps", edgecolor='none', format='eps', bbox_inches='tight')
    plt.show()
    plt.close()

if __name__ == '__main__': 
    main()

##x = np.arange(0, 5, 0.1);
##y = np.sin(x)
##plt.plot(x, y)
##plt.show()
