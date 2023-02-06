####
## import module
# listmatfile
from os import listdir
from scipy import stats

# loatmatvar
from hdf5storage import loadmat
# avref, allspectra_z
import numpy as np
# fcoefLR
from scipy import signal 
from scipy.fftpack import fft
# calparcoh
from sklearn.covariance import GraphicalLassoCV

# skggmparcoh
import sys
sys.path.append("/home/zhibinz2/Documents/GitHub/skggm/inverse_covariance/")
from inverse_covariance import (
    QuicGraphicalLassoCV)

# save_image
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

####

## class
class Zclass:
    def listmatfile(dir):
        matfile = list()
        for file in listdir(dir):
	        if file.endswith(".mat"):
		        matfile.append(file)
        return matfile

    def loadmatvar(dir,matfile,i):
        data=loadmat(dir+matfile[i]) # load one session
        bpchan = int(data['bpchan'][0][0])
        channels = data['channels'][0]
        conditionNames = data['conditionNames'][0]
        conditions = data['conditions'][0]
        dataL = data['dataL'][0]
        dataR = data['dataR'][0]
        intervals = data['intervals'][0]
        labels = data['labels'][0]
        samples = data['samples'][0]
        session = int(data['session'][0])
        sessionTypes = data['sessionTypes'][0]
        sr = int(data['sr'][0])
        return data, bpchan,channels,conditionNames,\
            conditions,dataL,dataR,intervals,labels,\
                samples,session,sessionTypes,sr

    def loadmatvar2(dir,filename):
        data=loadmat(dir+filename) # load one session
        bpchan = int(data['bpchan'][0][0])
        channels = data['channels'][0]
        conditionNames = data['conditionNames'][0]
        conditions = data['conditions'][0]
        dataL = data['dataL'][0]
        dataR = data['dataR'][0]
        intervals = data['intervals'][0]
        labels = data['labels'][0]
        samples = data['samples'][0]
        session = int(data['session'][0])
        sessionTypes = data['sessionTypes'][0]
        sr = int(data['sr'][0])
        return data, bpchan,channels,conditionNames,\
            conditions,dataL,dataR,intervals,labels,\
                samples,session,sessionTypes,sr

    def avref(datax):
        ref = np.mean(datax,axis = 1)
        refmat = np.tile(ref,(np.shape(datax)[1],1))
        datax = datax - np.transpose(refmat)
        return datax

    def fcoefLR(dataL,dataR,trl):
        # trl: select trial number
        # chan number
        nchan=32
        # sr
        sr=2000
        # zscore dataL and dataR
        dataL[trl] = stats.zscore(dataL[trl][:,0:nchan],axis=0)
        dataR[trl] = stats.zscore(dataR[trl][:,0:nchan],axis=0)
        # re-ref EEG
        ref_dataL=Zclass.avref(dataL[trl]) 
        ref_dataR=Zclass.avref(dataR[trl]) 
        ref_data=[ref_dataL, ref_dataR]
        # cut the longer eeg trial recording in the pair short (not optimal)
        nsampL = np.shape(ref_data[0])[0]
        nsampR = np.shape(ref_data[1])[0]
        nsamp_min=int(np.min((nsampL,nsampR)))        
        # decide length of epoch
        epoch=1 # 1 sec
        epoch = int(epoch*sr) # length of epoch
        # make the same number of nepoch for the pair
        nepoch = int(nsamp_min/epoch) 
        nsamp = nepoch*epoch  # length of trial data - cut the long one short (not good)
        # reshape into epoches x time x nchan
        ref_data[0] = np.reshape(ref_data[0][0:nsamp,0:nchan],(nepoch,epoch,nchan))
        ref_data[1] = np.reshape(ref_data[1][0:nsamp,0:nchan],(nepoch,epoch,nchan))
        # fft
        fcoef = dict() # create a dictionary to store the fcoef from the pair
        nbin = 50
        # detrend along the time dimension
        dataL = signal.detrend(ref_data[0],axis =1) # first subject L
        dataR = signal.detrend(ref_data[1],axis =1) # second subject R
        # fft along the time dimension
        datafL = fft(dataL,axis = 1)
        fcoef["L"] = datafL[:,0:nbin,:]/epoch
        datafR = fft(dataR,axis = 1)
        fcoef["R"] = datafR[:,0:nbin,:]/epoch
        del(dataL,dataR,ref_dataL,ref_dataR,ref_data,datafL,datafR) # delete the detrended data from workspace
        return fcoef

    def allspectra_z(fcoef,subj):
        nbin=50; nchan=32
        # L: subj=0; R: subj=L/R
        # compute power based on variance of fcoef
        # output power is 50 x 32 np ndarray
        power = np.var(fcoef[subj],axis = 0) # variance along epoch dimension (for L subject only)
        # initialize correlation of amplitue 
        ampcorr = np.zeros((nbin,nchan,nchan))
        # initialize cross spectrum (i.e. covariance) 
        cspectrum = np.zeros((nbin,nchan,nchan),dtype = complex)
        # initialize coherence
        coh = np.zeros((nbin,nchan,nchan))
        # loop over each freq
        for f in range(1,nbin): # one subject in each freq
            x_foef = np.transpose(np.squeeze(fcoef[subj][:,f,:])) # transpose into 32 chans x 128 epoch
            y_coh = np.abs(np.corrcoef(x_foef))**2 # coherence <- correlation**2 (32x32)
            z_cov = np.cov(x_foef) # covariance (32x32)
            a_ampcorr = np.corrcoef(np.abs(x_foef)) # amplitude correlation (32x32)
            ampcorr[f,:,:] = a_ampcorr
            coh[f,:,:] = y_coh
            cspectrum[f,:,:] = z_cov
        return power, ampcorr, cspectrum, coh

    def calparcoh(fcoef,subj,freq):
        # select the freq: e.g. freq=5Hz
        # create the object model
        cov_lasso = GraphicalLassoCV(cv=2)
        # compute amplitude in each epoch
        xx_amp= np.abs(np.squeeze(fcoef[subj][:,freq,:])) # select one frequency

        # # mean and std of amplitude across epochs
        # xx_amp_mean = np.mean(xx_amp,axis =  0) # 32 values 
        # # loop through each epoch to replace amplitue xx_amp as deviance
        # for j in range(np.shape(xx_amp)[0]):
        #     # compute the deviance of amplitude in each epoch
        #     xx_amp[j,:] = xx_amp[j,:] - xx_amp_mean # 148x32

        # round to the 6th decimal
        # xx_amp=np.round(xx_amp,6)

        # output in 148 epoch x 32 chan
        # Fit the GraphicalLasso covariance model to xx_amp (deviance).
        cov_lasso.fit(xx_amp)
        # return covariance and precision matrix
        cov_ = cov_lasso.covariance_ # covariance 32x32
        prec_ = cov_lasso.precision_ # precision 32x32

        # initialize the diagonal normalized precision - partial cohrelation matrix
        parcoh = np.zeros((32,32))
        # normalize each element in the precision matrix
        for j in range(32):
            for k in range(32):
                parcoh[j,k] = prec_[j,k]/np.sqrt(prec_[j,j]*prec_[k,k])
        return cov_, prec_, parcoh
    
    def skggmParcoh(fcoef,subj,freq):
        # compute amplitude in each epoch
        xx_amp= np.abs(np.squeeze(fcoef[subj][:,freq,:])) # select one frequency

        # create the object model
        modelquic = QuicGraphicalLassoCV(
        cv=2,  # cant deal w more folds at small size
        n_refinements=6,
        n_jobs=16,
        init_method="cov")

        modelquic.fit(xx_amp)

        print("   len(cv_lams): {}".format(len(modelquic.cv_lams_)))
        print("   lam_scale_: {}".format(modelquic.lam_scale_))
        print("   lam_: {}".format(modelquic.lam_))
        precision_quic = modelquic.precision_ # quic estimate cross validation
        covariance_quic = modelquic.covariance_

        # initialize the diagonal normalized precision - partial cohrelation matrix
        parcoh = np.zeros((32,32))
        # normalize each element in the precision matrix
        for j in range(32):
            for k in range(32):
                parcoh[j,k] = precision_quic[j,k]/np.sqrt(precision_quic[j,j]*precision_quic[k,k])
        return covariance_quic, precision_quic, parcoh


    def save_image(filename):
        # https://www.geeksforgeeks.org/save-multiple-matplotlib-figures-in-single-pdf-file-using-python/
        
        # PdfPages is a wrapper around pdf 
        # file so there is no clash and create
        # files with no error.
        p = PdfPages(filename)
        
        # get_fignums Return list of existing 
        # figure numbers
        fig_nums = plt.get_fignums()  
        figs = [plt.figure(n) for n in fig_nums]
        
        # iterating over the numbers in list
        for fig in figs: 
            
            # and saving the files
            fig.savefig(p, format='pdf') 
        
        # close the object
        p.close()  
