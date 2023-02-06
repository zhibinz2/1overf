
#%% import modules
from os import listdir

from hdf5storage import loadmat, savemat 
import numpy as np 
from matplotlib import pyplot as plt

from scipy import signal, stats
from scipy.fftpack import fft

#%% 
states4names=['Uncoupled','Leading','Following','Mutual']
from zzpackage.zzmodule import Zclass
dir='Cleaned_data/'

# %% file loading sequence
# organize conditions in all sessions into a vector of 192
filedates=[20220713,20220721,20220804,20220808,20220810,20220811,20220815,20220816,20221003,2022100401,2022100402,20221005]
numSes=len(filedates)

#%% Organize indicies of the conditions sequentially
# conditions_all=[]
# from zzpackage.zzmodule import Zclass
# for i in range(len(filedates)): # in time sequence    
#     filename='clean_'+str(filedates[i])+'.mat' 
#     _, _,_,_,conditions,_,_,_,_,_,_,_,sr = Zclass.loadmatvar2(dir,filename) 
#     conditions_all.append(conditions)

# %% save conditions_all and load it
# outdict=dict()
# outdict['conditions_all']=conditions_all
# savemat('conditions_all',outdict,store_python_metadata=True)
outdict=loadmat('conditions_all.mat')
conditions_all=outdict["conditions_all"]

# %% Organize indicies for the 4 states
Uncoupled_Ind=np.zeros((12,3))
L_Lead_Ind=np.zeros((12,3)) # they were leading indicies for L subject but following for R
R_Lead_Ind=np.zeros((12,3))
Mutual_Ind=np.zeros((12,3))
for ses in range(12):
    Uncoupled_Ind[ses]=np.asarray(np.where(conditions_all[ses]==1))
    L_Lead_Ind[ses]=np.asarray(np.where(conditions_all[ses]==2))
    R_Lead_Ind[ses]=np.asarray(np.where(conditions_all[ses]==3))
    Mutual_Ind[ses]=np.asarray(np.where(conditions_all[ses]==4))

# %% load pacorr_mat
# pacorr_all=np.zeros((12,2,12,30,32,32))
# for ses in range(12):
# 	filename='clean_'+str(filedates[ses])+'_pcorr.mat'
# 	outdict=loadmat('pcorr/'+filename)
# 	partial_correlation=outdict['partial_correlation']
# 	for subj in range(2):
# 		for trl in range(12):			
# 			for freq in range(30):
# 				pacorr_all[ses][subj][trl][freq,:,:]=partial_correlation[subj][trl][freq,:,:]

#%% save pacorr_all and load it
# outdict=dict()
# outdict['pacorr_all']=pacorr_all
# savemat('pacorr_all.mat',outdict,store_python_metadata=True)			
dict=loadmat('pacorr_all.mat')
pacorr_all=dict['pacorr_all']

# %% extract non zero matricies by thresholding
# nonzerosmat=np.zeros((12,2,12,30,32,32))
# for ses in range(12):
# 	for subj in range(2):
# 		for trl in range(12):			
# 			for freq in range(30):
# 				nonzerosmat[ses][subj][trl][freq,:,:]\
# 					=(pacorr_all[ses][subj][trl][freq,:,:]>0.1)*1\
#                         +(pacorr_all[ses][subj][trl][freq,:,:]<-0.1)*-1

# %% no need for thresholding
nonzerosmat=pacorr_all
			
# %% append the same state together
# append all the Uncoupled_Ind nozeromat in the same state
append_Uncoupled=list(); append_Mutual=list()
# append_L_Lead=list(); append_R_Lead=list(); 
append_Leading=list(); append_Following=list()
np.append
for ses in range(12):
    for subj in range(2):
        for trl in range(3):
            if subj == 0:
                append_Uncoupled.append(nonzerosmat[ses][subj][int(Uncoupled_Ind[ses][trl])])
                append_Leading.append(nonzerosmat[ses][subj][int(L_Lead_Ind[ses][trl])])
                append_Following.append(nonzerosmat[ses][subj][int(R_Lead_Ind[ses][trl])])
                append_Mutual.append(nonzerosmat[ses][subj][int(Mutual_Ind[ses][trl])])
            else:
                append_Uncoupled.append(nonzerosmat[ses][subj][int(Uncoupled_Ind[ses][trl])])
                append_Leading.append(nonzerosmat[ses][subj][int(R_Lead_Ind[ses][trl])])
                append_Following.append(nonzerosmat[ses][subj][int(L_Lead_Ind[ses][trl])])
                append_Mutual.append(nonzerosmat[ses][subj][int(Mutual_Ind[ses][trl])])

    
# %%
# Correct looping, loop only number of trl
# sum all the Uncoupled_Ind nonzeromat
sum_Uncoupled=np.zeros((30,32,32)); sum_Mutual=np.zeros((30,32,32))
sum_Leading=np.zeros((30,32,32));sum_Following=np.zeros((30,32,32))
for trl in range(72):
    sum_Uncoupled[:,:,:]=sum_Uncoupled[:,:,:]+append_Uncoupled[trl][:,:,:]
    sum_Leading[:,:,:]=sum_Leading[:,:,:]+append_Leading[trl][:,:,:]
    sum_Following[:,:,:]=sum_Following[:,:,:]+append_Following[trl][:,:,:]
    sum_Mutual[:,:,:]=sum_Mutual[:,:,:]+append_Mutual[trl][:,:,:]


# %%  names for the 4 states  
sum_4states=[sum_Uncoupled,sum_Leading,sum_Following,sum_Mutual]

# %% plot the figures and save in pdf
# in matplotlib.rcParams
plt.rcParams['figure.figsize']=[30, 100]

# visualization of parcor selection
figure,axs = plt.subplots(30,4,constrained_layout=True)
for state in range(4):
	for freq in range(30):
		im=axs[freq,state].imshow(sum_4states[state][freq,:,:],\
					vmin=-2,vmax=2,cmap='jet')#RdBu_r viridis
		axs[freq,state].set_title(states4names[state]+' freq: '+str(freq+1)+' Hz')
		figure.colorbar(im, ax=axs[freq,state])

# name your Pdf file
paffilename = "state4sum.pdf"
# call the function
from zzpackage.zzmodule import Zclass
Zclass.save_image(paffilename)


#%% save sum of nonzeros in table and save as txt file
del(dict)
sum_4states_list=dict()
for state in range(4):
    tmp=np.zeros((30))
    for freq in range(30):    
        tmp[freq]=np.sum(sum_4states[state][freq,:,:])-(72*32)
    sum_4states_list[states4names[state]]=tmp


# %%
from prettytable import PrettyTable

columns=['Frequency (Hz)']
columns.extend(states4names)

myTable = PrettyTable()

# Add Columns
myTable.add_column(columns[0], np.arange(1,31,1).tolist())
col=1
for key,value in sum_4states_list.items():
    value=value.tolist()
    myTable.add_column(columns[col], value)
    col=col+1


print(myTable)

# print to txt file
with open('state4sum.txt', 'w') as w:
    w.write(str(myTable))


#%% plot verus freq
plt.rcParams['figure.figsize']=[3, 3]
for state,color in enumerate(['green','red','blue','black']):
    plt.plot(range(30),sum_4states_list[states4names[state]],color)
    plt.xlim(1,30);plt.ylim(-1250,-600)
    plt.xlabel("frequency (Hz)");plt.ylabel("non zero sum")
plt.legend(states4names)
plt.show()

