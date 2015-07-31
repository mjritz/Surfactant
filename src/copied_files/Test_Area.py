import numpy as np

def compute_corr_time(data,block_sizes=np.array([100,300,500,700,1000],dtype=int)):
#Block averaging of data to determine 2*corr time
        n_samples=int(np.trunc(np.size(data,axis=0)/np.max(block_sizes)))*np.max(block_sizes)
        data_var=np.var(data[0:n_samples])
        data_var_b=np.zeros_like(block_sizes,dtype=float)
        data_var_t=np.zeros_like(block_sizes,dtype=float)
        for j in range(np.size(block_sizes)):
                block=block_sizes[j]
                for i in range(n_samples/block):
                        data_var_b[j]+=np.var(data[i*block:(i+1)*block])
                data_var_b[j]=data_var_b[j]/(n_samples/block)
                data_var_t[j]=data_var/block
        p=np.polyfit(data_var_t,data_var_b,1)
	return abs(p[0])

#set constants
kT_kcal=0.0019872041*298.15
kT_mJ=1.3806488*(10**-20)*298.15
#dA in meters squared
dA=2*(60.0**2)*((10**-10)**2)*0.0005


#calculate gamma by test area method
dE_data=np.loadtxt('dE.txt')
Up=np.exp(-(dE_data[:,2]-dE_data[:,1])/kT_kcal)
Um=np.exp(-(dE_data[:,3]-dE_data[:,1])/kT_kcal)
Up_avg=np.mean(Up)
Um_avg=np.mean(Um)
gamma=-kT_mJ*(np.log(Up_avg)-np.log(Um_avg))/(2*dA)
#print 'gamma (mJ/m^2): ',gamma

#estimate standard deviation of averages
Up_corr_time=compute_corr_time(Up)
Up_samples=np.size(Up,axis=0)/Up_corr_time
Up_stdev=np.std(Up)
Up_error=Up_stdev/np.sqrt(Up_samples)

Um_corr_time=compute_corr_time(Um)
Um_samples=np.size(Um,axis=0)/Um_corr_time
Um_stdev=np.std(Um)
Um_error=1.0*Um_stdev/np.sqrt(Um_samples)

print Up_error, Um_error, Up_corr_time, Um_corr_time

	
error=(kT_mJ/(2*dA))*np.sqrt((Up_error**2)/((Up_avg-Um_avg)**2)+(Um_error**2)/((Up_avg-Um_avg)**2))

print 'gamma (mJ/m^2): ',gamma,' +/- ',error
