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

#calculate packing fraction average and error
PF_data=np.loadtxt('vol_frac.txt')
PF=PF_data[:,1]
PF_avg=np.mean(PF)

#estimate standard deviation of averages
PF_corr_time=compute_corr_time(PF)
PF_samples=np.size(PF,axis=0)/PF_corr_time
PF_stdev=np.std(PF)
PF_error=PF_stdev/np.sqrt(PF_samples)

print "Packing fraction average:", PF_avg,",","Number of samples:", PF_samples,",","Standard deviation:", PF_stdev,",","Error", PF_error
