import os,stat


overallSubmitFile = open("submitThresholdJobs.sh",'w+')
biases = [1,100,10000,100000]
clusterTypes = ["RHG","XZZX"]
dxList = range(2,15)
pList = [[.01,.015,.015,.015],[.01,.03,.05,.05]]

for biasIndex in range(len(biases)): 
    for clusterTypeIndex in range(len(clusterTypes)): 
        for dx in dxList:
            p = pList[clusterTypeIndex][biasIndex]
            clusterType = clusterTypes[clusterTypeIndex]
            eta = biases[biasIndex]
            fileName = "submitSimulateJob_"+str(eta)+"_"+str(dx)+"_"+clusterType+".sh"
            overallSubmitFile.write("sbatch "+fileName+"; ")
            submitFile = open(fileName,"w+")
            submitFile.write("#!/bin/bash\n")
            submitFile.write("#SBATCH --job-name="+fileName+"\n")
            submitFile.write("#SBATCH --time=24:00:00\n\n")
        
            submitFile.write("module load matplotlib/3.3.4-foss-2020b\n")
            submitFile.write("for i in {1..1000}; do python SimulateThresholds.py "+str(eta)+" "+str(p)+" "+str(dx)+" "+clusterType+" ; done")
            submitFile.close()
overallSubmitFile.close()
os.chmod("submitThresholdJobs.sh",stat.S_IRWXU)
