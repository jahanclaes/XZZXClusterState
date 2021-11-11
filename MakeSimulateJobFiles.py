import os,stat


overallSubmitFile = open("submitThresholdJobs.sh",'w+')
biases = [1,100,10000,100000]
clusterTypes = ["RHG","XZZX"]
dxList = range(2,16)
pList = [[.0008,.003,.003,.003],[.0005,.006,.01,.01]]

for biasIndex in [1,2]:#range(len(biases)): 
    for clusterTypeIndex in [1]:#range(len(clusterTypes)): 
        for dx in dxList[4:]:
            p = pList[clusterTypeIndex][biasIndex]
            clusterType = clusterTypes[clusterTypeIndex]
            eta = biases[biasIndex]
            fileName = "submitSimulateJob_"+str(eta)+"_"+str(dx)+"_"+clusterType+".sh"
            overallSubmitFile.write("sbatch "+fileName+"; ")
            submitFile = open(fileName,"w+")
            submitFile.write("#!/bin/bash\n")
            submitFile.write("#SBATCH --job-name=Cluster_"+str(eta)+"_"+str(dx)+"_"+clusterType+"\n")
            submitFile.write("#SBATCH --requeue\n")
            submitFile.write("#SBATCH --partition scavenge\n")
            submitFile.write("#SBATCH --time=24:00:00\n\n")
            if dx==15 and (clusterType=="RHG" or eta==1):
                submitFile.write("#SBATCH --mem=20G\n\n")
            if dx>8 and clusterType=="XZZX" and eta>1:
                submitFile.write("#SBATCH --mem=20G\n\n")
                
        
            submitFile.write("module load matplotlib/3.3.4-foss-2020b\n")
            submitFile.write("for i in {1..1000}; do python SimulateThresholds.py "+str(eta)+" "+str(p)+" "+str(dx)+" "+clusterType+" ; done")
            submitFile.close()
overallSubmitFile.close()
os.chmod("submitThresholdJobs.sh",stat.S_IRWXU)
