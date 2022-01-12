import os,stat


biases = [1,5,10,50,100,500,1000,5000,10000,10001]
clusterTypes = ["RHG","XZZX"]
pList = [[.016,.026,.027,.03,.03,.03,.03,.03,.03,.03],[.016,.036,.043,.059,.062,.075,.08,.081,.082,.1]]

overallSubmitFile = open("submitThresholdJobs.sh",'w+')
overallSubmitFile2D = open("submitThresholdJobs2D.sh",'w+')

for biasIndex in range(4,len(biases)): 
    for clusterTypeIndex in range(len(clusterTypes)): 
        if clusterTypeIndex == 1 and biasIndex>3:
            dxList = range(4,9)
        else:
            dxList = range(12,16)
        for dx in dxList:
            p = pList[clusterTypeIndex][biasIndex]
            clusterType = clusterTypes[clusterTypeIndex]
            eta = biases[biasIndex]
            fileName = "submitSimulateJob_Z_"+str(eta)+"_"+str(dx)+"_"+clusterType+".sh"
            submitFile = open(fileName,"w+")
            submitFile.write("#!/bin/bash\n")
            submitFile.write("#SBATCH --job-name=Cluster_Z_"+str(eta)+"_"+str(dx)+"_"+clusterType+"\n")
            submitFile.write("#SBATCH --requeue\n")
            submitFile.write("#SBATCH --partition scavenge\n")
            submitFile.write("#SBATCH --time=24:00:00\n\n")
            if dx==15 and (clusterType=="RHG" or eta==1):
                submitFile.write("#SBATCH --mem=20G\n\n")
            if dx>8 and clusterType=="XZZX" and eta>1:
                submitFile.write("#SBATCH --mem=20G\n\n")
            submitFile.write("module load matplotlib/3.3.4-foss-2020b\n")
            if clusterType=="RHG":
                submitFile.write("for i in {1..1000}; do python SimulateThresholds.py "+str(eta)+" "+str(p)+" "+str(dx)+" "+clusterType+" ; done")
            else:
                submitFile.write("for i in {1..1000}; do python SimulateThresholds.py "+str(eta)+" "+str(p)+" "+str(dx)+" ; done")
            submitFile.close()
            if eta<=10000:
                overallSubmitFile.write("sbatch "+fileName+"; ")
            else:
                overallSubmitFile2D.write("sbatch "+fileName+"; ")

overallSubmitFile.close()
os.chmod("submitThresholdJobs.sh",stat.S_IRWXU)
os.chmod("submitThresholdJobs2D.sh",stat.S_IRWXU)
