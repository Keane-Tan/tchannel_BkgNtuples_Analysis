import os
import time

os.system("tar -zcvf tchannelCom.tgz tchannelCom")
os.system("xrdcp -f tchannelCom.tgz root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannelCom2020/tchannelCom.tgz")

fileIDList = ["16 QCD", "16 TTJets", "16 WJets", "16 ZJets", "17 QCD", "17 TTJets", "17 WJets", "17 ZJets",
"18PRE QCD", "18PRE TTJets", "18PRE WJets", "18PRE ZJets", "18POST QCD", "18POST TTJets", "18POST WJets", "18POST ZJets"]

for i in range(len(fileIDList)): # changing Arguments in submit.jdl
	rfile = open('submit.jdl','r+')
	f1 = rfile.readlines()

	fileID = fileIDList[i]

	f1[10] = list(f1[10])[:-1]
	f1[10] = f1[10][:12] + list(fileID) + list("\n")
	f1[10] = ''.join(f1[10])
	rfile.seek(0)
	rfile.writelines(f1)
	rfile.truncate()
	rfile.close()

	os.system("condor_submit submit.jdl")
