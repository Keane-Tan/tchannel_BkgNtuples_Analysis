import os
import ROOT as rt
import numpy as np
import sys

rt.gROOT.LoadMacro("Interface/lester_mt2_bisect.h")

def vetoPhiSpike(etaHSL,phiHSL,rad,ai,bi):
	for iep in range(len(etaHSL)):
		if ((etaHSL[iep] - ai)**2 + (phiHSL[iep] - bi)**2 < rad):
			return True
			break

def ST_Val(METv,jets):
	evtST = METv
	for jet in jets:
		evtST += jet.Pt()
	return evtST

def deltaPhi(jetphiL,metPhiL):
	dphi = jetphiL - metPhiL
	if dphi < -np.pi:
		dphi += 2*np.pi
	if dphi > np.pi:
		dphi -= 2*np.pi
	return abs(dphi)

def deltaRforTL(jet1,jet2):
	phi1 = jet1.Phi()
	phi2 = jet2.Phi()
	eta1 = jet1.Eta()
	eta2 = jet2.Eta()
	dp = abs(phi1 - phi2)
	if dp > np.pi:
		dp -= 2*np.pi
	deltaR2 = (eta1 - eta2) * (eta1 - eta2) + dp * dp
	return np.sqrt(deltaR2)

def trans_mass_Njet(jet0, jet1, met, metPhi):
	visible = jet0 + jet1
	jetMass2 = visible.M2()
	term1 = rt.TMath.Sqrt(jetMass2 + visible.Pt()**2) * met
	term2 = rt.TMath.Cos(metPhi-visible.Phi())*visible.Pt()*met
	mT_first2 = rt.TMath.Sqrt(jetMass2 + 2*(term1 - term2))
	delR = deltaRforTL(jet0,jet1)
	delPhi = deltaPhi(jet0.Phi(),jet1.Phi())
	delEta = abs(jet0.Eta() - jet1.Eta())
	return [mT_first2,delR,delPhi,delEta]

def M_2J(j1,j2):
	totJets = rt.TLorentzVector(j1.Px(),j1.Py(),j1.Pz(),j1.Energy())
	totJets += rt.TLorentzVector(j2.Px(),j2.Py(),j2.Pz(),j2.Energy())
	return totJets.M()

def MPCom(comList):
	diffList = []
	for com in comList:
		jc1 = com[0]
		jc2 = com[1]
		jc3 = com[2]
		jc4 = com[3]
		diffList.append(M_2J(jc1,jc2) - M_2J(jc3,jc4))
	diffList = np.array(diffList)
	mdI = np.argmin(diffList)
	return comList[mdI]

def MT2DRCal(FDjet0,FSMjet0,FDjet1,FSMjet1,METx,METy):
	Fjet0 = FDjet0 + FSMjet0
	Fjet1 = FDjet1 + FSMjet1

	MT2v = rt.asymm_mt2_lester_bisect.get_mT2(
	Fjet0.M(), Fjet0.Px(), Fjet0.Py(),
	Fjet1.M(), Fjet1.Px(), Fjet1.Py(),
	METx, METy, 0.0, 0.0, 0
	)
	delR = deltaRforTL(Fjet0,Fjet1)
	delPhi = deltaPhi(Fjet0.Phi(),Fjet1.Phi())
	delEta = abs(Fjet0.Eta() - Fjet1.Eta())

	return [MT2v,delR,delPhi,delEta]

# years = ["16","17","18PRE","18POST"]
# bkgs = ["QCD","TTJets","WJets","ZJets"]
years = [sys.argv[1]]
bkgs = [sys.argv[2]]

ffile = open('fileIDList.txt','r+')
fileIDList = ffile.readlines()
baseloc = "root://cmseos.fnal.gov//store/user/keanet/CondorOutput/tchannel/BkgNtuples/MCRoot/"

# check and make sure certain folders exist
foldList = []
for bkg in bkgs:
	for year in years:
		fold = "20" + year + "/" + bkg + year
		foldList.append(fold)

for fold in foldList:
	yearPart = fold.find("/")
	var_year = fold[yearPart+1:]

	print var_year

	listOfFiles = []

	for fileID in fileIDList:
		fileID = fileID[:-1] # gets rid of the "\n" character

		if var_year in fileID:
			listOfFiles.append(baseloc + fileID + ".root")


	# working with the tree of the combined root files
	msd_max = 200
	Bkg_Pt = rt.TH1F(var_year + "_Pt",var_year + "_Pt;p_{T}(j) [GeV];Events",40, 0, 3500)
	Bkg_Eta = rt.TH1F(var_year + "_Eta",var_year + "_Eta;#eta(j);Events",15, -7, 7)
	Bkg_Phi = rt.TH1F(var_year + "_Phi",var_year + "_Phi;#phi(j);Events",15, -3.2, 3.2)
	Bkg_DPhi = rt.TH1F(var_year + "_DPhi",var_year + "_DPhi; #Delta#phi(j,#slash{E}_{T});Events",30,0,3.5)
	Bkg_M = rt.TH1F(var_year + "_M",var_year + "_M; M(j);Events",40,0,1200)
	Bkg_axisMajor = rt.TH1F(var_year + "_axisMajor",var_year + "_axisMajor;#sigma_{major}(j);Events",40, 0, 0.5)
	Bkg_axisMinor = rt.TH1F(var_year + "_axisMinor",var_year + "_axisMinor;#sigma_{minor}(j);Events",40, 0, 0.3)
	Bkg_momentGirth = rt.TH1F(var_year + "_momentGirth",var_year + "_momentGirth;girth(j);Events",40, 0, 0.5)
	Bkg_ptD = rt.TH1F(var_year + "_ptD",var_year + "_ptD;p_{T}D(j);Events",40, 0, 1.2)
	Bkg_tau1 = rt.TH1F(var_year + "_tau1",var_year + "_tau1;#tau_{1}(j);Events",40, 0, 0.8)
	Bkg_tau2 = rt.TH1F(var_year + "_tau2",var_year + "_tau2;#tau_{2}(j);Events",40, 0, 0.65)
	Bkg_tau3 = rt.TH1F(var_year + "_tau3",var_year + "_tau3;#tau_{3}(j);Events",40, 0, 0.35)
	Bkg_tau21 = rt.TH1F(var_year + "_tau21",var_year + "_tau21;#tau_{21}(j);Events",40, 0, 1.3)
	Bkg_tau32 = rt.TH1F(var_year + "_tau32",var_year + "_tau32;#tau_{32}(j);Events",40, 0, 1.3)
	Bkg_SDM = rt.TH1F(var_year + "_SDM",var_year + "_SDM; m_{SD}(j);Events",40,0,msd_max)
	Bkg_pmult = rt.TH1F(var_year + "_pmult",var_year + "_pmult;Number of Particles in a Jet;Events",40,0,400)
	# ST variables
	Bkg_ST = rt.TH1F(var_year + "_ST",var_year + "_ST;S_{T} [GeV];Events",50, 0, 10000)
	Bkg_MET = rt.TH1F(var_year + "_MET",var_year + "_MET;#slash{E}_{T} [GeV];Events",50, 0, 10000)
	# Mass variables
	## MT
	Bkg_MT_2j = rt.TH1F(var_year + "_MT_2j",var_year + "_MT_2j;m_{T} (GeV);Events",40, 0, 9000)
	Bkg_MT_3j = rt.TH1F(var_year + "_MT_3j",var_year + "_MT_3j;m_{T} (GeV);Events",40, 0, 9000)
	Bkg_MT_4j = rt.TH1F(var_year + "_MT_4j",var_year + "_MT_4j;m_{T} (GeV);Events",40, 0, 9000)
	Bkg_MT_4pj = rt.TH1F(var_year + "_MT_4pj",var_year + "_MT_4pj;m_{T} (GeV);Events",40, 0, 9000)
	## MT2 : 4 highest pT jets, most probable combination? The 2 pairs of jets are similar in their mass?
	Bkg_MT2_mp = rt.TH1F(var_year + "_MT2_mp",var_year + "_MT2_mp; M_{T2} (GeV);Events",40, 0, 9000)


	tr = rt.TChain("bkg_tree")
	for fileName in listOfFiles:
		tr.Add(fileName)

	nEvents = tr.GetEntries()

	for count in range(nEvents):
		if count % 1000 == 0:
			print "Event " + str(count)
		tr.GetEntry(count)
		nAK8Jets = tr.nAK8Jets
		AK8Jets = tr.AK8Jets
		MET = tr.MET
		METPhi = tr.METPhi
		pmult = tr.pmult
		girth = tr.girth
		axisMj = tr.axisMj
		axisMn = tr.axisMn
		ptD = tr.ptD
		tau1 = tr.tau1
		tau2 = tr.tau2
		tau3 = tr.tau3
		SDM = tr.SDM
		jetID = tr.jetID
		dPhiMin = tr.dPhiMin
		NElectrons = tr.NElectrons
		NMuons = tr.NMuons
		Muons_MiniIso = tr.Muons_MiniIso
		MT_AK8 = tr.MT_AK8
		weightedLumi = tr.weightedLumi

		if len(AK8Jets) < 2:
			continue

		if AK8Jets[0].Pt() < 200 or AK8Jets[1].Pt() < 200:
			continue

		# phi spike filters


		rad = 0.028816 # half the length of the diagonal of the eta-phi rectangular cell
		rad = rad * 0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study

		veto = 0

		# hot spots for leading jets
		eta16lead = [0.048,0.24,1.488,1.584,-1.008]
		phi16lead = [-0.35,-0.35,-0.77,-0.77,-1.61]

		eta17lead = [0.144,1.488,1.488,1.584,-0.624]
		phi17lead = [-0.35,-0.77,-0.63,-0.77,0.91]

		eta18lead = [1.488,1.488,1.584]
		phi18lead = [-0.77,-0.63,-0.77]

		# hot spots for subleading jets
		eta16 = [-1.2,-0.912,-0.912,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,
0.912,0.912,1.008,1.296,-1.584,-0.816,-0.72,-0.144,-0.048,-0.048,0.048,1.104,1.488]
		phi16 = [-1.19,2.03,3.01,-1.75,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,
-0.49,0.63,1.47,-2.31,0.07,-2.59,0.77,0.91,-3.15,2.73]

		eta17 = [-0.912,-0.912,-0.816,-0.72,-0.528,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,
-1.2,-0.72,-0.72,-0.432,0.336,0.624,1.104,1.296]
		phi17 = [2.03,3.01,-1.75,-0.77,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-1.19,-2.31,-2.17,
2.73,-0.77,-0.77,-3.15,-0.49]

		eta18PRE = [-1.584,-1.2,-0.912,-0.912,-0.816,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,
-0.144,-0.048,-0.048,0.144,0.912,0.912,1.008,1.296,-0.72,1.104,1.488,1.776]
		phi18PRE = [0.63,-1.19,2.03,3.01,-1.75,-0.77,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.07,0.77,0.77,0.91,
1.75,1.75,2.87,0.63,-0.49,-2.31,-3.15,-0.21,0.77]

		eta18POST = [-1.2,-0.912,-0.912,-0.816,-0.72,-0.528,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,
1.008,1.296,-1.584,-0.816,-0.72,-0.72,-0.432,-0.144,-0.048,1.104,1.488,1.776]
		phi18POST = [-1.19,2.03,3.01,-1.75,-0.77,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-0.49,0.63,
-0.77,-2.31,-2.17,2.73,0.07,0.77,-3.15,-0.21,0.77]

		if "16" in fileID:
			if vetoPhiSpike(eta16,phi16,rad,AK8Jets[1].Eta(),AK8Jets[1].Phi()) or vetoPhiSpike(eta16lead,phi16lead,rad,AK8Jets[0].Eta(),AK8Jets[0].Phi()):
				veto = 1

		if "17" in fileID:
			if vetoPhiSpike(eta17,phi17,rad,AK8Jets[1].Eta(),AK8Jets[1].Phi()) or vetoPhiSpike(eta17lead,phi17lead,rad,AK8Jets[0].Eta(),AK8Jets[0].Phi()):
				veto = 1

		if "18PRE" in fileID:
			if vetoPhiSpike(eta18PRE,phi18PRE,rad,AK8Jets[1].Eta(),AK8Jets[1].Phi()) or vetoPhiSpike(eta18lead,phi18lead,rad,AK8Jets[0].Eta(),AK8Jets[0].Phi()):
				veto = 1

		if "18POST" in fileID:
			if vetoPhiSpike(eta18POST,phi18POST,rad,AK8Jets[1].Eta(),AK8Jets[1].Phi()) or vetoPhiSpike(eta18lead,phi18lead,rad,AK8Jets[0].Eta(),AK8Jets[0].Phi()):
				veto = 1

		if veto == 1:
#			vcount += 1
			continue

		# Calculating the ST variable
		Bkg_MET.Fill(MET,weightedLumi)
		Bkg_ST.Fill(ST_Val(MET,AK8Jets),weightedLumi)

		if len(AK8Jets) >= 2:
			MTv = trans_mass_Njet(AK8Jets[0],AK8Jets[1],MET,METPhi)
			if len(AK8Jets) == 2 and AK8Jets[0].Pt() > 200 and AK8Jets[1].Pt() > 200:
				Bkg_MT_2j.Fill(MTv[0],weightedLumi)
			if len(AK8Jets) == 3 and AK8Jets[0].Pt() > 200 and AK8Jets[1].Pt() > 200:
				Bkg_MT_3j.Fill(MTv[0],weightedLumi)
			if len(AK8Jets) == 4 and AK8Jets[0].Pt() > 200 and AK8Jets[1].Pt() > 200:
				Bkg_MT_4j.Fill(MTv[0],weightedLumi)
			if len(AK8Jets) > 4 and AK8Jets[0].Pt() > 200 and AK8Jets[1].Pt() > 200:
				Bkg_MT_4pj.Fill(MTv[0],weightedLumi)

		if len(AK8Jets) >= 4:
			j1 = AK8Jets[0]
			j2 = AK8Jets[1]
			j3 = AK8Jets[2]
			j4 = AK8Jets[3]

			com1 = [j1,j2,j3,j4]
			com2 = [j1,j3,j2,j4]
			com3 = [j1,j4,j2,j3]

			mpCom = MPCom([com1,com2,com3])
			MT2 = MT2DRCal(mpCom[0],mpCom[1],mpCom[2],mpCom[3],MET*np.cos(METPhi),MET*np.sin(METPhi))
			Bkg_MT2_mp.Fill(MT2[0],weightedLumi)

		for ij in range(len(AK8Jets)):
			Bkg_Pt.Fill(AK8Jets[ij].Pt(),weightedLumi)
			Bkg_DPhi.Fill(deltaPhi(AK8Jets[ij].Phi(),METPhi),weightedLumi)
			Bkg_M.Fill(AK8Jets[ij].M(),weightedLumi)
			Bkg_axisMajor.Fill(axisMj[ij],weightedLumi)
			Bkg_axisMinor.Fill(axisMn[ij],weightedLumi)
			Bkg_momentGirth.Fill(girth[ij],weightedLumi)
			Bkg_ptD.Fill(ptD[ij],weightedLumi)
			Bkg_tau1.Fill(tau1[ij],weightedLumi)
			Bkg_tau2.Fill(tau2[ij],weightedLumi)
			Bkg_tau3.Fill(tau3[ij],weightedLumi)
			if tau1[ij] > 0:
				Bkg_tau21.Fill(tau2[ij]/tau1[ij],weightedLumi)
			if tau2[ij] > 0:
				Bkg_tau32.Fill(tau3[ij]/tau2[ij],weightedLumi)
			Bkg_pmult.Fill(pmult[ij],weightedLumi)

		for sji in range(len(SDM)):
			Bkg_SDM.Fill(SDM[sji],weightedLumi)

	olist = [Bkg_Pt,Bkg_Eta,Bkg_Phi,Bkg_DPhi,Bkg_M,Bkg_axisMajor,Bkg_axisMinor,
	Bkg_momentGirth,Bkg_ptD,Bkg_tau1,Bkg_tau2,Bkg_tau3,Bkg_tau21,Bkg_tau32,
	Bkg_SDM,Bkg_pmult,Bkg_ST,Bkg_MET,Bkg_MT_2j,Bkg_MT_3j,Bkg_MT_4j,Bkg_MT_4pj,
	Bkg_MT2_mp]

	rootOutFile = rt.TFile.Open("tCom_" + var_year + ".root", "recreate")
	rootOutFile.cd()
	for histo in olist:
		histo.Write()
	rootOutFile.Close()
