import ROOT as R
import math as M
from operator import itemgetter

f=R.TFile("/tmp/DumpHits_sm10_2.root")

dh=f.Get("DumpHits")

histos = {}
histos["track_pt"]=R.TH1F("track_pt","track_pt",100,0.,10.)
histos["track_eta"]=R.TH1F("track_eta","track_eta",100,0.,1.5)
histos["track_phi"]=R.TH1F("track_phi","track_phi",100,-M.pi,M.pi)
histos["track_path_vs_pt"]=R.TProfile("track_path_vs_pt","track_path_vs_pt",100,0.,10.)
histos["matchedTrack_pt"]=R.TH1F("matchedTrack_pt","matchedTrack_pt",100,0.,10.)
histos["matchedTrack_eta"]=R.TH1F("matchedTrack_eta","matchedTrack_eta",100,0.,1.5)
histos["matchedTrack_phi"]=R.TH1F("matchedTrack_phi","matchedTrack_phi",100,-M.pi,M.pi)
histos["matchedTrack_nCluster"]=R.TH1F("matchedTrack_nCluster","matchedTrack_nCluster",10,0.,10.)
histos["matchedTrack_nHits"]=R.TH1F("matchedTrack_nHits","matchedTrack_nHits",10,0.,10.)
histos["cluster_energy"]=R.TH1F("cluster_energy","cluster_energy",100,0.,20.)
histos["cluster_time"]=R.TH1F("cluster_time","cluster_time",100,0.,25.)
histos["cluster_size"]=R.TH1F("cluster_size","cluster_size",20,0.,20.)
histos["cluster_sizeX"]=R.TH1F("cluster_sizeX","cluster_sizeX",20,0.,20.)
histos["cluster_sizeY"]=R.TH1F("cluster_sizeY","cluster_sizeY",20,0.,20.)
histos["cluster_seedEnergyRatio"]=R.TH1F("cluster_seedEnergyRatio","cluster_seedEnergyRatio",110,0.,1.1)
histos["matchedCluster_energy"]=R.TH1F("matchedCluster_energy","matchedCluster_energy",100,0.,20.)
histos["matchedCluster_time"]=R.TH1F("matchedCluster_time","matchedCluster_time",100,0.,25.)
histos["matchedCluster_DR"]=R.TH1F("matchedCluster_DR","matchedCluster_DR",100,0.,0.05)
histos["matchedCluster_DEta"]=R.TH1F("matchedCluster_DEta","matchedCluster_DEta",100,-0.05,0.05)
histos["matchedCluster_DPhi"]=R.TH1F("matchedCluster_DPhi","matchedCluster_DPhi",100,-0.05,0.05)
histos["matchedCluster_size"]=R.TH1F("matchedCluster_size","matchedCluster_size",20,0.,20.)
histos["matchedCluster_size_simHits"]=R.TH2F("matchedCluster_size_simHits","matchedCluster_size_simHits",20,0.,20.,20,0.,20.)
histos["matchedCluster_sizeX"]=R.TH1F("matchedCluster_sizeX","matchedCluster_sizeX",20,0.,20.)
histos["matchedCluster_sizeY"]=R.TH1F("matchedCluster_sizeY","matchedCluster_sizeY",20,0.,20.)
histos["matchedCluster_size_vs_pt"]=R.TH2F("matchedCluster_size_vs_pt","matchedCluster_size_vs_pt",100,0.,10.,20,-0.5,19.5)
histos["matchedCluster_size_vs_eta"]=R.TH2F("matchedCluster_size_vs_eta","matchedCluster_size_vs_eta",100,0.,1.5,20,-0.5,19.5)
histos["matchedMultipleCluster_map"]=R.TH2F("matchedMultipleCluster_map","matchedCluster_map",100,0.,1.5,100,-M.pi,M.pi)
histos["matchedMultipleCluster_local_x"]=R.TH1F("matchedMultipleCluster_local_x","matchedMultipleCluster_local_x",1000,-10,10)
histos["matchedMultipleCluster_local_y"]=R.TH1F("matchedMultipleCluster_local_y","matchedMultipleCluster_local_y",1000,-10,10)
histos["matchedMultipleCluster_local_z"]=R.TH1F("matchedMultipleCluster_local_z","matchedMultipleCluster_local_z",1000,-10,10)
histos["matchedMultipleCluster_row"]=R.TH1F("matchedMultipleCluster_row","matchedMultipleCluster_row",100,0,100)
histos["matchedMultipleCluster_col"]=R.TH1F("matchedMultipleCluster_col","matchedMultipleCluster_col",100,0,100)
histos["matchedTrack_local_x"]=R.TH1F("matchedTrack_local_x","matchedTrack_local_x",1000,-10,10)
histos["matchedTrack_local_y"]=R.TH1F("matchedTrack_local_y","matchedTrack_local_y",1000,-10,10)
histos["matchedTrack_local_z"]=R.TH1F("matchedTrack_local_z","matchedTrack_local_z",1000,-10,10)
histos["matchedTrack_row"]=R.TH1F("matchedTrack_row","matchedTrack_row",100,0,100)
histos["matchedTrack_col"]=R.TH1F("matchedTrack_col","matchedTrack_col",100,0,100)
histos["hitMap"]=R.TProfile2D("hitMap","hitMap",3,-1.5,1.5,3,-1.5,1.5)

for event in dh:
    for iclus in range(0,event.clusters_n):
        histos["cluster_energy"].Fill(event.clusters_energy[iclus])
        histos["cluster_time"].Fill(event.clusters_time[iclus])
        histos["cluster_size"].Fill(event.clusters_size[iclus])
        histos["cluster_sizeX"].Fill(event.clusters_size_x[iclus])
        histos["cluster_sizeY"].Fill(event.clusters_size_y[iclus])
        histos["cluster_seedEnergyRatio"].Fill(event.clusters_seed_energy[iclus]/event.clusters_energy[iclus])
    for itrack in range(0,len(event.track_idx)):
        histos["track_pt"].Fill(event.track_pt[itrack])
        histos["track_eta"].Fill(abs(event.track_eta_atBTL[itrack]))
        histos["track_phi"].Fill(event.track_phi_atBTL[itrack])
        histos["matchedTrack_nCluster"].Fill(event.matchedClusters_n[itrack])
        histos["matchedTrack_nHits"].Fill(event.matchedRecHits_n[itrack])
        firstSimHit=-1
        minR=99999.
        #try to find the first simHit
        hitMap = {}
        nhitMap = {}
        totEnergy=0
        #fill empty map
        for ihit in range(0,event.matchedSimHits_n[itrack]):
            hitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]=0
            nhitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]=0
        for ihit in range(0,event.matchedSimHits_n[itrack]):
            hitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]=event.matchedSimHits_energy[itrack][ihit]+hitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]
            nhitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]=1+nhitMap["%d_%d"%(event.matchedSimHits_ieta[itrack][ihit],event.matchedSimHits_iphi[itrack][ihit])]
            totEnergy=totEnergy+event.matchedSimHits_energy[itrack][ihit]
        sortedHitMap=sorted(hitMap.items(), key=itemgetter(1),reverse=True)
        if (len(sortedHitMap)>0):
            seed_ieta=int(sortedHitMap[0][0].split("_")[0])
            seed_iphi=int(sortedHitMap[0][0].split("_")[1])
            for (key,val) in sortedHitMap:
                ieta=int(key.split("_")[0])
                iphi=int(key.split("_")[1])
                histos["hitMap"].Fill(ieta-seed_ieta,iphi-seed_iphi,val/totEnergy)
        for ihit in range(0,event.matchedSimHits_n[itrack]):
            if (event.matchedSimHits_entry_global_R[itrack][ihit]<minR):
                firstSimHit=ihit
                minR=event.matchedSimHits_entry_global_R[itrack][ihit]
        if (firstSimHit!=-1):
            histos["matchedTrack_local_x"].Fill(event.matchedSimHits_entry_local_x[itrack][firstSimHit])
            histos["matchedTrack_local_y"].Fill(event.matchedSimHits_entry_local_y[itrack][firstSimHit])
            histos["matchedTrack_local_z"].Fill(event.matchedSimHits_entry_local_z[itrack][firstSimHit])
            histos["matchedTrack_row"].Fill(event.matchedSimHits_row[itrack][firstSimHit])
            histos["matchedTrack_col"].Fill(event.matchedSimHits_col[itrack][firstSimHit])
        if (event.matchedClusters_n[itrack]>0):
            histos["matchedTrack_pt"].Fill(event.track_pt[itrack])
            histos["matchedTrack_eta"].Fill(abs(event.track_eta_atBTL[itrack]))
            histos["matchedTrack_phi"].Fill(event.track_phi_atBTL[itrack])
        if (event.matchedClusters_n[itrack]>1):
            histos["matchedMultipleCluster_map"].Fill(abs(event.track_eta_atBTL[itrack]),event.track_phi_atBTL[itrack])
            if (firstSimHit!=-1):
                histos["matchedMultipleCluster_local_x"].Fill(event.matchedSimHits_entry_local_x[itrack][firstSimHit])
                histos["matchedMultipleCluster_local_y"].Fill(event.matchedSimHits_entry_local_y[itrack][firstSimHit])
                histos["matchedMultipleCluster_local_z"].Fill(event.matchedSimHits_entry_local_z[itrack][firstSimHit])
                histos["matchedMultipleCluster_row"].Fill(event.matchedSimHits_row[itrack][firstSimHit])
                histos["matchedMultipleCluster_col"].Fill(event.matchedSimHits_col[itrack][firstSimHit])

        for iclus in range(0,event.matchedClusters_n[itrack]):
            histos["matchedCluster_energy"].Fill(event.matchedClusters_energy[itrack][iclus])
            histos["matchedCluster_time"].Fill(event.matchedClusters_time[itrack][iclus])
            histos["matchedCluster_DR"].Fill(event.matchedClusters_track_DR[itrack][iclus])
            histos["matchedCluster_DEta"].Fill(event.matchedClusters_track_Deta[itrack][iclus])
            histos["matchedCluster_DPhi"].Fill(event.matchedClusters_track_Dphi[itrack][iclus])
            histos["matchedCluster_size"].Fill(event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_simHits"].Fill(len(sortedHitMap),event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_pt"].Fill(event.track_pt[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_eta"].Fill(event.track_eta_atBTL[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_sizeX"].Fill(event.matchedClusters_size_x[itrack][iclus])
            histos["matchedCluster_sizeY"].Fill(event.matchedClusters_size_y[itrack][iclus])

histos["effCluster_pt"]=R.TGraphAsymmErrors(histos["matchedTrack_pt"],histos["track_pt"])
histos["effCluster_eta"]=R.TGraphAsymmErrors(histos["matchedTrack_eta"],histos["track_eta"])
histos["effCluster_phi"]=R.TGraphAsymmErrors(histos["matchedTrack_phi"],histos["track_phi"])
histos["effMultipleCluster_local_x"]=R.TGraphAsymmErrors(histos["matchedMultipleCluster_local_x"],histos["matchedTrack_local_x"])
histos["effMultipleCluster_local_y"]=R.TGraphAsymmErrors(histos["matchedMultipleCluster_local_y"],histos["matchedTrack_local_y"])
histos["effMultipleCluster_row"]=R.TGraphAsymmErrors(histos["matchedMultipleCluster_row"],histos["matchedTrack_row"])
histos["effMultipleCluster_col"]=R.TGraphAsymmErrors(histos["matchedMultipleCluster_col"],histos["matchedTrack_col"])

fOut=R.TFile("clusterPlots_sm10_2.root","RECREATE")
for hn, histo in histos.iteritems():
    if isinstance(histo,R.TH1F):
        histo.SetMinimum(0.)
    histo.Write()
fOut.Close()
