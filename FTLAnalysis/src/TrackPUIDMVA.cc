#include "PrecisionTiming/FTLAnalysis/interface/TrackPUIDMVA.h"

TrackPUIDMVA::TrackPUIDMVA(std::string weights_file, bool is4D):
    is4D_(is4D)
{
    //---3D
    vars_.push_back(std::make_tuple("pt", 0.));
    vars_.push_back(std::make_tuple("eta", 0.));
    vars_.push_back(std::make_tuple("phi", 0.));
    // vars_.push_back(std::make_tuple("dxy", 0.));
    // vars_.push_back(std::make_tuple("dz", 0.));
    // vars_.push_back(std::make_tuple("dzErr", 0.));
    // vars_.push_back(std::make_tuple("dxyErr", 0.));
    vars_.push_back(std::make_tuple("chi2", 0.));
    vars_.push_back(std::make_tuple("ndof", 0.));
    vars_.push_back(std::make_tuple("numberOfValidHits", 0.));
    vars_.push_back(std::make_tuple("numberOfValidPixelBarrelHits", 0.));
    vars_.push_back(std::make_tuple("numberOfValidPixelEndcapHits", 0.));        
    //---4D
    if(is4D)
    {
        //vars_.push_back(std::make_tuple("dt", 0.));
        //vars_.push_back(std::make_tuple("sigmat0", 0.));
        vars_.push_back(std::make_tuple("btlMatchChi2", 0.));
        vars_.push_back(std::make_tuple("btlMatchTimeChi2", 0.));        
        vars_.push_back(std::make_tuple("etlMatchChi2", 0.));
        vars_.push_back(std::make_tuple("etlMatchTimeChi2", 0.));        
        vars_.push_back(std::make_tuple("mtdt", 0.));
        vars_.push_back(std::make_tuple("path_len", 0.));
    }

    mva_ = MVAComputer(&vars_, weights_file);
}

float TrackPUIDMVA::operator() (reco::TrackBaseRef& trk, reco::Vertex& vtx)
{
    const auto& pattern = trk->hitPattern();                
    
    std::get<1>(vars_[0]) = trk->pt();
    std::get<1>(vars_[1]) = trk->eta();
    std::get<1>(vars_[2]) = trk->phi();
    std::get<1>(vars_[3]) = trk->dxy(vtx.position());
    std::get<1>(vars_[4]) = trk->dz(vtx.position());                                     
    std::get<1>(vars_[5]) = trk->dzError();
    std::get<1>(vars_[6]) = trk->dxyError();
    std::get<1>(vars_[7]) = trk->chi2();
    std::get<1>(vars_[8]) = trk->ndof();
    std::get<1>(vars_[9]) = trk->numberOfValidHits();
    std::get<1>(vars_[10]) = pattern.numberOfValidPixelBarrelHits();
    std::get<1>(vars_[11]) = pattern.numberOfValidPixelEndcapHits();

    return mva_();
}

float TrackPUIDMVA::operator() (reco::TrackBaseRef& trk, reco::TrackBaseRef& ext_trk, reco::Vertex& vtx,
                                edm::ValueMap<float>& t0s,
                                edm::ValueMap<float>& sigma_t0s,
                                edm::ValueMap<float>& btl_chi2s,
                                edm::ValueMap<float>& btl_time_chi2s,                      
                                edm::ValueMap<float>& etl_chi2s,
                                edm::ValueMap<float>& etl_time_chi2s,
                                edm::ValueMap<float>& tmtds,                                
                                edm::ValueMap<float>& trk_lengths)
{
    const auto& pattern = ext_trk->hitPattern();                
    
    std::get<1>(vars_[0]) = trk->pt();
    std::get<1>(vars_[1]) = trk->eta();
    std::get<1>(vars_[2]) = trk->phi();
    std::get<1>(vars_[3]) = trk->chi2();
    std::get<1>(vars_[4]) = trk->ndof();
    std::get<1>(vars_[5]) = trk->numberOfValidHits();
    std::get<1>(vars_[6]) = pattern.numberOfValidPixelBarrelHits();
    std::get<1>(vars_[7]) = pattern.numberOfValidPixelEndcapHits();
    std::get<1>(vars_[8]) = btl_chi2s.contains(ext_trk.id()) ? btl_chi2s[ext_trk] : -1;
    std::get<1>(vars_[9]) = btl_time_chi2s.contains(ext_trk.id()) ? btl_time_chi2s[ext_trk] : -1;    
    std::get<1>(vars_[10]) = etl_chi2s.contains(ext_trk.id()) ? etl_chi2s[ext_trk] : -1;
    std::get<1>(vars_[11]) = etl_time_chi2s.contains(ext_trk.id()) ? etl_time_chi2s[ext_trk] : -1;    
    std::get<1>(vars_[12]) = tmtds[ext_trk];
    std::get<1>(vars_[13]) = trk_lengths[ext_trk];

    // std::get<1>(vars_[0]) = trk->pt();
    // std::get<1>(vars_[1]) = trk->eta();
    // std::get<1>(vars_[2]) = trk->phi();
    // std::get<1>(vars_[3]) = trk->dxy(vtx.position());
    // std::get<1>(vars_[4]) = trk->dz(vtx.position());                                     
    // std::get<1>(vars_[5]) = trk->dzError();
    // std::get<1>(vars_[6]) = trk->dxyError();
    // std::get<1>(vars_[7]) = trk->chi2();
    // std::get<1>(vars_[8]) = trk->ndof();
    // std::get<1>(vars_[9]) = trk->numberOfValidHits();
    // std::get<1>(vars_[10]) = pattern.numberOfValidPixelBarrelHits();
    // std::get<1>(vars_[11]) = pattern.numberOfValidPixelEndcapHits();
    // std::get<1>(vars_[12]) = t0s.contains(trk.id()) ? t0s[trk]-vtx.t() : -1-vtx.t();
    // std::get<1>(vars_[13]) = sigma_t0s.contains(trk.id()) ? sigma_t0s[trk] : -1;    
    // std::get<1>(vars_[14]) = btl_chi2s.contains(ext_trk.id()) ? btl_chi2s[ext_trk] : -1;
    // std::get<1>(vars_[15]) = btl_time_chi2s.contains(ext_trk.id()) ? btl_time_chi2s[ext_trk] : -1;    
    // std::get<1>(vars_[16]) = etl_chi2s.contains(ext_trk.id()) ? etl_chi2s[ext_trk] : -1;
    // std::get<1>(vars_[17]) = etl_time_chi2s.contains(ext_trk.id()) ? etl_time_chi2s[ext_trk] : -1;    
    // std::get<1>(vars_[18]) = tmtds[ext_trk];
    // std::get<1>(vars_[19]) = trk_lengths[ext_trk];

    return mva_();
}
