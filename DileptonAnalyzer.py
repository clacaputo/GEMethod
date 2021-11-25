import uproot
import numpy as np
import math
from utils import get_pileup

from ROOT import TH1F, TH2F, TFile
from root_numpy import fill_hist

def deltaPhi(p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while (res > math.pi):
        res -= 2*math.pi
    while (res < -math.pi):
        res += 2*math.pi
    return res

filename = "/eos/cms/store/group/phys_muon/moh/TnP_ntuples/muon/Z/Run2018_UL/MINIAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_TnP_ntuplizer_muon_Z_Run2018_UL_MINIAOD_DY_madgraph/210625_155340/0000/output_46.root"

config = {"selection": "HLT_Mu50_v==1 and tag_TuneP_pt>53 and tag_abseta<2.4 and tag_isHighPt==1 and tag_momentum_unc and tag_track_isolated and tag_dxy_abs < 0.02 and probe_TuneP_pt > 25 and probe_isHighPt==1 and probe_momentum_unc and probe_track_isolated and probe_dxy_abs < 0.02",
            "definitions": {"probe_abseta": "abs(probe_eta)",
             "tag_abseta": "abs(tag_eta)",
             "tag_momentum_unc": "tag_TuneP_pterr < (0.3*tag_TuneP_pt)",
             "probe_momentum_unc": "probe_TuneP_pterr < (0.3*probe_TuneP_pt)",
             "tag_track_isolated": "tag_iso03_sumPt < (0.1*tag_TuneP_pt)",
             "probe_track_isolated": "probe_iso03_sumPt < (0.1*probe_TuneP_pt)",
             "tag_kappa": "tag_charge/(tag_TuneP_pt/1000)",
             "probe_kappa": "probe_charge/(probe_TuneP_pt/1000)",
             "tag_dxy_abs": "abs(tag_dxy)",
             "probe_dxy_abs": "abs(probe_dxy)"},
        }

dPhi = np.frompyfunc(deltaPhi,2,1)

with uproot.open(filename) as f:
    print(f.keys())
    print(f['muon'].keys())
    #print(f['muon']['Events'].keys())
    tree = f['muon']['Events'].arrays(namedecode="utf-8")

    outputfile = uproot.recreate("tmp.root", compression=uproot.ZLIB(4))
    rootoutputfile = TFile("tmp_root.root","recreate")
    rootoutputfile.cd()

    pileup_ratio, pileup_edges = get_pileup("Z", "Run2018_UL","")
    pileupMap = {e: r for e, r in zip(pileup_edges[:-1], pileup_ratio)}
    #print("PU: \n", pileupMap)
    #print(pileupMap[np.around(tree["nTrueInteractions"][1])])
    puMap = lambda x: pileupMap[x]
    puMapV = np.vectorize(puMap)
    puWeights = puMapV(np.around(tree["nTrueInteractions"]))
    #print(puMapV(np.around(tree["nTrueInteractions"]))[:10])
    #for newvar, formula in config['definitions'].items():
    #    print(newvar, formula)
    #    exec(f"{newvar} = ")
    
    HTL_Mu50_path = tree["HLT_Mu50_v"]

    tag_eta = tree['tag_eta']
    tag_phi = tree['tag_phi']
    tag_charge = tree["tag_charge"]
    tag_tuneP_pterr = tree["tag_tuneP_pterr"]
    tag_tuneP_pt    = tree["tag_tuneP_pt"]
    tag_iso03_sumPt = tree["tag_iso03_sumPt"]
    tag_dxy = tree["tag_dxy"]
    tag_isHighPt = tree["tag_isHighPt"]
    
    probe_eta = tree['probe_eta']
    probe_phi = tree['probe_phi']
    probe_charge = tree["probe_charge"]
    probe_tuneP_pterr = tree["probe_tuneP_pterr"]
    probe_tuneP_pt    = tree["probe_tuneP_pt"]
    probe_iso03_sumPt = tree["probe_iso03_sumPt"]
    probe_dxy = tree["probe_dxy"]
    probe_isHighPt = tree["probe_isHighPt"]


    tag_abseta =  np.abs( tag_eta )
    tag_momentum_unc  =  tag_tuneP_pterr   < ( 0.3*tag_tuneP_pt )
    probe_momentum_unc=  probe_tuneP_pterr < ( 0.3*probe_tuneP_pt )
    tag_track_isolated=  tag_iso03_sumPt   < ( 0.1*tag_tuneP_pt )
    probe_track_isolated  =  probe_iso03_sumPt < (0.1*probe_tuneP_pt)
    tag_kappa =  tag_charge/(tag_tuneP_pt/1000)
    probe_kappa   =  probe_charge/(probe_tuneP_pt/1000)
    tag_dxy_abs   =  np.abs( tag_dxy )
    probe_dxy_abs =  np.abs( probe_dxy )


    triggerMask    = ( HTL_Mu50_path == 1)
    leadingLepMask = ( (tag_tuneP_pt>53) & (tag_abseta<2.4) & 
                       (tag_isHighPt==1) & (tag_momentum_unc) & 
                       (tag_track_isolated) & (tag_dxy_abs < 0.02) 
                      )
    subleadingLepMask = ( (probe_tuneP_pt > 25) & (probe_isHighPt==1) & 
                          (probe_momentum_unc) & (probe_track_isolated) & 
                          (probe_dxy_abs < 0.02) )

    selectionMask = ( triggerMask & leadingLepMask & subleadingLepMask)

    zMass2 = np.sqrt(2*tag_tuneP_pt*probe_tuneP_pt*
                       ( np.cosh( tag_eta-probe_eta ) - 
                         np.cos( dPhi(tag_phi,probe_phi).astype(float) )
                        )
                    
                    )

    ## Masks for differential analysis
    zMass2_highMask = (zMass2 > 200.)
    zMass2_LowMask = (zMass2 > 55.)
    zMass2_peakMask = ((zMass2 < 120.) & (zMass2 > 65.))

    tag_eta = tag_eta [selectionMask]
    probe_eta = probe_eta[selectionMask]
    tag_phi = tag_phi[selectionMask]
    probe_phi = probe_phi[selectionMask]


    etaRegions_lambdas_dict = {
    # Eta regions
    "eta_m2p4_m2p1" : lambda x: ((x > -2.4) & (x <= -2.1)),
    "eta_m2p1_m1p2" : lambda x: ((x > -2.1) & (x <= -1.2)),
    "eta_m1p2_0p"   : lambda x: ((x > - 1.2) & (x <= 0)),
    "eta_0p_1p2"    : lambda x: ((x > 0) & (x <= 1.2)),
    "eta_1p2_2p1"   : lambda x: ((x > 1.2) & (x <= 2.1)),
    "eta_2p1_2p4"   : lambda x: ((x > 2.1) & (x <= 2.4))
    }
    phiRegions_lambdas_dict = {
    # Phi regions
    "phi_mPI_mPI3"  : lambda x: ((x > -math.pi) & (x <= -math.pi/3.)),
    "phi_mPI3_PI3"  : lambda x: ((x > -math.pi/3.) & (x <= math.pi/3.)),
    "phi_PI3_PI"    : lambda x: ((x > math.pi/3.) & (x <= math.pi))
    }

    #mu2_etaRegions_dict = {
    ## Eta regions
    #"eta_m2p4_m2p1" : ((probe_eta > -2.4) & (probe_eta <= -2.1)),
    #"eta_m2p1_m1p2" : ((probe_eta > -2.1) & (probe_eta <= -1.2)),
    #"eta_m1p2_0p"   : ((probe_eta > - 1.2) & (probe_eta <= 0)),
    #"eta_0p_1p2"    : ((probe_eta > 0) & (probe_eta <= 1.2)),
    #"eta_1p2_2p1"   : ((probe_eta > 1.2) & (probe_eta <= 2.1)),
    #"eta_2p1_2p4"   : ((probe_eta > 2.1) & (probe_eta <= 2.4))
    #}
    #mu2_phiRegions_dict = {
    ## Phi regions
    #"phi_mPI_mPI3"  : ((probe_phi > -math.pi) & (probe_phi <= -math.pi/3.)),
    #"phi_mPI3_PI3"  : ((probe_phi > -math.pi/3.) & (probe_phi <= math.pi/3.)),
    #"phi_PI3_PI"    : ((probe_phi > math.pi/3.) & (probe_phi <= math.pi))
    #}

    InjScaleVector = np.arange(-1.0,1.001, 0.01)
    InjScaleMap = {}
    bin_kappa = np.array(np.linspace(-20, 20, 100))

    ## Weights for MC
    # Lumi Weigth
    #PUwei = cms.GetPuWeight(genInfo.trueNumberOfInteractions);
	#weight = weight*cms.GetPuWeight(genInfo.trueNumberOfInteractions);
	#weight = weight*genInfo.MCweight;
    #
    #// ID Iso
	#cms.Apply_lepton_SF(mu1, weight);
	#cms.Apply_lepton_SF(mu2, weight);
	#
	#// trigger
	#cms.Apply_dimuon_triggerSF(mu1, mu2, weight);
    #cms.Apply_lepton_PtZ(pt_Z, weight);

    for scaleBias in InjScaleVector:
        tmp_tag_kappa   = tag_kappa[selectionMask] + scaleBias
        tmp_probe_kappa = probe_kappa[selectionMask] + scaleBias
        puWeights_sel = puWeights[selectionMask]
        for etaReg in etaRegions_lambdas_dict.keys():
            for phiReg in phiRegions_lambdas_dict.keys():
        #InjScaleMap[f"{i:1.2f}"] = tmp_tag_kappa
                etaRegionSelector = etaRegions_lambdas_dict[etaReg]
                phiRegionSelector = phiRegions_lambdas_dict[phiReg]
                mu1_eta_phi_mask = (etaRegionSelector(tag_eta) & phiRegionSelector(tag_phi))
                mu2_eta_phi_mask = (etaRegionSelector(probe_eta) & phiRegionSelector(probe_phi))
                #print( len(mu1_etaRegions_dict[etaReg]) , len(mu1_phiRegions_dict[phiReg]) )
                #print( len(mu2_etaRegions_dict[etaReg]) , len(mu2_phiRegions_dict[phiReg]) )
                #print( len(mu1_eta_phi_mask) , len(mu2_eta_phi_mask) )
                differential_tag_kappa = tmp_tag_kappa[mu1_eta_phi_mask]
                differential_probe_kappa = tmp_probe_kappa[mu2_eta_phi_mask]
                #print(f"hKMuon_{scaleBias:1.2f}_{etaReg}_{phiReg}")
                histname = f"hKMuon_{scaleBias:1.2f}_{etaReg}_{phiReg}"
                histname_weigthed = f"hKMuon_{scaleBias:1.2f}_{etaReg}_{phiReg}_w"
                tmp_kappa = np.concatenate((differential_tag_kappa, differential_probe_kappa))

                hist = TH1F(histname, 'title', 100, -20, 20)
                histW = TH1F(histname_weigthed, 'title', 100, -20, 20)
                fill_hist(histW, tmp_kappa, weights=np.concatenate((puWeights_sel[mu1_eta_phi_mask],
                                                                    puWeights_sel[mu2_eta_phi_mask])))
                fill_hist(hist, tmp_kappa)
                outputfile[histname] = np.histogram( tmp_kappa, bins=bin_kappa)
                hist.Write()
                histW.Write()
                del hist
                del histW
    
    rootoutputfile.Close()
    outputfile["hKMuon"] = np.histogram( np.concatenate( ( tag_kappa[(selectionMask)], probe_kappa[(selectionMask)]) ), 
                                                        bins=bin_kappa )
    outputfile["hKMuon_puW"] = np.histogram( np.concatenate( ( tag_kappa[(selectionMask)], probe_kappa[(selectionMask)]) ), 
                                                        bins=bin_kappa , 
                                                        weights=np.concatenate( (puWeights[selectionMask], puWeights[selectionMask])) )
    outputfile["hKMuon_zMass2_LowMask"] = np.histogram( np.concatenate( (
                                                                        tag_kappa[((selectionMask) & (zMass2_LowMask))],
                                                                        probe_kappa[((selectionMask) & (zMass2_LowMask))]
                                                                        ) 
                                                                       ), 
                                                        bins=bin_kappa
                                                        )
    outputfile["hKMuon_zMass2_highMask"] = np.histogram( np.concatenate( (
                                                                        tag_kappa[((selectionMask) & (zMass2_highMask))],
                                                                        probe_kappa[((selectionMask) & (zMass2_highMask))]
                                                                        ) 
                                                                       ), 
                                                        bins=bin_kappa
                                                        )
    outputfile["hKMuon_zMass2_peakMask"] = np.histogram( np.concatenate( (
                                                                        tag_kappa[((selectionMask) & (zMass2_peakMask))],
                                                                        probe_kappa[((selectionMask) & (zMass2_peakMask))]
                                                                        ) 
                                                                       ), 
                                                        bins=bin_kappa
                                                        )