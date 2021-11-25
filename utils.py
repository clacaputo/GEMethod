## From spark_tnp by MuonPOG
import os
import uproot
import itertools


def get_pileup(resonance, era, subEra):
   '''
   Get the pileup distribution scalefactors to apply to simulation
   for a given era.
   '''
   # get the pileup
   dataPileup = {
       # Note: for now use ReReco version of pileup
       # TODO: need to redo splitting by 2016 B-F/F-H
       'Run2016_UL_HIPM': 'pileup/data/Run2016.root',
       'Run2016_UL': 'pileup/data/Run2016.root',
       'Run2017_UL': 'pileup/data/Run2017.root',
       'Run2018_UL': 'pileup/data/Run2018.root',
       # Double muon PD
       'Run2016_UL_HIPM_DM': 'pileup/data/Run2016.root',
       'Run2016_UL_DM': 'pileup/data/Run2016.root',
       'Run2017_UL_DM': 'pileup/data/Run2017.root',
       'Run2018_UL_DM': 'pileup/data/Run2018.root',
       'Run2016': 'pileup/data/Run2016.root',
       'Run2017': 'pileup/data/Run2017.root',
       'Run2018': 'pileup/data/Run2018.root'
   }
   mcPileup = {
       # TODO: do the two eras have different profiles?
       'Run2016_UL_HIPM': 'pileup/mc/Run2016_UL.root',
       'Run2016_UL': 'pileup/mc/Run2016_UL.root',
       'Run2017_UL': 'pileup/mc/Run2017_UL.root',
       'Run2018_UL': 'pileup/mc/Run2018_UL.root',
       # Double muon PD
       'Run2016_UL_HIPM_DM': 'pileup/mc/Run2016_UL.root',
       'Run2016_UL_DM': 'pileup/mc/Run2016_UL.root',
       'Run2017_UL_DM': 'pileup/mc/Run2017_UL.root',
       'Run2018_UL_DM': 'pileup/mc/Run2018_UL.root',
       'Run2016': 'pileup/mc/Run2016.root',
       'Run2017': 'pileup/mc/Run2017.root',
       'Run2018': 'pileup/mc/Run2018.root'
   }
   # get absolute path
   baseDir = os.path.dirname(__file__)
   dataPileup = {k: os.path.join(baseDir, dataPileup[k]) for k in dataPileup}
   mcPileup = {k: os.path.join(baseDir, mcPileup[k]) for k in mcPileup}
   with uproot.open(dataPileup[era]) as f:
       data_edges = f['pileup'].edges
       data_pileup = f['pileup'].values
       data_pileup /= sum(data_pileup)
   with uproot.open(mcPileup[era]) as f:
       mc_edges = f['pileup'].edges
       mc_pileup = f['pileup'].values
       mc_pileup /= sum(mc_pileup)
   pileup_edges = data_edges if len(data_edges) < len(mc_edges) else mc_edges
   pileup_ratio = [d/m if m else 1.0 for d, m in zip(
       data_pileup[:len(pileup_edges)-1], mc_pileup[:len(pileup_edges)-1])]

   return pileup_ratio, pileup_edges