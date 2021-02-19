#ifndef Physics_Analysis_KsKpi_H
#define Physics_Analysis_KsKpi_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"


class KsKpi : public Algorithm {

public:
  KsKpi(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:

  double m_Ecms;
  // Declare r0, z0 cut for charged tracks
  double m_vr0cut;
  double m_vz0cut;
  double m_mdccos;

  //Declare energy, dphi, dthe cuts for fake gamma's
  double m_energyThreshold1;
  double m_energyThreshold2;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaAngleCut;

  int m_test4C;
  int m_test5C;
  int m_checkDedx;
  int m_checkTof;

  // define Ntuples here
  NTuple::Tuple*  m_tuple1;
  NTuple::Matrix<double>  mcp3_true;
  NTuple::Array<double>  mcp_true;
  NTuple::Array<double>  mce_true;
  NTuple::Array<double>  mcm_true;
  NTuple::Item<double>  mc_cos;
  NTuple::Item<double>  mc_pt;
  NTuple::Item<double>  mc_fpvr;
  NTuple::Item<double>  mc_fpvz;
  NTuple::Array<double> m_Rz;
  NTuple::Array<double> m_Rxy;
  NTuple::Matrix<double>  m_vfit_p3;
  NTuple::Array<double>   m_vfit_p;
  NTuple::Array<double>   m_vfit_e;
  NTuple::Item<double>   m_vfit_m;
  NTuple::Item<double>   m_vfit_cos;
  NTuple::Item<double>   m_vfit_phi;
  NTuple::Item<double>   m_vfit_pt;
  NTuple::Item<double>   m_vfit_pkaon;
  NTuple::Array<double>   m_vfit_mkstar;
  NTuple::Item<double>   m_var_kstar;
  NTuple::Item<double>  m_vfits_chi;
  NTuple::Item<double>  m_vfit2_mks;
  NTuple::Item<double>  m_vfit2_chi;
  NTuple::Item<double>  m_vfit2_ct;
  NTuple::Item<double>  m_vfit2_dl;
  NTuple::Item<double>  m_vfit2_dle;
  NTuple::Item<int>  m_mode;
  NTuple::Item<long>  m_nGood;
  NTuple::Item<long>  m_nkcond;
  NTuple::Array<double> m_angle;
  NTuple::Item<long>  m_nMatch;
  NTuple::Array<double> m_angle_match;
  NTuple::Array<double> m_pid_dedx_K;
  NTuple::Array<double> m_pid_tof1_K;
  NTuple::Array<double> m_pid_tof2_K;
  NTuple::Array<double> m_pid_prob_K;
  NTuple::Array<double> nM_tag;
  NTuple::Array<double> m_charge;
  NTuple::Array<double> m_krvxy0;
  NTuple::Array<double> m_krvz0;
  NTuple::Array<double> nG_tag;
  NTuple::Item<long>  m_pid_absolute;
  NTuple::Item<long>  m_pid_relative;
  NTuple::Item<long>  m_pid_relative_b;
  NTuple::Item<long>  m_pid_misid_pi;
  NTuple::Item<long>  m_pid_misid_proton;
  NTuple::Item<long>  m_pid_dedx;
  NTuple::Item<long>  m_pid_tof;
  NTuple::Item<long>  m_pid_tof_e;
  NTuple::Item<double>  m_final_jpsi;
  NTuple::Item<double>  m_final_kaonp;
  NTuple::Item<double>  m_final_kaonphi;
  NTuple::Item<double>  m_final_kaoncos;
  NTuple::Item<long> m_dedx_pid_relative;
  NTuple::Item<long> m_tof_pid_relative;
  NTuple::Item<int>  m_nrun;
  NTuple::Item<int>  m_nrec;
  NTuple::Item<int>   m_idxmc;
  NTuple::Array<int>  m_pdgid;
  NTuple::Array<int>  m_motheridx;
  NTuple::Array<int>  m_drank;
  NTuple::Array<int>  m_motherid;
  NTuple::Array<long>  m_trkidx;
};

#endif 
