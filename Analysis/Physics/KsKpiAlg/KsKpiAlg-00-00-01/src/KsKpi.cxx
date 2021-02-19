#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "KsKpiAlg/KsKpi.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"   
#include "ParticleID/ParticleID.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "McTruth/McParticle.h"
#include <vector>
/*************************************************************************/
/**************** variable defination*************************************/
/*************************************************************************/
const double mpion = 0.13957;
const double mkaon=0.493677;
const double mks0 = 0.497614;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int igg=0;
int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,NcutX,NcutX1;
KsKpi::KsKpi(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
		declareProperty("Ecms",   m_Ecms=3.0971873);
		declareProperty("Vr0cut", m_vr0cut=1.0);
		declareProperty("Vz0cut", m_vz0cut=10.0);
		declareProperty("mdccos", m_mdccos=0.93);
		declareProperty("EnergyThreshold1", m_energyThreshold1=0.025);
		declareProperty("EnergyThreshold2", m_energyThreshold2=0.050);
		declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
		declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
		declareProperty("GammaAngleCut", m_gammaAngleCut=20.0);
	}
/*************************************************************************/
/***************root book*************************************************/
/*************************************************************************/
StatusCode KsKpi::initialize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;

	NTuplePtr nt1(ntupleSvc(), "FILE1/vfit");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/vfit", CLID_ColumnWiseTuple, 
				"ks N-Tuple example");
		if ( m_tuple1 ) {  
			status = m_tuple1->addItem("vfitp3",              5, 3, m_vfit_p3); 
			status = m_tuple1->addItem("vfitp",               5, m_vfit_p);
			status = m_tuple1->addItem("vfite",               5, m_vfit_e);
			status = m_tuple1->addItem("vfitm",               m_vfit_m);
			status = m_tuple1->addItem("vfitpkaon",           m_vfit_pkaon);
			status = m_tuple1->addItem("vfitcos",             m_vfit_cos);
			status = m_tuple1->addItem("vfitphi",             m_vfit_phi);
			status = m_tuple1->addItem("vfitpt",              m_vfit_pt);
			status = m_tuple1->addItem("var_kstar",           m_var_kstar);
			status = m_tuple1->addItem("vfitmkstar",          2, m_vfit_mkstar);
			status = m_tuple1->addItem("vfits_chi",           m_vfits_chi); 
			status = m_tuple1->addItem("vfit2_mks",           m_vfit2_mks);
			status = m_tuple1->addItem("vfit2_chi",           m_vfit2_chi);
			status = m_tuple1->addItem("vfit2_ct",            m_vfit2_ct);
			status = m_tuple1->addItem("vfit2_dl",            m_vfit2_dl);
			status = m_tuple1->addItem("vfit2_dle",           m_vfit2_dle);
			status = m_tuple1->addItem("mode",                m_mode);
			status = m_tuple1->addItem("nGood",               m_nGood, 0, 30);
			status = m_tuple1->addItem("nkcond",              m_nkcond, 0, 30);
			status = m_tuple1->addIndexedItem("mangle",       m_nkcond, m_angle);
			status = m_tuple1->addIndexedItem("KRz",          m_nkcond, m_Rz);
			status = m_tuple1->addIndexedItem("KRxy",         m_nkcond, m_Rxy);
			status = m_tuple1->addIndexedItem("nG_tag",       m_nkcond, nG_tag);
			status = m_tuple1->addItem("nMatch",              m_nMatch, 0, 30);
			status = m_tuple1->addIndexedItem("mangle_match", m_nMatch, m_angle_match);
			status = m_tuple1->addIndexedItem("nM_tag",       m_nMatch, nM_tag);
			status = m_tuple1->addIndexedItem("pid_dedx_K",   m_nMatch, m_pid_dedx_K);
			status = m_tuple1->addIndexedItem("pid_tof1_K",   m_nMatch, m_pid_tof1_K);
			status = m_tuple1->addIndexedItem("pid_tof2_K",   m_nMatch, m_pid_tof2_K);
			status = m_tuple1->addIndexedItem("pid_prob_K",   m_nMatch, m_pid_prob_K);
			status = m_tuple1->addItem("pid_absolute",        m_pid_absolute);
			status = m_tuple1->addItem("pid_relative",        m_pid_relative);
			status = m_tuple1->addItem("pid_relative_b",      m_pid_relative_b);
			status = m_tuple1->addItem("pid_misid_pi",        m_pid_misid_pi);
			status = m_tuple1->addItem("pid_misid_proton",    m_pid_misid_proton);
			status = m_tuple1->addItem("pid_dedx",            m_pid_dedx);
			status = m_tuple1->addItem("pid_tof",             m_pid_tof);
			status = m_tuple1->addItem("pid_tof_e",           m_pid_tof_e);
			status = m_tuple1->addItem("final_jpsi",          m_final_jpsi);
			status = m_tuple1->addItem("final_kaonp",         m_final_kaonp);
			status = m_tuple1->addItem("final_kaonphi",       m_final_kaonphi);
			status = m_tuple1->addItem("final_kaoncos",       m_final_kaoncos);

			status = m_tuple1->addItem("nrun",             m_nrun);
			status = m_tuple1->addItem("nrec",             m_nrec);
			status = m_tuple1->addItem("mcp",              6, mcp_true);
			status = m_tuple1->addItem("mcp3",             6, 3, mcp3_true);
			status = m_tuple1->addItem("mce",              6, mce_true);
			status = m_tuple1->addItem("mcm",              6, mcm_true);
			status = m_tuple1->addItem("mccos",            mc_cos);
			status = m_tuple1->addItem("mcpt",             mc_pt);
			status = m_tuple1->addItem("indexmc",          m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("drank",     m_idxmc, m_drank);
			status = m_tuple1->addIndexedItem("trackidx",  m_idxmc, m_trkidx);
			status = m_tuple1->addIndexedItem("motherpid", m_idxmc, m_motherid);
			status = m_tuple1->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
		}
		else    { 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) 
				<< endmsg;
			return StatusCode::FAILURE;
		}
	}

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;
}

/*************************************************************************/
/*****************Main program********************************************/
/*************************************************************************/
StatusCode KsKpi::execute() {
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),
			"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	m_nrun=runNo;
	m_nrec= event;
	log << MSG::DEBUG <<"run, evtnum = "
		<< runNo << " , "
		<< event <<endreq;
	Ncut0++;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::
			EvtRecEvent);
	log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::
			EvtRecTrackCol);
	SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
			"/Event/MC/McParticleCol");
	/***********************************************************************/
	/********************mctrue information*********************************/
	/***********************************************************************/
	double temp_p_kp,temp_p_km,temp_p_pp,temp_p_pm;
	if(m_nrun<0) {
		if(!mcParticleCol){
			cout<<"Could not retrieve McParticelCol"<<endl;
			return( StatusCode::FAILURE);
		}
		int temp_ngam=0;
		double temp_ppi0;

		for(int imc=0;imc<6;imc++){
			mcp_true[imc]=0;
			mce_true[imc]=0;
		}

		Event::McParticle temp;

		bool JpsiDecay = false;
		int rootIndex = -1;
		HepLorentzVector p_true,p_pp,p_pm;
		HepLorentzVector mc_p_temp;
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++){
			HepLorentzVector true_4mom=(*iter_mc)->initialFourMomentum();
			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;

			int pdgid = (*iter_mc)->particleProperty();
			int pdgid_mother = ((*iter_mc)->mother()).particleProperty(); 

			if(pdgid==321){
				mc_p_temp = true_4mom;

				mcp3_true[1][0]=mc_p_temp.px();
				mcp3_true[1][1]=mc_p_temp.py();
				mcp3_true[1][2]=mc_p_temp.pz();
				mce_true[1]=mc_p_temp.e();

				mc_cos=mc_p_temp.cosTheta();
				mc_pt=sqrt(mc_p_temp.px()*mc_p_temp.px()+mc_p_temp.py()*mc_p_temp.py());
			}
			if(pdgid==310){
				mc_p_temp = true_4mom;

				mcp3_true[2][0]=mc_p_temp.px();
				mcp3_true[2][1]=mc_p_temp.py();
				mcp3_true[2][2]=mc_p_temp.pz();
				mce_true[2]=mc_p_temp.e();
			}
			if(pdgid==-211&&pdgid_mother!=310){

				mc_p_temp = true_4mom;

				mcp3_true[3][0]=mc_p_temp.px();
				mcp3_true[3][1]=mc_p_temp.py();
				mcp3_true[3][2]=mc_p_temp.pz();
				mce_true[3]=mc_p_temp.e();
			}	 	  
			if(pdgid==-211&&pdgid_mother==310){      
				mc_p_temp = true_4mom;
				mcp3_true[4][0]=mc_p_temp.px();
				mcp3_true[4][1]=mc_p_temp.py();
				mcp3_true[4][2]=mc_p_temp.pz();
				mce_true[4]=mc_p_temp.e();
			}	 	  

			if(pdgid==211&&pdgid_mother==310){      
				mc_p_temp = true_4mom;

				mcp3_true[5][0]=mc_p_temp.px();
				mcp3_true[5][1]=mc_p_temp.py();
				mcp3_true[5][2]=mc_p_temp.pz();
				mce_true[5]=mc_p_temp.e();
			}
		}// end of mc particle loop
	}

	/***********************************************************************/
	/***********************Good track selection****************************/
	/***********************************************************************/
	Vint iGood,ipionp,ipionm;

	iGood.clear();
	ipionp.clear();
	ipionm.clear();

	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		if(fabs(cos(mdcTrk->theta()))>m_mdccos) continue;
		iGood.push_back(i);
	}

	int nGood = iGood.size();
	log << MSG::DEBUG << "ngood" << nGood << endreq;
	if(nGood <3) return StatusCode::SUCCESS;
	m_nGood=nGood;
	Ncut1++;
	/***********************************************************************/
	/*********************PID***********************************************/
	/***********************************************************************/
	int npionp=0;
	int npionm=0;
	int ippionm=9999;
	int ippionp=9999;

	HepLorentzVector ppionp, ppionm;
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i]; 
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() );
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() );   
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		if( ( pid->probPion() > pid->probKaon()) && (pid->probPion() > pid->probProton())) { 
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
			if (!(*itTrk)->isMdcKalTrackValid()) continue;
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
			if(mdcKalTrk->charge()>0){
				npionp++;
				ippionp = iGood[i]; 
				ipionp.push_back(ippionp);
			}
			else{
				npionm++;
				ippionm = iGood[i];
				ipionm.push_back(ippionm);
			}
		}
	} 
	if (npionp*npionm != 2) return StatusCode::SUCCESS;
    if ((npionp==1)&&(npionm==2)) m_mode = -1;
    else if ((npionp==2)&&(npionm==1)) m_mode = 1;
    else return StatusCode::SUCCESS;
	Ncut2++;

	/***********************************************************************/
	/**********************Vertex fit **************************************/
	/***********************************************************************/
	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;
	VertexParameter vxpar;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);
	bool okloop=false; 
	bool okvs=false;
	VertexFit *vtxfit_s = VertexFit::instance(); // S second vertex 
	SecondVertexFit *vtxfit2 = SecondVertexFit::instance();
	int ipipks= 9999;
	int ipimks= 9999;
	int ipi = 9999;
	double mks_dist=9.0;
	double mks_temp=9.0;
	int n=999;
	WTrackParameter  wvpipTrk, wvpimTrk;
	/***************K0_S vertex fit***********************************/
	/**********step 1 sort two ks daughter candidate******************/
	for(int i=0; i<npionp;i++){
		RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ipionp[i]))
			->mdcKalTrack();
		wvpipTrk = WTrackParameter(mpion, pipTrk->getZHelix(), 
				pipTrk->getZError());
		for(int j = 0; j < npionm; j++) {
			RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+ipionm[j]))->mdcKalTrack();
			wvpimTrk = WTrackParameter(mpion, pimTrk->getZHelix(), pimTrk->getZError());
			VertexParameter vxpar1;
			vxpar1.setVx(vx);
			vxpar1.setEvx(Evx);
			vtxfit_s->init();
			vtxfit_s->AddTrack(0, wvpipTrk);
			vtxfit_s->AddTrack(1, wvpimTrk);
			vtxfit_s->AddVertex(0, vxpar1, 0, 1);
			okvs = vtxfit_s->Fit();
			if(!okvs) continue;
			vtxfit_s->Swim(0);
			HepLorentzVector p1 = wvpipTrk.p();
			HepLorentzVector p2 = wvpimTrk.p();
			HepLorentzVector p2pi = p1 + p2;
			mks_temp = p2pi.m();
			if(fabs(mks_temp-mks0)<mks_dist){
				ipipks = ipionp[i];
				ipimks = ipionm[j];	
                if (m_mode == -1) {
				    n = 1 - j;
				    ipi = ipionm[n];
                }
                if (m_mode == 1) {
				    n = 1 - i;
				    ipi = ipionp[n];
                }
				mks_dist=fabs(mks_temp-mks0);
			}
		}
	}

	if(!(mks_temp<9.0)) return StatusCode::SUCCESS;
	WTrackParameter wvpipksTrk,wvpimksTrk,wvpiTrk,wks;
	RecMdcKalTrack *pipksTrk = (*(evtRecTrkCol->begin()+ipipks))->mdcKalTrack();
	RecMdcKalTrack *pimksTrk = (*(evtRecTrkCol->begin()+ipimks))->mdcKalTrack();
	RecMdcKalTrack *piTrk = (*(evtRecTrkCol->begin()+ipi))->mdcKalTrack();
	wvpipksTrk = WTrackParameter(mpion, pipksTrk->getZHelix(), pipksTrk->getZError());
	wvpimksTrk = WTrackParameter(mpion, pimksTrk->getZHelix(), pimksTrk->getZError());
	wvpiTrk = WTrackParameter(mpion, piTrk->getZHelix(), piTrk->getZError());
	/**************step two vertex fit************************************/  
	vtxfit_s->init();
	vtxfit_s->AddTrack(0, wvpipksTrk);
	vtxfit_s->AddTrack(1, wvpimksTrk);
	vtxfit_s->AddVertex(0, vxpar, 0, 1);
	if(!(vtxfit_s->Fit(0))) return StatusCode::SUCCESS;
	vtxfit_s->Swim(0);
	vtxfit_s->BuildVirtualParticle(0);
	m_vfits_chi = vtxfit_s->chisq(0);
	WTrackParameter wvkshort = vtxfit_s->wVirtualTrack(0);
	VertexParameter vparks  = vtxfit_s->vpar(0);
	/***************Primary vertex fit **********************************/
	HepPoint3D newvx(0., 0., 0.);
	HepSymMatrix newEvx(3, 0);
	VertexParameter primaryVpar;
	IVertexDbSvc* vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(!(vtxsvc->isVertexValid())) return StatusCode::SUCCESS;
	double* db_vx = vtxsvc->PrimaryVertex();
	double* db_vx_err = vtxsvc->SigmaPrimaryVertex();
	newvx.setX(db_vx[0]);
	newvx.setY(db_vx[1]);
	newvx.setZ(db_vx[2]);
	newEvx[0][0] = db_vx_err[0]*db_vx_err[0];
	newEvx[1][1] = db_vx_err[1]*db_vx_err[1];
	newEvx[2][2] = db_vx_err[2]*db_vx_err[2];
	primaryVpar.setVx(newvx);
	primaryVpar.setEvx(newEvx);
	/*************second vertex fit ************************************/
	vtxfit2->init();
	vtxfit2->setPrimaryVertex(primaryVpar);
	vtxfit2->AddTrack(0, wvkshort);
	vtxfit2->setVpar(vtxfit_s->vpar(0));
	if(!vtxfit2->Fit()) return StatusCode::SUCCESS;

	Ncut3++;

	m_vfit2_chi = vtxfit2->chisq();
	wks         = vtxfit2->wpar();
	m_vfit2_mks = (wks.p()).m(); 
	m_vfit2_ct  = vtxfit2->ctau();
	m_vfit2_dl  = vtxfit2->decayLength();
	m_vfit2_dle = vtxfit2->decayLengthError();
	HepLorentzVector v4ecms(0.011*m_Ecms, 0, 0, m_Ecms);

	HepLorentzVector p2=wvpiTrk.p();
	HepLorentzVector p3=wvkshort.p();
	HepLorentzVector p1=v4ecms-p2-p3;
	HepLorentzVector p4=wvpipksTrk.p();
	HepLorentzVector p5=wvpimksTrk.p();

	HepLorentzVector ptempv;

	for (int ii=0; ii<5; ii++) {
		if (ii==0) {  
			ptempv = p1; 
		}
		if (ii==1) {    
			ptempv = p2;
		}
		if (ii==2) {    
			ptempv = p3;
		}
		if (ii==3) {
			ptempv = p4;
		}
		if (ii==4) {
			ptempv = p5;
		}

		m_vfit_p3[ii][0]   = ptempv.px();
		m_vfit_p3[ii][1]   = ptempv.py();
		m_vfit_p3[ii][2]   = ptempv.pz();

		m_vfit_p[ii]   = ptempv.rho();
		m_vfit_e[ii]   = ptempv.e();
	}
	m_vfit_m = p1.m();
	double m_temp=p1.m(); 
	if(m_temp>0.7||m_temp<0.3)return StatusCode::SUCCESS;
	m_vfit_cos = p1.cosTheta();
	m_vfit_phi = p1.phi();
	m_vfit_pt = sqrt(p1.px()*p1.px()+p1.py()*p1.py());
	m_vfit_pkaon = p1.rho();
	m_vfit_mkstar[0] = (p1+p2).m();
	m_vfit_mkstar[1] = (p2+p3).m();
	Vint final_iGood;
	final_iGood.clear();
	m_nkcond=0;
	m_nMatch=0;
	int nMatch_temp=0;
	int nkcond_temp=0;
	Vint ipfinal;
	ipfinal.clear();
	int finaltemp=999;
	Hep3Vector xorigin(0,0,0);
	double* dbv = vtxsvc->PrimaryVertex();
	double* vv  = vtxsvc->SigmaPrimaryVertex();
	xorigin.setX(dbv[0]);
	xorigin.setY(dbv[1]);
	xorigin.setZ(dbv[2]);
	/****************************************************************/
	/*****************Tracking Efficiency****************************/
	/****************************************************************/
	for(int i = 0; i < nGood; i++){
		EvtRecTrackIterator final_itTrk=evtRecTrkCol->begin() + iGood[i];
		if(!(*final_itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *final_mdcTrk = (*final_itTrk)->mdcTrack();
		HepVector a = final_mdcTrk->helix();
		HepSymMatrix Ea = final_mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);     
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]); 
		double  Rvz0=vecipa[3];   
		double  Rvphi0=vecipa[1];
        if (iGood[i]==ipi) {
		    if(fabs(Rvz0)  >= m_vz0cut){return StatusCode::SUCCESS;} 
		    if(fabs(Rvxy0) >= m_vr0cut){return StatusCode::SUCCESS;}
        }
		if (!(iGood[i]==ipipks||iGood[i]==ipimks||iGood[i]==ipi)) {
			/***************Kaon Candidates************************************/
			HepLorentzVector kcond_ptemp;
			kcond_ptemp.setPx(final_mdcTrk->px());
			kcond_ptemp.setPy(final_mdcTrk->py());
			kcond_ptemp.setPz(final_mdcTrk->pz());
			m_angle[nkcond_temp] =p1.angle(kcond_ptemp)*180/(CLHEP::pi);
			m_Rz[nkcond_temp] = fabs(Rvz0);
			m_Rxy[nkcond_temp] = fabs(Rvxy0);
			if((*final_itTrk)->isMdcKalTrackValid()) nG_tag[nkcond_temp]=1;
			else nG_tag[nkcond_temp]=0;
			nkcond_temp++;
            if((fabs(Rvz0)<=m_vz0cut)&&(fabs(Rvxy0)<= m_vr0cut)){
			    HepLorentzVector finalptemp;
			    finalptemp.setPx(final_mdcTrk->px());
			    finalptemp.setPy(final_mdcTrk->py());
			    finalptemp.setPz(final_mdcTrk->pz());
			    m_angle_match[nMatch_temp] =p1.angle(finalptemp)*180/(CLHEP::pi);
			    if((*final_itTrk)->isMdcKalTrackValid()) nM_tag[nMatch_temp]=1;
			    else nM_tag[nMatch_temp]=0;
			    nMatch_temp++;
			    finaltemp =iGood[i];
			    ipfinal.push_back(finaltemp);
            }
		}
	}
	m_nkcond=nkcond_temp;
	m_nMatch=nMatch_temp;
	// if(m_vfit2_dl/m_vfit2_dle<2.3) return StatusCode::SUCCESS;
	// if(fabs(m_vfit2_mks-0.497614)>0.008) return StatusCode::SUCCESS;
	m_var_kstar=fabs(m_vfit_mkstar[0]-0.89594);
	// if(var_kstar>=0.12) return StatusCode::SUCCESS;
    // if(m_vfit_pt>0.1) return StatusCode::SUCCESS;
	Ncut4++;

	/****************************************************************/
	/*****************PID Efficiency*********************************/
	/****************************************************************/
	m_pid_absolute = 0;
	m_pid_relative = 0;
	m_pid_relative_b = 0;
	m_pid_misid_pi = 0;
	m_pid_misid_proton =0;
	m_pid_dedx = 0;
	m_pid_tof =0 ;
	m_final_jpsi = 9999;
	m_final_kaonp = 9999;
	m_final_kaoncos = 9999;
	m_final_kaonphi = 9999;
	for(int i=0;i<m_nMatch;i++){
		EvtRecTrackIterator itTrk_final = evtRecTrkCol->begin() + ipfinal[i];
		RecMdcTrack *kaonTrk = (*(evtRecTrkCol->begin()+ipfinal[i]))->mdcTrack();
		HepLorentzVector finalkaon;
		finalkaon.setPx(kaonTrk->px());
		finalkaon.setPy(kaonTrk->py());
		finalkaon.setPz(kaonTrk->pz());
		double pktemp = finalkaon.mag();
		finalkaon.setE(sqrt(pktemp*pktemp+mkaon*mkaon));
		m_final_jpsi = (finalkaon+p2+p3).m();
		m_final_kaonp = finalkaon.rho();
		m_final_kaoncos = finalkaon.cosTheta();
		m_final_kaonphi = finalkaon.phi();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk_final);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2());
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() );
		pid->calculate();
		m_pid_dedx_K[i] = pid->chiDedx(3);
		m_pid_tof1_K[i] = pid->chiTof1(3);
		m_pid_tof2_K[i] = pid->chiTof2(3);
		m_pid_prob_K[i] = pid->probKaon();
		if(pid->IsPidInfoValid()){
			if(pid->probKaon() > 0.001 && (pid->probKaon() > pid->probPion())
					&& (pid->probKaon() > pid->probProton())){
				m_pid_absolute++;
			}
			if( (pid->probKaon() > pid->probPion())
					&& (pid->probKaon() > pid->probProton())){
				m_pid_relative++;
			}
			if( (pid->probKaon() > pid->probPion())){
				m_pid_relative_b++;
			}
			if( (pid->probPion() > pid->probKaon())
					&& (pid->probPion() > pid->probProton())){
				m_pid_misid_pi++;
			}
			if( (pid->probProton() > pid->probKaon())
					&& (pid->probProton() > pid->probPion())){
				m_pid_misid_proton++;
			}
		}
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk_final);
		pid->usePidSys(pid->useDedx());
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() );
		pid->calculate();
		if(pid->IsPidInfoValid()){
			if( (pid->probKaon() > pid->probPion())
					&& (pid->probKaon() > pid->probProton())){
				m_pid_dedx++;
			}
		}
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk_final);
		pid->usePidSys(pid->useTof1() | pid->useTof2());
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() );
		pid->calculate();
		if(pid->IsPidInfoValid()){
			if( (pid->probKaon() > pid->probPion())
					&& (pid->probKaon() > pid->probProton())){
				m_pid_tof++;
			}
		}
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk_final);
		pid->usePidSys(pid->useTof1() | pid->useTof2() | pid->useTofE());
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton() );
		pid->calculate();
		if(pid->IsPidInfoValid()){
			if( (pid->probKaon() > pid->probPion())
					&& (pid->probKaon() > pid->probProton())){
				m_pid_tof_e++;
			}
		}
	}
	if (m_nrun<0){
		Event::McParticle temp;
		if (!mcParticleCol) {
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		int rootIndex(-1), num_mc(0);
		bool switchDecay(false);
		int c_pdt_num=443;
		Event::McParticleCol::iterator iter_mc;
		for (iter_mc=mcParticleCol->begin(); iter_mc!=mcParticleCol->end(); iter_mc++){
			int motherpdg=((*iter_mc)->mother()).particleProperty();
			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;
			if ((*iter_mc)->particleProperty()==c_pdt_num){
				switchDecay = true;
				rootIndex=(*iter_mc)->trackIndex();
			}
			if(!switchDecay) continue;
			int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex; //change long to int
			int trkidx= (*iter_mc)->trackIndex() - rootIndex;
			int pdgid = (*iter_mc)->particleProperty();
			if(num_mc>=100) break;
			m_pdgid[num_mc]=pdgid;
			m_motheridx[num_mc]=mcidx;
			m_trkidx[num_mc]=rootIndex;
			m_motherid[num_mc]=motherpdg;
			temp = (*iter_mc)->mother();
			if (pdgid == 443) {
				m_drank[num_mc] = 0;
			}else{
				for (int i = 1; i < 100; i++){
					if (temp.particleProperty () == 443) {
						m_drank[num_mc] = i;
						break;
					}
					temp = temp.mother ();
				}
			}
			num_mc+= 1;
		}

		m_idxmc=num_mc;
	}

	m_tuple1->write();
return StatusCode::SUCCESS;
}


/************************************************************************/  
StatusCode KsKpi::finalize() {
	cout<<"total number:         "<<Ncut0<<endl;
	cout<<"nGood>=3              "<<Ncut1<<endl;
	cout<<"Pass Pid:             "<<Ncut2<<endl;
	cout<<"Pass Vertex Fit:      "<<Ncut3<<endl;
	cout<<"Pass final cut        "<<Ncut4<<endl;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
