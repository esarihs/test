//#include "/group/had/g-2/share/g2esoft/g-2_master/g2esoft/app/include/objects.h"
//#pragma cling load("/group/had/g-2/share/g2esoft/g-2_master/g2esoft/lib/libg2esoftCommon.so")

#include "/home/g-2/g2esoft-v01-02-00/app/include/objects.h"
#pragma cling load("/home/g-2/g2esoft-v01-02-00/lib/libg2esoftCommon.so")

#include <stdio.h>
#include <fstream>

//ZZZ
int kickInput(const char *kickerFileName, const double facB, vector<double> &Ikicker, vector<double> &posXc, vector<double> &posZc);
int kickSpatial2(double x[3],double &BRkick0,double &BYkick0, vector<double> Ikicker, vector<double> posXc, vector<double> posZc);
void kickMag2(double time, double BRkick0, double BYkick0, double &BRkick, double &BYkick, double &km2);

void cp3DrawMuonOrbit(const Char_t *inputFileName, /*const Int_t event,*/ const char *kickerFileName, const Bool_t drawCanvas=/*kFALSE*/kTRUE)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat("");

  TGraph *gr_z_zp = new TGraph;
  TGraph *gr_t_z = new TGraph;
  TGraph *gr_x_y = new TGraph;
  TGraph *gr_x_z = new TGraph;
  TGraph2D *a2D = new TGraph2D;
  //ZZZ
  TGraph *gr_t_br = new TGraph;
  TGraph *gr_t_by = new TGraph;
  TGraph *gr_t_br_kick = new TGraph;
  TGraph *gr_t_by_kick = new TGraph;
  TGraph *gr_t_mag2 = new TGraph;

  TH2F *hxzgood = new TH2F("goodxz","goodxz",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hxzbad = new TH2F("badxz","badxz",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hyzgood = new TH2F("goodyz","goodyz",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hyzbad = new TH2F("badyz","badyz",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hyxgood = new TH2F("goodyx","goodyx",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hyxbad = new TH2F("badyx","badyx",30,-1.5,1.5,30,-1.5,1.5);
  TH2F *hrzgood = new TH2F("goodrz","goodrz",30,0,1.5,30,-1.5,1.5);
  TH2F *hrzbad = new TH2F("badrz","badrz",30,0,1.5,30,-1.5,1.5);
  TH2F *hysgood = new TH2F("goodys","goodys",30,-1.5,1.5,30,0,1.5);//s=sqrt(x^2+z^2)
  TH2F *hysbad = new TH2F("badys","badys",30,-1.5,1.5,30,0,1.5);  
  TH2F *hx2y2good = new TH2F("goodx2y2","goodx2y2",30,0,1.5,30,0,1.5);
  TH2F *hx2y2bad = new TH2F("badx2y2","badx2y2",30,0,1.5,30,0,1.5);
  TH2F *hpy_ygood = new TH2F("goodpy_y","goodpy_y",100,-1000,1000,20,-1000,1000);
  TH2F *hpy_ybad = new TH2F("badpy_y","badpy_y",100,-1000,1000,20,-1000,1000);

  TGraph *txpxgood = new TGraph;
  TGraph *txpxbad = new TGraph;

  //cout<<"hxzgood_mae"<<hxzgood<<endl;

  TFile *inputFile = new TFile(inputFileName, "READ");
  TTree *tree = (TTree*)inputFile->Get("trk");
  vector<const g2esoft::MCParticle *> *p_mcps = 0;
  tree->SetBranchAddress("MCParticles", &p_mcps);

  Long64_t nentries = tree->GetEntries();
  Int_t ipoint = 0; 
  TVector3 pos;
  TLorentzVector mom;
  Double_t time, tmax;
  //ZZZ
  const double facB = 3.45;

  double x[3];
  double BRkick0,BYkick0,BRkick,BYkick;
  double yg[128*100],brg_kick[128*100];
  double byg_kick[128*100];
  double RM,YM;
  double km2;
  vector<double> Ikicker;
  vector<double> posXc;
  vector<double> posZc;

  //ZZZ
  kickInput(kickerFileName, facB, Ikicker, posXc, posZc);

  Long64_t ndiv = 10;
  if( nentries<10 ) ndiv = 1;

  TH1F *h_zamp = new TH1F("zamp", "zamp;Vertical Amplitude [mm];Number of muons", 100, 0, 300);
  vector<TGraph *> grVec;
  Int_t ntotal = 0;
  Int_t nstored = 0;
  Int_t ngood = 0;
  Int_t nbad = 0;

  const Double_t tmax_threshold = 790.; // ns
  const Double_t tmin_fit = 200.; // ns
  const Double_t tmax_fit = 800.; // ns
  const Double_t betatronPeriod = 600.; // ns
  const Double_t z_threshold = 100.; // mm

  //--------fout_2_in-------
  ofstream fout;//ファイル書き出し用ストリームオブジェクトfoutの定義
  fout.open("good.txt");
  ofstream fout2;//ファイル書き出し用ストリームオブジェクトfoutの定義
  fout2.open("bad.txt");
  //--------fout_2_out-------

  TF1 *f1 = new TF1("f1","[0]*sin([1]*x+[2])",tmin_fit,tmax_fit);
  Int_t event = 20;
  for(Long64_t ievent=0; ievent<nentries; ievent++){
    //if(ievent!=event) continue;
    cout << "Draw event : " << ievent << endl;
    tree->GetEvent(ievent);
    if(ievent%(nentries/ndiv)==0) cout << "event " << ievent << endl;

    for(UInt_t imcp=0; imcp<p_mcps->size(); imcp++){
      const g2esoft::MCParticle *mcp = p_mcps->at(imcp);
      if( mcp->_pdg==-13 ){
	TGraph *grTZ = new TGraph;    
	ipoint = 0;
	tmax = 0.;

	grTZ->Set(0);
	ntotal++;      

	for(UInt_t istep=0; istep<mcp->_steps.size(); istep++){
	  const g2esoft::MCStep *mcstep = mcp->_steps.at(istep);
	  //cout<<"event"<<istep<<endl;
	  pos = mcstep->_pos;
	  mom = mcstep->_p;
	  time = mcstep->_time;
	  //cout<<"time "<<time<<endl;
	  gr_z_zp->SetPoint(ipoint, pos.Z()*1e-3, mom.Z()/mom.P());
	  gr_t_z->SetPoint(ipoint, time, pos.Z()*1e-3);
	  gr_x_y->SetPoint(ipoint, pos.X()*1e-3, pos.Y()*1e-3);
	  gr_x_z->SetPoint(ipoint, pos.X()*1e-3, pos.Z()*1e-3);
	  a2D->SetPoint(ipoint,pos.X()*1e-3, pos.Y()*1e-3, pos.Z()*1e-3);

	  //ZZZ
	  RM=sqrt(pow(pos.X()*1e-3,2)+pow(pos.Y()*1e-3,2));
	  //	  cout << "RM" << RM << endl;
	  YM=pos.Z()*1e-3;
	  //	  cout << "YM" << YM << endl;
	  x[1]=YM;
	  x[0]=RM;
	  x[2]=0;
	  yg[istep]=YM;
	  /////////////Get kick-B-spatial/////////
	  int kick_spatial=0;
	  if(sqrt(pow(x[1],2))<1.0){
	    kick_spatial=kickSpatial2(x,BRkick0,BYkick0,Ikicker,posXc,posZc);//Kicker-2D-cir
	  }
	  //	  cout << "BRkick0" << BRkick0 << endl;
	  brg_kick[istep]=BRkick0;
	  byg_kick[istep]=BYkick0;
	  gr_t_br->SetPoint(ipoint, time, brg_kick[istep]); 
	  gr_t_by->SetPoint(ipoint, time, brg_kick[istep]); 

	  kickMag2(time, BRkick0, BYkick0, BRkick, BYkick, km2);
	  brg_kick[istep]=BRkick;
	  byg_kick[istep]=BYkick;
	  //	  cout <<"BRkick; " <<BRkick <<endl;
	  gr_t_br_kick->SetPoint(ipoint, time, brg_kick[istep]); 
	  gr_t_by_kick->SetPoint(ipoint, time, byg_kick[istep]); 
	  gr_t_mag2->SetPoint(ipoint, time, km2); 

	  grTZ->SetPoint(ipoint, time, pos.Z());

	  ipoint++;
	  if(time>tmax) tmax = time;

	}
      
	grVec.push_back(grTZ);
	if( tmax>tmax_threshold ){
	  nstored++;
	  f1->SetParameters(50., TMath::TwoPi()/betatronPeriod, 0.);
	  grTZ->Fit("f1","RQ0");
	  h_zamp->Fill(fabs(f1->GetParameter(0)));
	  if( fabs(f1->GetParameter(0))<z_threshold ){

	    //--------fout_3_in------

	    cout << "ievent = " << ievent << endl;
	    fout << ievent << endl;

	    //cout << "pos.Z()*1e-3 = " << pos.Z()*1e-3 << endl;
	    hxzgood->Fill(pos.X()*1e-3,pos.Z()*1e-3);
	    hyzgood->Fill(pos.Y()*1e-3,pos.Z()*1e-3);
	    hyxgood->Fill(pos.Y()*1e-3,pos.X()*1e-3);
	    hrzgood->Fill(sqrt(pow(pos.X()*1e-3,2)+pow(pos.Y()*1e-3,2)),pos.Z()*1e-3);
	    hysgood->Fill(pos.Y()*1e-3,sqrt(pow(pos.X()*1e-3,2)+pow(pos.Z()*1e-3,2)));
	    hx2y2good->Fill(pow(pos.X()*1e-3,2),pow(pos.Y()*1e-3,2));
	    hpy_ygood->Fill(mom.Y(),pos.Y()*1e-3);
	    txpxgood->SetPoint(ipoint,time,pos.X()*1e-3);
	    ngood++;
	  }else{
	    cout << "ievent = " << ievent << "\t"
		 << "impc = " << imcp << endl;
	    fout2 << ievent << endl;
	    hxzbad->Fill(pos.X()*1e-3,pos.Z()*1e-3);
	    hyzbad->Fill(pos.Y()*1e-3,pos.Z()*1e-3);
	    hyxbad->Fill(pos.Y()*1e-3,pos.X()*1e-3);
	    hrzbad->Fill(sqrt(pow(pos.X()*1e-3,2)+pow(pos.Y()*1e-3,2)),pos.Z()*1e-3);
	    hysbad->Fill(pos.Y()*1e-3,sqrt(pow(pos.X()*1e-3,2)+pow(pos.Z()*1e-3,2)));
	    hx2y2bad->Fill(pow(pos.X()*1e-3,2),pow(pos.Y()*1e-3,2));
	    hpy_ybad->Fill(mom.Y(),pos.Y()*1e-3);
	    //txpxbad->SetPoint(time,pos.X()*1e-3,mom.X());
	    nbad++;
	  }
	}else{
	  cout << "ievent = " << ievent << "\t"
	       << "impc = " << imcp << endl;
	  fout2 << ievent << endl;
	  hxzbad->Fill(pos.X()*1e-3,pos.Z()*1e-3);
	  hyzbad->Fill(pos.Y()*1e-3,pos.Z()*1e-3);
	  hyxbad->Fill(pos.Y()*1e-3,pos.X()*1e-3);
	  hrzbad->Fill(sqrt(pow(pos.X()*1e-3,2)+pow(pos.Y()*1e-3,2)),pos.Z()*1e-3);
	  hysbad->Fill(pos.Y()*1e-3,sqrt(pow(pos.X()*1e-3,2)+pow(pos.Z()*1e-3,2)));
	  hx2y2bad->Fill(pow(pos.X()*1e-3,2),pow(pos.Y()*1e-3,2));
	  hpy_ybad->Fill(mom.Y(),pos.Y()*1e-3);
	  //txpxbad->SetPoint(time,pos.X()*1e-3,mom.X());
	  nbad++;
        }
      	//--------fout_3_out-------
      }
    } // mcp
  }

  //--------fout_4_in-------
  fout.close();
  fout2.close();
  //--------fout_4_out-------

  cout << "good/total = " << ngood << "/" << ntotal << " = " << (Double_t)ngood/ntotal << "+-" << sqrt(ngood * (1-ngood/(Double_t)ntotal))/ntotal << endl;

  cout <<"nbad = " << nbad << endl;

  //txpxgood->Print("all");

  /*
  int N = tmax;
  for(int i=0; i<N; i++ ){
  }
  */

  /*
  for(UInt_t istep=0; istep<mcp->_steps.size(); istep++){
    //const g2esoft::MCStep *mcstep = mcp->_steps.at(istep);
    cout<<"time"<<txpxgood->GetX()[i]<<endl;
  }
  */

  if( drawCanvas ){
    /*
    Char_t c1Name[1024] = "OneMuonOrbit";
    sprintf(c1Name, "%s_%d", c1Name, event);
    TCanvas *c1 = new TCanvas(c1Name, c1Name, 1200, 1200);
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gr_t_z->GetYaxis()->SetTitle("z [m]");
    gr_t_z->GetXaxis()->SetTitle("t [ns]");
    gr_t_z->Draw("APC");

    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gr_z_zp->GetYaxis()->SetTitle("z'");
    gr_z_zp->GetXaxis()->SetTitle("z [m]");
    gr_z_zp->Draw("APC");

    c1->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gr_x_y->GetYaxis()->SetTitle("y [m]");
    gr_x_y->GetXaxis()->SetTitle("x [m]");
    gr_x_y->Draw("APC");

    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gr_x_z->GetYaxis()->SetTitle("z [m]");
    gr_x_z->GetXaxis()->SetTitle("x [m]");
    gr_x_z->Draw("APC");

    c1->SaveAs(Form("%s.png",c1->GetName()));

    Char_t c2Name[1024] = "OneMuonOrbit";
    sprintf(c2Name, "%s_%d_3D", c2Name, event);
    TCanvas *c2 = new TCanvas(c2Name, c2Name, 600, 600);
    c2->Draw();
    a2D->SetTitle(";x [cm];y [cm]; z [cm]");
    a2D->SetLineColor(2);
    gPad->SetLeftMargin(0.15);    
    gPad->SetRightMargin(0.05);    
    a2D->GetXaxis()->SetTitleOffset(1.8);    
    a2D->GetYaxis()->SetTitleOffset(2.2);    
    a2D->GetZaxis()->SetTitleOffset(1.8);

    a2D->Draw("line");
    c2->SaveAs(Form("%s.png",c2->GetName()));

    
    //ZZZ
    Char_t c3Name[1024] = "OneMuonOrbitMagnetField";
    sprintf(c3Name, "%s_%d_D", c3Name, event);
    TCanvas *c3 = new TCanvas(c3Name, c3Name, 1200, 1200);
    c3->Draw();
    c3->Divide(1,3);
    
    c3->cd(1);
    gr_t_br->GetYaxis()->SetTitle("Br [gauss]");
    gr_t_br->GetXaxis()->SetTitle("t [ns]");
    gr_t_br->SetLineColor(2);
    gr_t_br_kick->SetLineColor(4);  
    gr_t_br->SetMarkerColor(2);
    gr_t_br_kick->SetMarkerColor(4);
    gr_t_br->Draw();
    gr_t_br_kick->Draw("same pl");

    c3->cd(2);
    gr_t_by->GetYaxis()->SetTitle("Bz [gauss]");
    gr_t_by->GetXaxis()->SetTitle("t [ns]");
    gr_t_by->SetLineColor(2);
    gr_t_by_kick->SetLineColor(4);
    gr_t_by->SetMarkerColor(2);
    gr_t_by_kick->SetMarkerColor(4);
    gr_t_by->Draw();
    gr_t_by_kick->Draw("same pl");

    c3->cd(3);
    gr_t_mag2->GetYaxis()->SetTitle("mag2 [gauss]");
    gr_t_mag2->GetXaxis()->SetTitle("t [ns]");
    gr_t_mag2->Draw();

    c3->SaveAs(Form("%s.png",c3->GetName()));
    */

    Char_t c4Name[1024] = "goodbad1";
    sprintf(c4Name, "%s_%d", c4Name, event);
    TCanvas *c4 = new TCanvas(c4Name, c4Name, 1200, 1200);

    //cout<<"hxzgood_ato"<<hxzgood<<endl;
    //hxzgood->Print();

    c4->Divide(2,2);
    c4->Draw();
    c4->cd(1);
    hxzgood->SetXTitle("x[cm]");
    hxzgood->SetYTitle("z[cm]");
    hxzgood->Draw("colz");

    c4->cd(2);
    hxzbad->SetXTitle("x[cm]");
    hxzbad->SetYTitle("z[cm]");
    hxzbad->Draw("colz");

    c4->cd(3);
    hyzgood->SetXTitle("y[cm]");
    hyzgood->SetYTitle("z[cm]");
    hyzgood->Draw("colz");

    c4->cd(4);
    hyzbad->SetXTitle("y[cm]");
    hyzbad->SetYTitle("z[cm]");
    hyzbad->Draw("colz");

    c4->SaveAs(Form("%s.png",c4->GetName()));

    Char_t c5Name[1024] = "goodbad2";
    sprintf(c5Name, "%s_%d", c5Name, event);
    TCanvas *c5 = new TCanvas(c5Name, c5Name, 1200, 1200);
    c5->Divide(2,2);
    c5->Draw();

    c5->cd(1);
    hyxgood->SetXTitle("y[cm]");
    hyxgood->SetYTitle("x[cm]");
    hyxgood->Draw("colz");

    c5->cd(2);
    hyxbad->SetXTitle("y[cm]");
    hyxbad->SetYTitle("x[cm]");
    hyxbad->Draw("colz");

    c5->cd(3);
    hrzgood->SetXTitle("r[cm]");
    hrzgood->SetYTitle("z[cm]");
    hrzgood->Draw("colz");

    c5->cd(4);
    hrzbad->SetXTitle("r[cm]");
    hrzbad->SetYTitle("z[cm]");
    hrzbad->Draw("colz");

    /*
    c5->cd(3);
    hpy_ygood->SetXTitle("py[MeV/c]");
    hpy_ygood->SetYTitle("y[cm]");
    hpy_ygood->Draw("colz");

    c5->cd(4);
    hpy_ybad->SetXTitle("py[MeV/c]");
    hpy_ybad->SetYTitle("y[cm]");
    hpy_ybad->Draw("colz");
    */

    Char_t c6Name[1024] = "goodbad3";
    sprintf(c6Name, "%s_%d", c6Name, event);
    TCanvas *c6 = new TCanvas(c6Name, c6Name, 1200, 1200);
    c6->Divide(2,2);
    c6->Draw();

    c6->cd(1);
    hysgood->SetXTitle("y[cm]");
    hysgood->SetYTitle("s[cm]");
    hysgood->Draw("colz");

    c6->cd(2);
    hysbad->SetXTitle("y[cm]");
    hysbad->SetYTitle("s[cm]");
    hysbad->Draw("colz");

    c6->cd(3);
    hx2y2good->SetXTitle("x^2[cm^2]");
    hx2y2good->SetYTitle("y^2[cm^2]");
    hx2y2good->Draw("colz");

    c6->cd(4);
    hx2y2bad->SetXTitle("x^2[cm^2]");
    hx2y2bad->SetYTitle("y^2[cm^2]");
    hx2y2bad->Draw("colz");

  }
  inputFile->Close();
}


//ZZZ
//---------------------------------------------------------------------  
int cep12d_n10(double rk,double i,double &ak,double &ae,double &ill){                                

  int ans=0;
  double a0=1.38629436111989e0; 
  double a1=0.0965735902811690e0;
  double a2=0.0308851465246711e0;
  double a3=0.0149380448916805e0;
  double a4=8.79078273952743e-3;
  double a5=6.18901033637687e-3;
  double a6=6.87489687449949e-3;
  double a7=9.85821379021226e-3;
  double a8=7.97404013220415e-3;
  double a9=2.28025724005875e-3;
  double a10=1.37982864606273e-4;
  
  double b0=0.5000e0;
  double b1=0.124999999999870e0;
  double b2=0.0703124996963957e0;
  double b3=0.0488280347570998e0;
  double b4=0.0373774314173823e0;
  double b5=0.0301204715227604e0;
  double b6=0.0239089602715924e0;
  double b7=0.0154850516649762e0;
  double b8=5.94058303753167e-3;
  double b9=9.14184723865917e-4;
  double b10=2.94078955048598e-5;
  
  double ea0=1.0000e0;
  double ea1=0.443147180560990e0;
  double ea2=0.0568051945617860e0;
  double ea3=0.0218317996015557e0;
  double ea4=0.0115688436810574e0;
  double ea5=7.58395289413514e-3;
  double ea6=7.77395492516787e-3;
  double ea7=0.0107350949056076e0;
  double ea8=8.68786816565889e-3;  
  double ea9=2.50888492163602e-3;
  double ea10=1.53552577301013e-4;
  
  double eb0=0.0e0; 
  double eb1=0.249999999999888e0;
  double eb2=0.0937499997197644e0;
  double eb3=0.0585936634471101e0;
  double eb4=0.0427180926518931e0; 
  double eb5=0.0334833904888224e0;
  double eb6=0.0261769742454493e0;
  double eb7=0.0168862163993311e0;
  double eb8=6.50609489976927e-3;  
  double eb9=1.00962792679356e-3;
  double eb10=3.27954898576485e-5;
  
  i=1;                                                              
  ill=0;                                                           
  double xm1,xm2,xm3,xm4,xm5,xm6,xm7,xm8,xm9,xm10,dalxm,alxm;
  double bzz;
  if(rk<0.0e0   ||   rk>1.0e0){               
    ill=1;                                                            
    ak=0.0;                                                            
    ae=0.0;                                                            
  }else{
    xm1=1.0e0-rk;                                                      
    xm2=xm1*xm1;                                                      
    xm3=xm2*xm1;                                                      
    xm4=xm3*xm1;                                                      
    xm5=xm4*xm1;                                                      
    xm6=xm5*xm1;                                                      
    xm7=xm6*xm1;                                                      
    xm8=xm7*xm1;                                                      
    xm9=xm8*xm1;                                                      
    xm10=xm9*xm1;                                                      
    dalxm=1.0e0/(double)xm1;                                                  
    alxm=log(dalxm);                                                
    bzz=b0 +b1*xm1+b2*xm2+b3*xm3+b4*xm4+b5*xm5+b6*xm6+b7*xm7+b8*xm8+b9*xm9+b10*xm10;             
    ak= a0 +a1*xm1+a2*xm2+a3*xm3+a4*xm4+a5*xm5+a6*xm6+a7*xm7+a8*xm8+a9*xm9+a10*xm10 +bzz*alxm;    
    
    bzz=eb0+eb1*xm1+eb2*xm2+eb3*xm3+eb4*xm4+eb5*xm5+eb6*xm6+eb7*xm7+eb8*xm8+eb9*xm9+eb10*xm10;          
    ae=ea0 +ea1*xm1+ea2*xm2+ea3*xm3+ea4*xm4+ea5*xm5+ea6*xm6+ea7*xm7+ea8*xm8+ea9*xm9+ea10*xm10+bzz*alxm;              
    ans=1;
  }
  
  return ans;
}                                                            


//---------------------------------------------------------            
int bfield(double CR,double CZ,double CI,double RI,double ZJ,double &BR,double &BZ,double &APHI){                      
  int cep=0;
  double XMU=4.E-7;                                                        
  double S =RI*RI+CR*CR+(ZJ-CZ)*(ZJ-CZ);                                  
  double P =RI*CR+RI*CR;                                                  
  double RK=(P+P)/(double)(S+P);  
  double PSI;
  double ELPK,ELPE,ILL; 
  cep=cep12d_n10(RK,1,ELPK,ELPE,ILL);                                  
  
  if(ILL==0){                                          
    RK=sqrt(RK);                                                    
    BZ =XMU*CI/(double)(2.E0* sqrt(S+P))*(ELPK-(S-2.E0*CR*CR)/(double)(S-P)*ELPE);    
    BR =XMU*CI/(double)(2.E0* sqrt(S+P))*(ZJ-CZ)/(double)RI*(-ELPK+ S/(double)(S-P)*ELPE);    
    PSI=CI/(double)RK* sqrt(RI*CR)*((1.E0-RK*RK/(double)2.E0)*ELPK-ELPE);            
    APHI=XMU*PSI/(double)RI;                                                  
  }
  return cep;
} 

//---------------------------------------------------------------
int bflfit2(double RM,double ZM,double &BR,double &BZ,double &APHI,vector<double> Ikicker, vector<double> posXc, vector<double> posZc){
  int ans=0;
  BR=0.0E0;
  BZ=0.0E0;
  APHI=0.0E0;
  double DDD,DBR,DBZ,DAPHI;
  for(int i=0;i<Ikicker.size();++i){
    DDD=sqrt(pow(posXc[i]-RM,2)+pow(posZc[i]-ZM,2));
    if(DDD!=0.0E0){
      int ans1= bfield(posXc[i],posZc[i],Ikicker[i],RM,ZM,DBR,DBZ,DAPHI);  
    }
    BR=BR+DBR;
    BZ=BZ+DBZ;
    APHI=APHI+DAPHI;
  }
  ans=1;
  return ans;
}


//////////////KICK-SPATIAL//2021FEB24
int kickInput(const char *kickerFileName, const double facB, vector<double> &Ikicker, vector<double> &posXc, vector<double> &posZc){
  FILE *infile;
  char sdmy[256];
  
  infile=fopen(kickerFileName,"r");
  fscanf(infile,"%s",sdmy);
  
  double I,Xc,Zc;
  while(fscanf(infile,"%lf,%lf,%lf",&Xc,&Zc,&I)==3){
    Ikicker.push_back(I*facB*2);
    posXc.push_back(Xc);
    posZc.push_back(Zc);
  }
  fclose(infile);
  return 0;
}

//////////////KICK-SPATIAL
int kickSpatial2(double x[3],double &BRkick0,double &BYkick0,vector<double> Ikicker, vector<double> posXc, vector<double> posZc){
  double YM,RM;
  double BR,BY,APHI;
  YM=x[1];
  RM=sqrt(pow(x[0],2)+pow(x[2],2));
  bflfit2(RM,YM,BR,BY,APHI,Ikicker,posXc,posZc);
  
  BRkick0=BR*1E4;//Gauss
  BYkick0=BY*1E4;
  
  //  cout << "BR" << BR << endl;
  return 1;
}

void kickMag2(double time, double BRkick0, double BYkick0, double &BRkick, double &BYkick, double &km2)
{

    // kicker parameter
  const double _kickerTK      = 261.7; // [ns], kicker duration time, which corresponds to the period of sinusoidal function
  const double _kickerStartT =   0.0; // [ns], kicker start timing
  const double _kickerTau    = 1.e29; // [ns], time constant for exponential (By default, it is set to too high values so that exponential term is negligible.)

  double BRkickTemp=BRkick0;
  double BYkickTemp=BYkick0;
  double theta=(time / _kickerTK)*TMath::TwoPi();
  km2 = sin(theta)*exp(-time/_kickerTau);  
  BRkick = BRkickTemp*sin(theta)*exp(-time/_kickerTau);
  BYkick = BYkickTemp*sin(theta)*exp(-time/_kickerTau);

}


