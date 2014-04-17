{
	
	/***********************
	 DEFINE THE CANVASSES
	 ***********************/
	
	TCanvas* c[100];
	
	for (unsigned int i=0; i<sizeof(c)/sizeof(c[0]);i++)
		c[i] = NULL;
		
		int w=1024;
		int h=768;
		
		c[0] = new TCanvas("c1","MC Control Plots",w,h);
		c[0]->Divide(2,2);
		c[1] = new TCanvas("c2","reweight Control Plots",w,h);
		c[1]->Divide(2,3);
		c[2] = new TCanvas("c3","Results without right-Reweighting",w,h);
		c[2]->Divide(2,2);
		c[3] = new TCanvas("c4","Results WITH right-Reweighting",w,h);
		c[3]->Divide(2,2);
		c[4] = new TCanvas("c5","Cross Section Measurement",w,h);
		c[4]->Divide(3,3);
		
		c[5] = new TCanvas("c6","Btag distribution measurement for all taggers",w,h);
		c[5]->Divide(2,4);
		
		c[6] = new TCanvas("c7","Btag eff measurement",w,h);
		c[6]->Divide(2,1);
		
		c[7] = new TCanvas("c8","Just a test",w,h);
		c[7]->Divide(1,1);
				
		int ccounter[8];
	
	for (unsigned int i=0; i<sizeof(ccounter)/sizeof(ccounter[0]);i++)
		ccounter[i]=1;
		
	/*********************
	 DEFINE THE PLOTS
	 *********************/
		
		string plots[100][100];
	
	int k=0;
	
	// mc control plots
	plots[k][0]="m_{lj} distribution"; // plot title
	plots[k][1]="0"; // canvas number
	plots[k][2]="0"; // logy?
	plots[k][3]="1"; // NORMALIZE TO 1?
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_Var0"; // plot name
	plots[k][5]="b jets"; // plot title
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0";
	plots[k][7]="Non-b jets";
	k++;
	
	plots[k][0]="TrackCountingHighEfficiency discriminant";
	plots[k][1]="0";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="b jets";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_BtagAll";
	plots[k][7]="Non-b jets";
	k++;
	
	plots[k][0]="Purity";
	plots[k][1]="0";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1SoverSB_Var0";
	plots[k][5]="";
	k++;
	
	plots[k][0]="#epsilon_{b} MC";
	plots[k][1]="0";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC";
	plots[k][5]="";
	k++;
	
	// reweighting control plots
	
	plots[k][0]="Tagger distribution (Data)";
	plots[k][1]="1";
	plots[k][2]="1";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY";
	plots[k][5]="Right";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DYReweigh";
	plots[k][7]="RightReweigh";
	k++;
	
	plots[k][0]="Tagger distribution (Data)";
	plots[k][1]="1";
	plots[k][2]="1";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagAll";
	plots[k][5]="All";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_Left1DY";
	plots[k][7]="Left";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_Right1DY";
	plots[k][9]="Right";
	k++;
	
	plots[k][0]="Left-Right p_{T} Reweighting";
	plots[k][1]="1";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX";
	plots[k][5]="";
	k++;
	
	plots[k][0]="(p_{T},#eta) Reweighting";
	plots[k][1]="1";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH2Data_SignalControlVar12_";
	plots[k][5]="";
	k++;
	
	plots[k][0]="(p_{T},#eta) Reweighting: m_{lj}";
	plots[k][1]="1";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar";
	plots[k][5]="Control sample all";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0";
	plots[k][7]="Signal sample non-b";
	k++;
	
	plots[k][0]="(p_{T},#eta) Reweighting: m_{lj}";
	plots[k][1]="1";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh";
	plots[k][5]="Reweighted Control sample all";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0";
	plots[k][7]="Signal sample non-b";
	k++;
	
	// non reweighted results
	
	plots[k][0]="Measured Btag Distribution";
	plots[k][1]="2";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasured";
	plots[k][5]="Measured";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasured";
	plots[k][7]="MC";
	k++;
	
	plots[k][0]="Btagging Efficiency";
	plots[k][1]="2";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasured";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasured";
	plots[k][9]="Method(F_{MC})";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC";
	plots[k][5]="MC";
	k++;
	
	plots[k][0]="Btagging Efficiency #Delta(Method(F_{data}),MC)";
	plots[k][1]="2";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredDiff";
	plots[k][5]="";
	k++;
	
	plots[k][0]="Btagging Efficiency #Delta(Method(F_{MC}),MC)";
	plots[k][1]="2";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredDiff";
	plots[k][5]="";
	k++;
	
	// reweighted results
	
	plots[k][0]="Measured Btag Distribution";
	plots[k][1]="3";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][5]="Measured";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="MC";
	k++;
	
	plots[k][0]="Btagging Efficiency";
	plots[k][1]="3";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRR";
	plots[k][9]="Method(F_{MC})";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC";
	plots[k][5]="MC";
	k++;
	
	plots[k][0]="Btagging Efficiency #Delta(Method(F_{data}),MC)";
	plots[k][1]="3";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRRDiff";
	plots[k][5]="";
	k++;
	
	plots[k][0]="Btagging Efficiency #Delta(Method(F_{MC}),MC)";
	plots[k][1]="3";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRRDiff";
	plots[k][5]="MC";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMeasuredRRDiff";
	plots[k][7]="Data";
	k++;
	
	// XS plots 
	
	plots[k][0]="Comparison M_{lj} VV template - Control sample";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC";
	plots[k][5]="VV MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	plots[k][7]="Control Sample";
	k++;
	
	/*plots[k][0]="Comparison M_{lj} VV template - t #bar{t} other template";
	 plots[k][1]="4";
	 plots[k][2]="0";
	 plots[k][3]="1";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVTTotherMC";
	 plots[k][5]="VV+ttother MC Template";
	 plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	 plots[k][7]="Control Sample";
	 k++;
	 
	 plots[k][0]="Shape ttother == shape control sample?";
	 plots[k][1]="4";
	 plots[k][2]="0"; 
	 plots[k][3]="1";
	 plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	 plots[k][7]="Control Sample";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbarOther";
	 plots[k][5]="t#bar{t} other MC Template";
	 k++;
	 
	 plots[k][0]="Comparison M_{lj} VV template - t #bar{t} template";
	 plots[k][1]="4";
	 plots[k][2]="0";
	 plots[k][3]="1";
	 plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC";
	 plots[k][7]="VV MC Template";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar";
	 plots[k][5]="t#bar{t} MC Template";
	 k++;*/
	
	plots[k][0]="Comparison M_{lj} VV template - Control sample";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC";
	plots[k][5]="VV MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagL";
	plots[k][7]="VV+btagL MC Template";
	//plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	//plots[k][9]="Control Sample";
	k++;
	
	plots[k][0]="Comparison M_{lj} VV template - Control sample";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagM";
	plots[k][5]="VV+btagM MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	plots[k][7]="Control Sample";
	k++;
	
	plots[k][0]="Comparison M_{lj} VV template - Control sample";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagT";
	plots[k][5]="VV+btagT MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	plots[k][7]="Control Sample";
	k++;
	
	plots[k][0]="t#bar{t} MC Template with BTAG";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar_bTagL";
	plots[k][5]="t#bar{t}+bTagL MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar_bTagT";
	plots[k][7]="t#bar{t}+bTagM MC Template";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar_bTagM";
	plots[k][9]="t#bar{t}+bTagT MC Template";
	k++;
	
	plots[k][0]="VV MC Template with BTAG";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagL";
	plots[k][5]="VV+btagL MC Template";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagT";
	plots[k][7]="VV+btagM MC Template";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC_bTagM";
	plots[k][9]="VV+btagT MC Template";
	k++;
	
	plots[k][0]="Var0 (Data) with BTAG";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_bTagL";
	plots[k][5]="Data + bTagL";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_bTagT";
	plots[k][7]="Data + bTagM";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_bTagM";
	plots[k][9]="Data + bTagT";
	k++;
	
	plots[k][0]="Compare data to MC templates";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVMC";
	plots[k][7]="VV MC Template";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar";
	plots[k][5]="t#bar{t} MC Template";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0";
	plots[k][9]="Data";  
	k++;
	
	plots[k][0]="Compare data to MC-Data templates";
	plots[k][1]="4";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVData";
	plots[k][7]="VV Data Template";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbar";
	plots[k][5]="t#bar{t} MC Template";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0";
	plots[k][9]="Data";  
	k++;
	
	/*plots[k][0]="Comparison M_{lj} VV template - t #bar{t} template";
	 plots[k][1]="4";
	 plots[k][2]="0";
	 plots[k][3]="1";
	 plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_VVTTotherMC";
	 plots[k][7]="VV+ttother MC Template";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_TTbarSemiMu";
	 plots[k][5]="t#bar{t}-semi#mu MC Template";
	 k++;*/
	
	/*  plots[k][0]="test1";
	 plots[k][1]="5";
	 plots[k][2]="0";
	 plots[k][3]="1";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar";
	 plots[k][5]="Data";
	 plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_ControlVar";
	 plots[k][7]="Signal";
	 plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_ControlVar";
	 plots[k][9]="Background";
	 k++;
	 
	 plots[k][0]="test2";
	 plots[k][1]="5";
	 plots[k][2]="0";
	 plots[k][3]="1";
	 plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Pt_Control";
	 plots[k][5]="";
	 k++;*/
	
	/*plots[k][0]="Compare data to MC-Data templates";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_Var0";
	plots[k][5]="Pure bjets SS";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_Var0_bTagL";
	plots[k][7]="pure bjets SS + btag";
	k++;
*/
	
	/*plots[k][0]="Check btag cut (loose)";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_Var0_bTagL";
	plots[k][5]="Sng";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0_bTagL";
	plots[k][7]="Bkg";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0_bTagL";
	plots[k][9]="Data";
	k++;
	
	plots[k][0]="Check btag cut (no cut)";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_Var0";
	plots[k][5]="Sng";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0";
	plots[k][7]="Bkg";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_Var0";
	plots[k][9]="Data";
	k++;
*/
	
	/*plots[k][0]="Compare m_{lj} shape t#bar{t} VS non-t#bar{t}";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar_TTBar";
	plots[k][5]="t#bar{t} Control Sample";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVar_VVMC";
	plots[k][7]="non-t#bar{t} Control Sample";
	k++;
	
	
	plots[k][0]="Compare m_{lj} shape t#bar{t} VS non-t#bar{t}";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh_TTBar";
	plots[k][5]="t#bar{t} Control Sample Reweighted";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh_VVMC";
	plots[k][7]="non-t#bar{t} Control Sample Reweighted";
	k++;*/
	
	plots[k][0]="TCHE";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	/*plots[k][0]="TCHP";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="SSVHE";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="SSVHP";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_3/nDisc_3_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_3/nDisc_3_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_3/nDisc_3_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="CSV";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_4/nDisc_4_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="CSVMVA";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_5/nDisc_5_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_5/nDisc_5_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_5/nDisc_5_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="JP";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_6/nDisc_6_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_6/nDisc_6_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_6/nDisc_6_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	
	plots[k][0]="JBP";
	plots[k][1]="5";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMCMeasuredRR";
	plots[k][7]="Method";
	plots[k][8]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagMeasuredRR";
	plots[k][9]="Method (Fdata)";
	k++;
	*/
	
	plots[k][0]="Btagging Efficiency TCHE";
	plots[k][1]="6";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC";
	plots[k][5]="True";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasured";
	plots[k][7]="Measured";
	
	k++;
	
	plots[k][0]="Bias";
	plots[k][1]="6";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredDiff";
	plots[k][5]="Bias";
	k++;

	/*plots[k][0]="Btagging Efficiency SSVHE";
	plots[k][1]="6";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffMC";
	plots[k][5]="True";
	plots[k][6]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRR";
	plots[k][7]="Measured";
	
	k++;
	
	plots[k][0]="Bias";
	plots[k][1]="6";
	plots[k][2]="0";
	plots[k][3]="0";
	plots[k][4]="Variable_1_0/Discriminator_2/nDisc_2_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRRDiff";
	plots[k][5]="Bias";
	k++;*/
	
	plots[k][0]="TCHE";
	plots[k][1]="7";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Sng_BtagEffAll";
	plots[k][5]="MC Truth";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_BtagEffMCMeasuredRR ";
	plots[k][7]="Method";
	k++;
	
	/*k--;plots[k][0]="(p_{T},#eta) Reweighting: m_{lj}";
	plots[k][1]="7";
	plots[k][2]="0";
	plots[k][3]="1";
	plots[k][6]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Data_ControlVarReweigh";
	plots[k][7]="Control sample";
	plots[k][4]="Variable_1_0/Discriminator_7/nDisc_7_ptbinlow_0_etabinlow_-9990_TH1Bkg_Var0";
	plots[k][5]="Signal sample non-b";
	k++;*/
	
	/*********************
	 DO THE PLOTTING 
	 *********************/
	
	TFile* f = new TFile("TTrees/DATA_6/NTupleAnalyzed.root","READ");
	
	for (unsigned int i=0; i<100; i++) {
		
		if (plots[i][0] == "") continue;
		
		int canv = atoi(plots[i][1].c_str());
		
		cout << "+> Plotting " << plots[i][0] << endl;
		
		TLegend* leg = new TLegend(0.5,0.65,0.8.,0.8);
		
		float y1=0.79;
		float y2=0.99;
		
		Int_t col = 1;
		
		bool firstPlot = true;
		for (unsigned int j=4; j<100; j+=2) {
			
			//cout << 
			
			if (plots[i][j] == "") continue;
			
			cout << "  +++> Processing plot " << plots[i][j] << " title: " << plots[i][j+1] << endl;
			
			c[canv]->cd(ccounter[canv]);
			
			gPad->SetLogy(atoi(plots[i][2].c_str()));
			
			gPad->SetGrid();
			
			TH1F* tmp = (TH1F*) f->Get((plots[i][j]).c_str());
			
			TH1F* plot = (TH1F*)tmp->Clone();
			
			//cout << f->Get((plots[i][j]).c_str())->GetName() << endl;
			
			bool emptyPlot = false;
			
			for (int b=0;b<plot->GetNbinsX()+1; b++) {
				
				//cout << "Bin " << b << " Content " << plot->GetBinContent(b) << endl;
				
				stringstream s; s << plot->GetBinContent(b);
				
				if (s.str().find("nan") == 0)
					emptyPlot=true;
				
				
			}
			cout << "Empty plot? : " << emptyPlot << endl;
			if (!emptyPlot && plot && plot->GetEntries() > 0) {
				
				//cout << plot->GetEntries() << endl;
				
				plot->SetTitle(plots[i][0].c_str());
				plot->SetLineColor(col);    
				plot->SetMarkerColor(col);
				
				if (plots[i][3] == "1")
					plot->Scale(1./plot->Integral());
				
				if (firstPlot) {
					plot->Draw("COLZ");
					firstPlot=false;
				}
				else
					plot->Draw("sames");
				
				if ((plots[i][j+1]).c_str() != "") {
					leg->AddEntry(plot,(plots[i][j+1]).c_str(),"lp");
					leg->Draw();
				}
				
				//if (canv==6) plot->GetXaxis()->SetTitle("TCHE btag discriminator cut");
				//if (canv==6) plot->GetYaxis()->SetTitle("#epsilon_{b}");
				//if (canv==6) plot->SetTitle("");
				
				c[canv]->Update();
				
				TPaveStats* sb = (TPaveStats*)plot->GetListOfFunctions()->FindObject("stats");
				sb->SetTextColor(col);
				sb->SetLineColor(col);
				sb->SetX1NDC(.83);
				sb->SetX2NDC(.98);
				sb->SetY1NDC(y1);
				sb->SetY2NDC(y2);
				
				y1-=0.2;
				y2-=0.2;
				
				c[canv]->Draw();
				
				col++;
			}
			
		}
		
		ccounter[canv]++;
		
		}
		
		for (unsigned int i=0; i<sizeof(c)/sizeof(c[0]);i++) {
			
			if (c[i] == NULL) continue;
			
			cout << c[i]->GetTitle() << endl;
			
			stringstream s; s << i;
			string title = "Results/Canvas_"+s.str()+".pdf";
			c[i]->SaveAs(title.c_str());
		}
		}
