rm outputs/out_*_zjetsall.root;
hadd -f outputs/out_0_zjetsall.root outputs/out_0_zjets*;
hadd -f outputs/out_1_zjetsall.root outputs/out_1_zjets*;
hadd -f outputs/out_2_zjetsall.root outputs/out_2_zjets*;

hadd -f outputs/out_0_twdr.root outputs/out_0_tw_dr.root outputs/out_0_atw_dr.root;
hadd -f outputs/out_1_twdr.root outputs/out_1_tw_dr.root outputs/out_1_atw_dr.root;
hadd -f outputs/out_2_twdr.root outputs/out_2_tw_dr.root outputs/out_2_atw_dr.root;

hadd -f outputs/out_0_twds.root outputs/out_0_tw_ds.root outputs/out_0_atw_ds.root;
hadd -f outputs/out_1_twds.root outputs/out_1_tw_ds.root outputs/out_1_atw_ds.root;
hadd -f outputs/out_2_twds.root outputs/out_2_tw_ds.root outputs/out_2_atw_ds.root;

hadd -f outputs/out_0_st.root outputs/out_0_t.root outputs/out_0_at.root outputs/out_0_s.root outputs/out_0_as.root;
hadd -f outputs/out_1_st.root outputs/out_1_t.root outputs/out_1_at.root outputs/out_1_s.root outputs/out_1_as.root;
hadd -f outputs/out_2_st.root outputs/out_2_t.root outputs/out_2_at.root outputs/out_2_s.root outputs/out_2_as.root;

hadd -f outputs/out_0_di.root outputs/out_0_ww.root outputs/out_0_wz.root outputs/out_0_zz.root ;
hadd -f outputs/out_1_di.root outputs/out_1_ww.root outputs/out_1_wz.root outputs/out_1_zz.root ;
hadd -f outputs/out_2_di.root outputs/out_2_ww.root outputs/out_2_wz.root outputs/out_2_zz.root ;

hadd -f outputs/out_0_others.root outputs/out_0_di.root outputs/out_0_st.root outputs/out_0_wjets.root  ;
hadd -f outputs/out_1_others.root outputs/out_1_di.root outputs/out_1_st.root outputs/out_1_wjets.root outputs/out_1_qcd_mu.root ;
hadd -f outputs/out_2_others.root outputs/out_2_di.root outputs/out_2_st.root outputs/out_2_wjets.root outputs/out_2_qcd_mu.root ;

hadd -f outputs/out_0_mc.root outputs/out_0_twdr.root outputs/out_0_tt.root outputs/out_0_others.root outputs/out_0_zjetsall.root ;
hadd -f outputs/out_1_mc.root outputs/out_1_twdr.root outputs/out_1_tt.root outputs/out_1_others.root outputs/out_1_zjetsall.root ; 
hadd -f outputs/out_2_mc.root outputs/out_2_twdr.root outputs/out_2_tt.root outputs/out_2_others.root outputs/out_2_zjetsall.root ;

