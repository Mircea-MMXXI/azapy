@echo off
set OUTDIR=..\docs\


set FILEIN=CVaR_th_doc_base.md ^
 + CVaRAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_CVaR_class.md ^
 + Port_CVaR_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%CVaR_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=MAD_th_doc_base.md ^
 + MADAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_MAD_class.md ^
 + Port_MAD_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%MAD_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=SMCR_th_doc_base.md ^
 + SMCRAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_SMCR_class.md ^
 + Port_SMCR_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%SMCR_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=LSSD_th_doc_base.md ^
 + LSSDAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_LSSD_class.md ^
 + Port_LSSD_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%LSSD_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=GINI_th_doc_base.md ^
 + GINIAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_GINI_class.md ^
 + Port_GINI_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%GINI_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=Omega_th_doc_base.md ^
 + OmegaAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_Omega_class.md ^
 + Port_Omega_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%Omega_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=MV_th_doc_base.md ^
 + MVAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_MV_class.md ^
 + Port_MV_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%MV_th_doc.md
copy %FILEIN% %FILEOUT%


set FILEIN=SD_th_doc_base.md ^
 + SDAnalyzer_class.md ^
 + Risk_getWeights.md + Risk_getRisk.md + Risk_getPositions.md + Risk_ViewFrontiers.md ^
 + Risk_set_mktdata.md + Risk_set_rrate.md + Risk_set_rtype.md ^
 + Port_SD_class.md ^
 + Port_SD_set_model.md + Port_port_view.md + Port_port_view_all.md ^
 + Port_port_drawdown.md + Port_port_perf.md + Port_port_annual_returns.md ^
 + Port_port_monthly_returns.md + Port_port_period_returns.md ^
 + Port_get_nshares.md + Port_get_account.md + Port_get_mktdata.md
set FILEOUT=%OUTDIR%SD_th_doc.md
copy %FILEIN% %FILEOUT%
