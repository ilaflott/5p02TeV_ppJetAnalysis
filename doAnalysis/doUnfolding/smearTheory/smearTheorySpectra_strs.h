const std::string SCRATCH_BASE ="/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";

const std::string _PY8_FILESTR="02.18.19_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_02-18-19_JERS_0.0eta3.0_unf/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
const std::string _PY8_DESCTXT="GEN PY8 CUETP8M1";

//PDF file strings + descriptive text, these are from raghav, in |eta| bins
const std::string _CT10FILESTR="fnl5020_LO2_R04Jets_modify_CT10nlo_HS";
const std::string _CT14FILESTR="fnl5020_LO2_R04Jets_modify_CT14nlo_HS";
const std::string _HERAFILESTR="fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS";
const std::string _MMHTFILESTR="fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS";
const std::string NNPDFFILESTR="fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC";


const std::string _CT10DESCTXT="CT10";
const std::string _CT14DESCTXT="CT14";
const std::string _HERADESCTXT="HERAPDF15";
const std::string _MMHTDESCTXT="MMHT14";
const std::string NNPDFDESCTXT="NNPDF30";

//PDF file strings + descriptive text, these are from Joao, in |y| bins
//const std::string _CT14NNLOFILESTR ="fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt_v4";
const std::string _CT14NNLOFILESTR2="fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt1_v4";
//const std::string _CT14NNLOFILESTR3="fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQHTp_v4";
//const std::string _CT14NNLODESCTXT ="CT14NNLO";

const std::string _CT14NNLOFILESTR ="fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt_v5";
//const std::string _CT14NNLOFILESTR2="fnl5020_R04Jets_ybins_CT14nlo_murEQmufEQpt1_v5";
const std::string _CT14NNLOFILESTR3="fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQHTp_v5";
const std::string _CT14NNLODESCTXT ="CT14NNLO";


//NP corr fits strings + descriptive text
const std::string HERWGEE4_NPS = "NPC_HerwigEE4C_R4";  //LO
const std::string HERWGEE5_NPS = "NPC_HerwigEE5C_R4";  //LO
const std::string _PYTHIA8_NPS = "NPC_PYTHIA8_R4";      //LO
const std::string _POW_PY8_NPS = "NPC_POWPY8_R4";  //NLO
const std::string POWPY8CT_NPS = "NPC_POWPY8CTEQ_R4";  //NLO

const std::string HERWGEE4_NPS_TXT = "Herwig EE4C R4";  //LO
const std::string HERWGEE5_NPS_TXT = "Herwig EE5C R4";  //LO
const std::string _PYTHIA8_NPS_TXT = "PYTHIA8 R4";      //LO
const std::string _AVERAGE_NPS_TXT = "AVG(PY8+HWG) R4";      //LO
const std::string _POW_PY8_NPS_TXT = "POW+PY8 R4";  //NLO
const std::string POWPY8CT_NPS_TXT = "POW+PY8 CTEQ R4";  //NLO


//descriptive text
std::string ddxsec_yax="#frac{d^{2}#sigma}{dp_{T}dy} [nb/GeV]";


//default file strings
const std::string in_JERFile_default="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots_JERS/output/ak4PF_MCJEC_00eta05.root";
//const std::string NPCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV.root";
const std::string NPCorrFile="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV_fromJohn_refit.root";

const std::string in_NLOFile_default=NNPDFFILESTR;
const std::string in_PY8File_default=_PY8_FILESTR;

std::string ybin_str_arr[]={
  "0.0_to_0.5",
  "0.5_to_1.0",
  "1.0_to_1.5",
  "1.5_to_2.0"
};


std::string etabin_altstr_arr[]={
  "00eta05",
  "15eta10",
  "10eta15",
  "15eta20"
};


bool str_replace(std::string& str, const std::string& from, const std::string& to) {
  //bool str_replace(std::string str, const std::string from, const std::string to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}
