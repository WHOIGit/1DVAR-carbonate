!----------------------------------------------------------------------------
!     CVS:$Id: eco_common.F90,v 1.9 2005/01/19 19:20:41 duse Exp $
!     CVS:$Name:  $
!     Yawei Luo revised, 2007/02/13
!----------------------------------------------------------------------------
module eco_common
  use eco_params, only : numstatevar,numdiagvar

  use types, only : ncvar_attrib,ncf_float

  implicit none
  save


! wnsv - flag sinking rate for each state variable (m/day) 
! 1.0 use wnsvo as sinking rate
! 0.0 disable sinking rate
   double precision, dimension(NumStateVar), parameter :: &
        wnsvflag=(/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                1.0, 1.0, 1.0, 0.0, 0.0, 0.0, &
                0.0, 0.0, 0.0 /)          ! **Required**  

! aeonsv - Aeolean flux flag,  0=off or 1=on
  integer, dimension(NumStateVar), parameter :: &
       aeonsv=0                               ! **Required**   (but not used)

! irradiance
  double precision :: rlt


!-----------------------------------------------------------------------
! ecosystem model parameters and namelist
!-----------------------------------------------------------------------
  double precision :: ae,mu_SP,alpha_SP,a_SP,v_SPn,k_nh4SP,k_no3SP,&
      v_SPp,k_po4SP, zeta, theta, r_excrSP_1,r_excrSP_2,r_pomSP,&
      mu_TR,alpha_TR,a_TR,k_nh4TR,v_TRn,k_no3TR,&
      v_TRp, k_po4TR,mu_pickTRpo4, zeta_nf, &
      r_excrTR_1,r_excrTR_n,r_excrTR_2,r_pomTR,&
      mu_UN,alpha_UN,k_DOM,r_SDOM,mu_BA,& !b_SDONlabi,b_SDOPlabi
      b_BAresp,r_BAadju,r_BAremi,r_BArefr,f_BAslct,r_BAresp_1, r_BAresp_min, &
      r_BAresp_max, r_BAmort, mu_PRT,g_sp,&
      g_ba,r_PRTex,f_exPRTldom,&
      r_PRTresp_1,r_PRTresp_2,r_PRTadju,r_PRTremi,&
      r_pomPRT,mu_MZ,g_prt,g_tr,r_MZex,&
      f_exMZldom,r_MZresp_1,r_MZresp_2,r_MZadju,r_MZremi,&
      r_MZpom,r_MZrefr,r_MZremv,f_HZsdom,f_HZpom,r_SDOMrefr, q_refrDOM_n, q_refrDOM_p, &
      q_POM_n,q_POM_p,r_nitrf, remin_prf_n,remin_prf_p,wnsvo,remin
  
  namelist /ecosys_parms_nml/ &
       ae, &
       mu_SP, &
       alpha_SP, &
       a_SP, &
       v_SPn, &
       k_nh4SP, &
       k_no3SP, &
       v_SPp, &
       k_po4SP, &
       zeta, &
       theta, &
       r_excrSP_1, &
       r_excrSP_2, &
       r_pomSP, &
       mu_TR, &
       alpha_TR, &
       a_TR, &
       v_TRn, &
       k_nh4TR, &
       k_no3TR, &
       v_TRp, &
       k_po4TR, &
       mu_pickTRpo4, &
       zeta_nf, &
       r_excrTR_1, &
       r_excrTR_n, &
       r_excrTR_2, &
       r_pomTR, &
       mu_UN, &
       alpha_UN, &
       k_DOM, &
       !b_SDONlabi, &
       !b_SDOPlabi, &
       r_SDOM , &
       mu_BA, &
       b_BAresp, &
       r_BAadju, &
       r_BAremi, &
       r_BArefr, &
       f_BAslct, &
       r_BAresp_1, &
       r_BAresp_min, &
       r_BAresp_max, &
       r_BAmort, &
       mu_PRT, &
       g_sp, &
       g_ba, &
       r_PRTex, &
       f_exPRTldom, &
       r_PRTresp_1, &
       r_PRTresp_2, &
       r_PRTadju, &
       r_PRTremi, &
       r_pomPRT, &
       mu_MZ, &
       g_prt, &
       g_tr, &
       r_MZex, &
       f_exMZldom, &
       r_MZresp_1, &
       r_MZresp_2, &
       r_MZadju, &
       r_MZremi, &
       r_MZpom, &
       r_MZrefr, &
       r_MZremv, &
       f_HZsdom, &
       f_HZpom, &
       r_SDOMrefr, &
       q_refrDOM_n, &
       q_refrDOM_p, &
       q_POM_n, &
       q_POM_p, &
       r_nitrf, &
       remin_prf_n, &
       remin_prf_p, &
       wnsvo, &
       remin

! unit number for ecosystem parameters I/O
  integer, parameter :: ecopar_unit=45





!-------------------------------------------------------------------------
! NetCDF variable attributes for biological scalars
! ***Required*** for netcdf output of state variables and diagnostics
!-------------------------------------------------------------------------
  type(ncvar_attrib), dimension(NumStateVar) :: bio_ncvaratts=(/ &
       ncvar_attrib('SPc','phytoplankton C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SPn','phytoplankton N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SPp','phytoplankton P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('TRc','Trichodesmium C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('TRn','Trichodesmium N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('TRp','Trichodesmium P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('UNc','Unicellular N2-fixers C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('UNn','Unicellular N2-fixers N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('UNp','Unicellular N2-fixers P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('BAc','bacterial C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('BAn','bacterial N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('BAp','bacterial P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('PRTc','protozoan C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('PRTn','protozoan N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('PRTp','protozoan P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('MZc','microzooplankton C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('MZn','microzooplankton N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('MZp','microzooplankton P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('LDOMc','labile DOM C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('LDOMn','labile DOM N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('LDOMp','labile DOM P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SDOMc','semi-labile DOM C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SDOMn','semi-labile DOM N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SDOMp','semi-labile DOM P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('DETc','detritus C','mmol/m3',-1,ncf_float), &
       ncvar_attrib('DETn','detritus N','mmol/m3',-1,ncf_float), &
       ncvar_attrib('DETp','detritus P','mmol/m3',-1,ncf_float), &
       ncvar_attrib('NH4','ammonium','mmol/m3',-1,ncf_float), &
       ncvar_attrib('NO3','nitrate','mmol/m3',-1,ncf_float), &
       ncvar_attrib('PO4','phosphate','mmol/m3',-1,ncf_float), &
       ncvar_attrib('SPchl','SP chlorophyll a','mg/m3',-1,ncf_float), &
       ncvar_attrib('TRchl','TR chlorophyll a','mg/m3',-1,ncf_float), &
       ncvar_attrib('UNchl','UN chlorophyll','mg/m3',-1,ncf_float) &
       /)

  type(ncvar_attrib), dimension(NumDiagVar) :: diag_ncvaratts=(/ &
       ncvar_attrib('pp','primary productivity','mgC/m3/sec',-1,ncf_float),&
       ncvar_attrib('prBAc','prBAc','mmol C/m3/sec',-1,ncf_float), &
       ncvar_attrib('growSPc','growSPc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growSPnh4','growSPnh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growSPno3','growSPno3','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growSPn','growSPn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growSPp','growSPp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_1c','excrSP_1c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_1n','excrSP_1n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_1p','excrSP_1p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_2c','excrSP_2c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_2n','excrSP_2n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrSP_2p','excrSP_2p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomSPc','pomSPc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomSPn','pomSPn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomSPp','pomSPp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazSPc','grazSPc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazSPn','grazSPn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazSPp','grazSPp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRc','growTRc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRnh4','growTRnh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRno3','growTRno3','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRnf','growTRnf','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRn','growTRn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRpo4','growTRpo4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pickTRpo4','pickTRpo4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growTRp','growTRp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_1c','excrTR_1c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_1n','excrTR_1n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_1p','excrTR_1p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_nh4','excrTR_nh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_2c','excrTR_2c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_2n','excrTR_2n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrTR_2p','excrTR_2p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomTRc','pomTRc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomTRn','pomTRn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomTRp','pomTRp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazTRc','grazTRc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazTRn','grazTRn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazTRp','grazTRp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNc','growUNc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNnh4','growUNnh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNno3','growUNno3','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNnf','growUNnf','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNn','growUNn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growUNp','growUNp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_1c','excrUN_1c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_1n','excrUN_1n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_1p','excrUN_1p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_nh4','excrUN_nh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_2c','excrUN_2c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_2n','excrUN_2n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrUN_2p','excrUN_2p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomUNc','pomUNc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomUNn','pomUNn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomUNp','pomUNp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazUNc','grazUNc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazUNn','grazUNn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazUNp','grazUNp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAldoc','growBAldoc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAldon','growBAldon','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAldop','growBAldop','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAsdoc','growBAsdoc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAsdon','growBAsdon','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAsdop','growBAsdop','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAnh4','growBAnh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAno3','growBAno3','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBApo4','growBApo4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAc','growBAc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAn','growBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growBAp','growBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('respBA','respBA','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrBAc','refrBAc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrBAn','refrBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrBAp','refrBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrBAc','excrBAc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrBAn','excrBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrBAp','excrBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiBAn','remiBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiBAp','remiBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazBAc','grazBAc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazBAn','grazBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazBAp','grazBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('mortBAc','mortBAc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('mortBAn','mortBAn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('mortBAp','mortBAp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('fluxBAnh4','fluxBAnh4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('fluxBApo4','fluxBApo4','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growPRTc','growPRTc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growPRTn','growPRTn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growPRTp','growPRTp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('respPRT','respPRT','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTldomc','excrPRTldomc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTldomn','excrPRTldomn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTldomp','excrPRTldomp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdomc','excrPRTsdomc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdomn','excrPRTsdomn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdomp','excrPRTsdomp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdom2c','excrPRTsdom2c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdom2n','excrPRTsdom2n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrPRTsdom2p','excrPRTsdom2p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiPRTn','remiPRTn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiPRTp','remiPRTp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomPRTc','pomPRTc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomPRTn','pomPRTn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomPRTp','pomPRTp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazPRTc','grazPRTc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazPRTn','grazPRTn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('grazPRTp','grazPRTp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growMZc','growMZc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growMZn','growMZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('growMZp','growMZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('respSP','respSP','mmol C/m3/sec',-1,ncf_float), &
       ncvar_attrib('respTR','respTR','mmol C/m3/sec',-1,ncf_float), &
       ncvar_attrib('respUN','respUN','mmol C/m3/sec',-1,ncf_float), & 
       ncvar_attrib('respMZ','respMZ','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZldomc','excrMZldomc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZldomn','excrMZldomn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZldomp','excrMZldomp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdomc','excrMZsdomc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdomn','excrMZsdomn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdomp','excrMZsdomp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdom2c','excrMZsdom2c','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdom2n','excrMZsdom2n','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrMZsdom2p','excrMZsdom2p','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiMZn','remiMZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiMZp','remiMZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrMZc','refrMZc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrMZn','refrMZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrMZp','refrMZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomMZc','pomMZc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomMZn','pomMZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomMZp','pomMZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remvMZc','remvMZc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remvMZn','remvMZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remvMZp','remvMZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomHZc','pomHZc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomHZn','pomHZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('pomHZp','pomHZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrHZsdomc','excrHZsdomc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrHZsdomn','excrHZsdomn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('excrHZsdomp','excrHZsdomp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiHZn','remiHZn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('remiHZp','remiHZp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('disDETc','disDETc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('disDETn','disDETn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('disDETp','disDETp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('nitrf','nitrf','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrSDOMc','refrSDOMc','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrSDOMn','refrSDOMn','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('refrSDOMp','refrSDOMp','mmol/m3/sec',-1,ncf_float), &
       ncvar_attrib('exportc','C export flux','mmol C/m2/sec',-1,ncf_float), &
       ncvar_attrib('exportn','N export flux','mmol C/m2/sec',-1,ncf_float), &
       ncvar_attrib('exportp','P export flux','mmol C/m2/sec',-1,ncf_float) /)
contains



! ***Required*** 
  subroutine read_eco_params(bioparams_default,bioparams)

    use common_mod, only : ecopar_fname
    use eco_params, only : iae, imu_SP,ialpha_SP,ia_SP,iv_SPn,ik_nh4SP,ik_no3SP,&
          iv_SPp,ik_po4SP,izeta,itheta,ir_excrSP_1,ir_excrSP_2,ir_pomSP,&
	      imu_TR,ialpha_TR,ia_TR,iv_TRn, ik_nh4TR,ik_no3TR,iv_TRp,ik_po4TR,&
	      imu_pickTRpo4, izeta_nf,&
          ir_excrTR_1,ir_excrTR_n,ir_excrTR_2,ir_pomTR,imu_UN,&
          ialpha_UN,ik_DOM,ir_SDOM,imu_BA,& !ib_SDONlabi,ib_SDOPlabi,
          ib_BAresp,ir_BAadju,ir_BAremi,ir_BArefr,&
          if_BAslct,ir_BAresp_1, ir_BAresp_min, ir_BAresp_max, ir_BAmort, & 
          imu_PRT,ig_sp,ig_ba,&
          ir_PRTex,if_exPRTldom,ir_PRTresp_1,&
          ir_PRTresp_2,ir_PRTadju,ir_PRTremi,&
          ir_pomPRT,imu_MZ,ig_prt,ig_tr,&
          ir_MZex,if_exMZldom,ir_MZresp_1,ir_MZresp_2,&
          ir_MZadju,ir_MZremi,ir_MZpom,ir_MZrefr,ir_MZremv,if_HZsdom,if_HZpom,ir_SDOMrefr,&
          iq_refrDOM_n, iq_refrDOM_p,iq_POM_n,iq_POM_p,&
          ir_nitrf,iremin_prf_n,iremin_prf_p,iwnsvo,iremin
    implicit none

!-----------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------
    double precision, dimension(:), intent(out) :: &
         bioparams,bioparams_default

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer, parameter :: ecopar_unit=45

!-----------------------------------------------------------------------
!   default namelist settings
!-----------------------------------------------------------------------

 AE = 4000.0
 MU_SP = 1.7096560001498988
 ALPHA_SP = 5.51273526304891
 A_SP = 5.0E-4
 V_SPN = 0.3
 K_NH4SP = 0.01
 K_NO3SP = 0.1
 V_SPP = 0.010133492869505855
 K_PO4SP = 4.8E-3
 ZETA = 2.0
 THETA = 2.1288532650970557
 R_EXCRSP_1 = 0.05
 R_EXCRSP_2 = 0.05
 R_POMSP = 0.038
 MU_TR = 3.0
 ALPHA_TR = 0.25
 A_TR = 2.0E-3
 V_TRN = 0.4
 K_NH4TR = 0.05
 K_NO3TR = 0.5
 V_TRP = 0.02
 K_PO4TR = 0.047
 MU_PICKTRPO4 = 1.0
 ZETA_NF = 2.5
 R_EXCRTR_1 = 0.02
 R_EXCRTR_N = 0.36
 R_EXCRTR_2 = 0.05
 R_POMTR = 0.08
 MU_UN = 1.6
 ALPHA_UN = 0.5
 K_DOM = 0.5
 R_SDOM = 9.586773754551453E-4
 MU_BA = 2.0
 B_BARESP = 0.28
 R_BAADJU = 2.0
 R_BAREMI = 6.0
 R_BAREFR = 0.018
 F_BASLCT = 0.2
 R_BARESP_1 = 0.01
 R_BARESP_MIN = 0.35
 R_BARESP_MAX = 0.6734616399721387
 R_BAMORT = 0.049
 MU_PRT = 4.0
 G_SP = 0.7512749623651648
 G_BA = 1.4299980380426615
 R_PRTEX = 0.2
 F_EXPRTLDOM = 0.9
 R_PRTRESP_1 = 0.01
 R_PRTRESP_2 = 0.7584721354023397
 R_PRTADJU = 2.0
 R_PRTREMI = 4.7
 R_POMPRT = 0.027
 MU_MZ = 1.3
 G_PRT = 1.4
 G_TR = 0.20072614502968935
 R_MZEX = 0.3
 F_EXMZLDOM = 0.75
 R_MZRESP_1 = 0.03
 R_MZRESP_2 = 0.22
 R_MZADJU = 2.0
 R_MZREMI = 4.0
 R_MZPOM = 0.15
 R_MZREFR = 0.02
 R_MZREMV = 0.81
 F_HZSDOM = 0.1
 F_HZPOM = 0.14
 R_SDOMREFR = 9.0E-4
 Q_REFRDOM_N = 0.05
 Q_REFRDOM_P = 7.0E-4
 Q_POM_N = 0.12
 Q_POM_P = 4.5E-3
 R_NITRF = 0.1
 REMIN_PRF_N = 1.1
 REMIN_PRF_P = 4.0
 WNSVO = 5.262475444554562
 REMIN = 0.011
 
!-----------------------------------------------------------------------
! set default parameter values
!-----------------------------------------------------------------------
    bioparams_default(iae          )=ae
    bioparams_default(imu_SP          )=mu_SP
    bioparams_default(ialpha_SP       )=alpha_SP
    bioparams_default(ia_SP       )=a_SP
    bioparams_default(iv_SPn       )=v_SPn
    bioparams_default(ik_nh4SP        )=k_nh4SP
    bioparams_default(ik_no3SP        )=k_no3SP
    bioparams_default(iv_SPp       )=v_SPp
    bioparams_default(ik_po4SP        )=k_po4SP
    bioparams_default(izeta       )=zeta
    bioparams_default(itheta       )=theta
    bioparams_default(ir_excrSP_1     )=r_excrSP_1
    bioparams_default(ir_excrSP_2     )=r_excrSP_2
    bioparams_default(ir_pomSP        )=r_pomSP
    bioparams_default(imu_TR          )=mu_TR
    bioparams_default(ialpha_TR       )=alpha_TR
    bioparams_default(ia_TR       )=a_TR
    bioparams_default(iv_TRn       )=v_TRn
    bioparams_default(ik_nh4TR        )=k_nh4TR
    bioparams_default(ik_no3TR        )=k_no3TR
    bioparams_default(iv_TRp       )=v_TRp
    bioparams_default(ik_po4TR        )=k_po4TR
    bioparams_default(imu_pickTRpo4) = mu_pickTRpo4
    bioparams_default(izeta_nf       )=zeta_nf
    bioparams_default(ir_excrTR_1     )=r_excrTR_1
    bioparams_default(ir_excrTR_n   )=r_excrTR_n
    bioparams_default(ir_excrTR_2     )=r_excrTR_2
    bioparams_default(ir_pomTR        )=r_pomTR
    bioparams_default(imu_UN          )=mu_UN
    bioparams_default(ialpha_UN       )=alpha_UN
    bioparams_default(ik_DOM        )=k_DOM
    !bioparams_default(ib_SDONlabi     )=b_SDONlabi
    !bioparams_default(ib_SDOPlabi     )=b_SDOPlabi
    bioparams_default(ir_SDOM     )=r_SDOM
    bioparams_default(imu_BA          )=mu_BA
    bioparams_default(ib_BAresp       )=b_BAresp
    bioparams_default(ir_BAadju       )=r_BAadju
    bioparams_default(ir_BAremi       )=r_BAremi
    bioparams_default(ir_BArefr       )=r_BArefr
    bioparams_default(if_BAslct) = f_BAslct
    bioparams_default(ir_BAresp_1) = r_BAresp_1
    bioparams_default(ir_BAresp_min) = r_BAresp_min
    bioparams_default(ir_BAresp_max) = r_BAresp_max
    bioparams_default(ir_BAmort) = r_BAmort
    bioparams_default(imu_PRT       )=mu_PRT
    bioparams_default(ig_sp        )=g_sp
    bioparams_default(ig_ba        )=g_ba
    bioparams_default(ir_PRTex        )=r_PRTex
    bioparams_default(if_exPRTldom    )=f_exPRTldom
    bioparams_default(ir_PRTresp_1    )=r_PRTresp_1
    bioparams_default(ir_PRTresp_2    )=r_PRTresp_2
    bioparams_default(ir_PRTadju      )=r_PRTadju
    bioparams_default(ir_PRTremi      )=r_PRTremi
    bioparams_default(ir_pomPRT       )=r_pomPRT
    bioparams_default(imu_MZ       )=mu_MZ
    bioparams_default(ig_prt        )=g_prt
    bioparams_default(ig_tr         )=g_tr
    bioparams_default(ir_MZex         )=r_MZex
    bioparams_default(if_exMZldom     )=f_exMZldom
    bioparams_default(ir_MZresp_1     )=r_MZresp_1
    bioparams_default(ir_MZresp_2     )=r_MZresp_2
    bioparams_default(ir_MZadju       )=r_MZadju
    bioparams_default(ir_MZremi       )=r_MZremi
    bioparams_default(ir_MZpom        )=r_MZpom
    bioparams_default(ir_MZrefr       )=r_MZrefr
    bioparams_default(ir_MZremv       )=r_MZremv
    bioparams_default(if_HZsdom       )=f_HZsdom
    bioparams_default(if_HZpom        )=f_HZpom
    bioparams_default(ir_SDOMrefr     )=r_SDOMrefr
    bioparams_default(iq_refrDOM_n) = q_refrDOM_n
    bioparams_default(iq_refrDOM_p) = q_refrDOM_p
    bioparams_default(iq_POM_n        )=q_POM_n
    bioparams_default(iq_POM_p        )=q_POM_p
    bioparams_default(ir_nitrf        )=r_nitrf
    bioparams_default(iremin_prf_n    )=remin_prf_n
    bioparams_default(iremin_prf_p    )=remin_prf_p
    bioparams_default(iwnsvo         )=wnsvo
    bioparams_default(iremin          )=remin


!-----------------------------------------------------------------------
!   read in parameter values
!-----------------------------------------------------------------------
    open(unit=ecopar_unit, file=trim(ecopar_fname), status='old')
    read(unit=ecopar_unit, nml=ecosys_parms_nml)
    close(unit=ecopar_unit)

    bioparams(iae          )=ae
    bioparams(imu_SP          )=mu_SP
    bioparams(ialpha_SP       )=alpha_SP
    bioparams(ia_SP       )=a_SP
    bioparams(iv_SPn       )=v_SPn
    bioparams(ik_nh4SP        )=k_nh4SP
    bioparams(ik_no3SP        )=k_no3SP
    bioparams(iv_SPp       )=v_SPp
    bioparams(ik_po4SP        )=k_po4SP
    bioparams(izeta       )=zeta
    bioparams(itheta       )=theta
    bioparams(ir_excrSP_1     )=r_excrSP_1
    bioparams(ir_excrSP_2     )=r_excrSP_2
    bioparams(ir_pomSP        )=r_pomSP
    bioparams(imu_TR          )=mu_TR
    bioparams(ialpha_TR       )=alpha_TR
    bioparams(ia_TR       )=a_TR
    bioparams(iv_TRn       )=v_TRn
    bioparams(ik_nh4TR        )=k_nh4TR
    bioparams(ik_no3TR        )=k_no3TR
    bioparams(iv_TRp       )=v_TRp
    bioparams(ik_po4TR        )=k_po4TR
    bioparams(imu_pickTRpo4) = mu_pickTRpo4
    bioparams(izeta_nf       )=zeta_nf
    bioparams(ir_excrTR_1     )=r_excrTR_1
    bioparams(ir_excrTR_n   )=r_excrTR_n
    bioparams(ir_excrTR_2     )=r_excrTR_2
    bioparams(ir_pomTR        )=r_pomTR
    bioparams(imu_UN          )=mu_UN
    bioparams(ialpha_UN       )=alpha_UN
    bioparams(ik_DOM        )=k_DOM
    !bioparams(ib_SDONlabi     )=b_SDONlabi
    !bioparams(ib_SDOPlabi     )=b_SDOPlabi
    bioparams(ir_SDOM     )=r_SDOM
    bioparams(imu_BA          )=mu_BA
    bioparams(ib_BAresp       )=b_BAresp
    bioparams(ir_BAadju       )=r_BAadju
    bioparams(ir_BAremi       )=r_BAremi
    bioparams(ir_BArefr       )=r_BArefr
    bioparams(if_BAslct) = f_BAslct
    bioparams(ir_BAresp_1) = r_BAresp_1
    bioparams(ir_BAresp_min) = r_BAresp_min
    bioparams(ir_BAresp_max) = r_BAresp_max
    bioparams(ir_BAmort) = r_BAmort
    bioparams(imu_PRT       )=mu_PRT
    bioparams(ig_sp        )=g_sp
    bioparams(ig_ba        )=g_ba
    bioparams(ir_PRTex        )=r_PRTex
    bioparams(if_exPRTldom    )=f_exPRTldom
    bioparams(ir_PRTresp_1    )=r_PRTresp_1
    bioparams(ir_PRTresp_2    )=r_PRTresp_2
    bioparams(ir_PRTadju      )=r_PRTadju
    bioparams(ir_PRTremi      )=r_PRTremi
    bioparams(ir_pomPRT       )=r_pomPRT
    bioparams(imu_MZ       )=mu_MZ
    bioparams(ig_prt        )=g_prt
    bioparams(ig_tr         )=g_tr
    bioparams(ir_MZex         )=r_MZex
    bioparams(if_exMZldom     )=f_exMZldom
    bioparams(ir_MZresp_1     )=r_MZresp_1
    bioparams(ir_MZresp_2     )=r_MZresp_2
    bioparams(ir_MZadju       )=r_MZadju
    bioparams(ir_MZremi       )=r_MZremi
    bioparams(ir_MZpom        )=r_MZpom
    bioparams(ir_MZrefr       )=r_MZrefr
    bioparams(ir_MZremv       )=r_MZremv
    bioparams(if_HZsdom       )=f_HZsdom
    bioparams(if_HZpom        )=f_HZpom
    bioparams(ir_SDOMrefr     )=r_SDOMrefr
    bioparams(iq_refrDOM_n) = q_refrDOM_n
    bioparams(iq_refrDOM_p) = q_refrDOM_p
    bioparams(iq_POM_n        )=q_POM_n
    bioparams(iq_POM_p        )=q_POM_p
    bioparams(ir_nitrf        )=r_nitrf
    bioparams(iremin_prf_n    )=remin_prf_n
    bioparams(iremin_prf_p    )=remin_prf_p
    bioparams(iwnsvo          )=wnsvo
    bioparams(iremin          )=remin
    
  end subroutine read_eco_params

  



  subroutine write_eco_params(bioparams,filename)
    use eco_params
    implicit none

!-----------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------
    double precision, dimension(:), intent(in) :: bioparams
    character(len=*) :: filename

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    ae          = bioparams(iae         )
    mu_SP          = bioparams(imu_SP          )
    alpha_SP       = bioparams(ialpha_SP       )
    a_SP       = bioparams(ia_SP       )
    v_SPn        = bioparams(iv_SPn        )
    k_nh4SP        = bioparams(ik_nh4SP        )
    k_no3SP        = bioparams(ik_no3SP        )
    v_SPp           = bioparams(iv_SPp       )
    k_po4SP        = bioparams(ik_po4SP        )
    zeta              = bioparams(izeta       )
    theta            = bioparams(itheta       )
    r_excrSP_1     = bioparams(ir_excrSP_1     )
    r_excrSP_2     = bioparams(ir_excrSP_2     )
    r_pomSP        = bioparams(ir_pomSP        )
    mu_TR          = bioparams(imu_TR          )
    alpha_TR       = bioparams(ialpha_TR       )
    a_TR           = bioparams(ia_TR       )
    v_TRn        = bioparams(iv_TRn        )
    k_nh4TR        = bioparams(ik_nh4TR        )
    k_no3TR        = bioparams(ik_no3TR        )
    v_TRp        = bioparams(iv_TRp        )
    k_po4TR        = bioparams(ik_po4TR        )
    mu_pickTRpo4 = bioparams(imu_pickTRpo4)
    zeta_nf        = bioparams(izeta_nf        )
    r_excrTR_1     = bioparams(ir_excrTR_1     )
    r_excrTR_n   = bioparams(ir_excrTR_n   )
    r_excrTR_2     = bioparams(ir_excrTR_2     )
    r_pomTR        = bioparams(ir_pomTR        )
    mu_UN          = bioparams(imu_UN          )
    alpha_UN       = bioparams(ialpha_UN       )
    k_DOM        = bioparams(ik_DOM        )
    !b_SDONlabi     = bioparams(ib_SDONlabi     )
    !b_SDOPlabi     = bioparams(ib_SDOPlabi     )
    r_SDOM     = bioparams(ir_SDOM     )
    mu_BA          = bioparams(imu_BA          )
    b_BAresp       = bioparams(ib_BAresp       )
    r_BAadju       = bioparams(ir_BAadju       )
    r_BAremi       = bioparams(ir_BAremi       )
    r_BArefr       = bioparams(ir_BArefr       )
    f_BAslct = bioparams(if_BAslct)
    r_BAresp_1 = bioparams(ir_BAresp_1)
    r_BAresp_min = bioparams(ir_BAresp_min)
    r_BAresp_max = bioparams(ir_BAresp_max)
    r_BAmort = bioparams(ir_BAmort)    
    mu_PRT       = bioparams(imu_PRT       )
    g_sp        = bioparams(ig_sp        )
    g_ba        = bioparams(ig_ba        )
    r_PRTex        = bioparams(ir_PRTex        )
    f_exPRTldom    = bioparams(if_exPRTldom    )
    r_PRTresp_1    = bioparams(ir_PRTresp_1    )
    r_PRTresp_2    = bioparams(ir_PRTresp_2    )
    r_PRTadju      = bioparams(ir_PRTadju      )
    r_PRTremi      = bioparams(ir_PRTremi      )
    r_pomPRT       = bioparams(ir_pomPRT       )
    mu_MZ       = bioparams(imu_MZ       )
    g_prt        = bioparams(ig_prt        )
    g_tr         = bioparams(ig_tr         )
    r_MZex         = bioparams(ir_MZex         )
    f_exMZldom     = bioparams(if_exMZldom     )
    r_MZresp_1     = bioparams(ir_MZresp_1     )
    r_MZresp_2     = bioparams(ir_MZresp_2     )
    r_MZadju       = bioparams(ir_MZadju       )
    r_MZremi       = bioparams(ir_MZremi       )
    r_MZpom        = bioparams(ir_MZpom        )
    r_MZrefr       = bioparams(ir_MZrefr       )
    r_MZremv       = bioparams(ir_MZremv       )
    f_HZsdom       = bioparams(if_HZsdom       )
    f_HZpom        = bioparams(if_HZpom        )
    r_SDOMrefr     = bioparams(ir_SDOMrefr     )
    q_refrDOM_n = bioparams(iq_refrDOM_n) 
    q_refrDOM_p = bioparams(iq_refrDOM_p)
    q_POM_n        = bioparams(iq_POM_n        )
    q_POM_p        = bioparams(iq_POM_p        )
    r_nitrf        = bioparams(ir_nitrf        )
    remin_prf_n    = bioparams(iremin_prf_n    )
    remin_prf_p    = bioparams(iremin_prf_p    )
    wnsvo          = bioparams(iwnsvo          )
    remin          = bioparams(iremin          )

    open(unit=ecopar_unit, file=trim(filename), status='replace')
    write(unit=ecopar_unit, nml=ecosys_parms_nml)
    close(unit=ecopar_unit)
    
  end subroutine write_eco_params




end module eco_common
