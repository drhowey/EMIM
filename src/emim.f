cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c EMIM, version 3.22
c Copyright 2012-2015,
c Heather Cordell
c Institute of Human Genetics, Newcastle University
c
c heather.cordell@ncl.ac.uk
c http://www.staff.ncl.ac.uk/heather.cordell/
c
c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program.  If not, see <http://www.gnu.org/licenses/>.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c    UPDATED VERSION USES 9 MATING TYPES (CPG) IF DESIRED

c
c     Prog to analyses ma/child duos and case/parent trios
c     a la Weinberg
c     But using multinomial max lik
c     rather than log linear model

c     Can also HWE and random mating in order
c     to improve identifiability
c
c     Actually use ln of RR params  R1, R2, S1, S2, Im, Ip, 
c     gam11, gam12, gam21, gam22, A1, mew1-6

c
c     Compile using the following:

c     gfortran emim.f maxfun.f -o emim




      module shared_var
             implicit none

        double precision casetriocount(17)
        double precision casemacount(9)
        double precision casepacount(9)
        double precision propNoPhTRIO, propPhTRIO 
        double precision cntNoPhTRIO, cntPhTRIO
        double precision propNoPhCAMD, propPhCAMD
        double precision cntNoPhCAMD, cntPhCAMD
        double precision propNoPhCAFD, propPhCAFD 
        double precision cntNoPhCAFD, cntPhCAFD

        integer conmacount(7), conpacount(7),
     +casecount(3),  conparcount(9), concount(3),
     +HWEIND,
     +ALSYMIND, fixAIND, CPGIND,
     +casetrioIND, casemaduoIND, casepaduoIND, caseIND,
     +maofcaseIND,paofcaseIND,parofcaseIND,
     +maofcasecount(3),paofcasecount(3),parofcasecount(9),
     +conparIND,conmaduoIND, conpaduoIND, conIND,
     +INDR1eqR2, INDR1sqeqR2, INDS1eqS2, INDS1sqeqS2,
     +INDgameq, INDweinIm, INDweinIp, INDsinsgam01, 
     +INDsinsgam21, INDpalmer, INDjc, INDIm, INDIp

      save
      end module shared_var


       PROGRAM EMIM
       use shared_var
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
       PARAMETER (NP=45,NPV=45)

c
c      COMMON /MAXFLB/ LABEL(NP)
      COMMON /MAXFOP/ IOUT,IDET,LPRT
      COMMON /MAXF1/ THIN(NP),THL(NP),THU(NP),STPIN(NP),EPSD,YOTA,EPST,
     $               EPSC1,EPSC2,EPSC3,ISTIN(NP),NT,MAXIT,METHOD,IXVC,
     $               IHIT
      COMMON /MAXF2/ THPR(NP),CTH(NP),STP(NP),G(NPV),H(NPV,NPV),
     $               V(NPV,NPV),AV(NP,NP),STDE(NP),ULT(NPV,NPV),
     $               DIAG(NPV),PDIR(NPV),ERM,FPR,FCH,DIFMAX,GTG,PTG,
     $               TSTEP,IST(NP),NE,ND,NI,NB,NV,IMPBND,IT,NSURF2,IGFL,
     $               IVFL,IGAGE,IVAGE,IDIF
	COMMON /XAMP2/ Y(1000), EAIJ(1000,5),PREVTR(11),RND2,CONST,N,K,KPU,
     $                                    KX
      COMMON /MKUND/ KP2,KK
C
C
C Local variables
C

      CHARACTER char1*4, char2*6, 
     +c01*5, c02*8, c03*8, c04*8,
     +c11*5, c12*8, c13*8, c14*8,
     +c21*5, c22*8, c23*8, c24*8,    
     +c31*5, c32*8, c33*8, c34*8,
     +c41*5, c42*8, c43*8, c44*8,
     +c51*5, c52*8, c53*8, c54*8, 
     +c61*8, c62*11, c63*11, c64*11,    
     +c71*8, c72*11, c73*11, c74*11,
     +c81*8, c82*11, c83*11, c84*11,
     +c91*8, c92*11, c93*11, c94*11, 
     +char3*10, char4*10, char5*10, charfr*5, char6*9
     

      INTEGER LFL, NFE, i, j, k, snp, nsnps,
     +INDR1, INDR2, INDS1, INDS2, 
     +INDgam11, INDgam12, INDgam21,
     +INDgam22, IND(10), numestparam,
     +sumall,sumcasetrio,sumcasemaduo,sumcasepaduo,summaofcase,
     +sumpaofcase,sumcase,sumconpar,sumconmaduo,sumconpaduo,sumcon
     +IOstatus, warning

      DIMENSION THETA(20)

      double precision R1, R2, S1, S2, Im, Ip, A1, 
     +gam11, gam12, gam21, gam22,
     +partotcount, totcount, fixA1,fixA2,
     +lnlikfull, lnliknull, start(17), mew1, mew2, 
     +mew3, mew4, mew5, mew6, sd(10), mew1dash, mew2dash, 
     +mew3dash, mew4dash, mew5dash, mew6dash, snpno, 
     +oldsnpno, totvar, cpgstart(20), mew2a, mew2b,
     +mew3a, mew3b, mew5a, mew5b

c        real 


      EXTERNAL FUNCTION
      EXTERNAL MAXFUN
      EXTERNAL DEPAR1
      
      INTEGER lenMkrFile
      CHARACTER markersFile*200, summaryFile*200,
     +resultsFile*200, caseParentTriosFile*200,
     +caseParentsFile*200,caseFatherDuosFile*200,
     +caseFathersFile*200,caseMotherDuosFile*200,
     +caseMothersFile*200,casesFile*200,
     +conFatherDuosFile*200,conMotherDuosFile*200,
     +conParentsFile*200,consFile*200,
     +snpNoStr*20, inputDir*200, outputDir*200
      INTEGER noOfArgs
     
     
      noOfArgs = IARGC()

      IF(noOfArgs .gt. 0) THEN
      
       CALL GETARG(1, snpNoStr)

       markersFile='emimmarkers'//TRIM(snpNoStr)//'.dat'
       summaryFile='emimsummary'//TRIM(snpNoStr)//'.out'
       resultsFile='emimresults'//TRIM(snpNoStr)//'.out'
       caseParentTriosFile='caseparenttrios'//TRIM(snpNoStr)//'.dat'
       caseparentsFile='caseparents'//TRIM(snpNoStr)//'.dat'
       caseFatherDuosFile='casefatherduos'//TRIM(snpNoStr)//'.dat'
       caseFathersFile='casefathers'//TRIM(snpNoStr)//'.dat'
       caseMotherDuosFile='casemotherduos'//TRIM(snpNoStr)//'.dat'
       caseMothersFile='casemothers'//TRIM(snpNoStr)//'.dat'
       casesFile='cases'//TRIM(snpNoStr)//'.dat'
       conFatherDuosFile='confatherduos'//TRIM(snpNoStr)//'.dat'
       conMotherDuosFile='conmotherduos'//TRIM(snpNoStr)//'.dat'
       conParentsFile='conparents'//TRIM(snpNoStr)//'.dat'
       consFile='cons'//TRIM(snpNoStr)//'.dat'
      ELSE
       markersFile='emimmarkers.dat'
       summaryFile='emimsummary.out'
       resultsFile='emimresults.out'
       caseParentTriosFile='caseparenttrios.dat'
       caseparentsFile='caseparents.dat'
       caseFatherDuosFile='casefatherduos.dat'
       caseFathersFile='casefathers.dat'
       caseMotherDuosFile='casemotherduos.dat'
       caseMothersFile='casemothers.dat'
       casesFile='cases.dat'
       conFatherDuosFile='confatherduos.dat'
       conMotherDuosFile='conmotherduos.dat'
       conParentsFile='conparents.dat'
       consFile='cons.dat'
      ENDIF

      IF(noOfArgs .eq. 3) THEN
      
       CALL GETARG(2, inputDir)
       CALL GETARG(3, outputDir)

       markersFile=TRIM(inputDir)//markersFile
       summaryFile=TRIM(outputDir)//summaryFile
       resultsFile=TRIM(outputDir)//resultsFile
       caseParentTriosFile=TRIM(inputDir)//caseParentTriosFile
       caseparentsFile=TRIM(inputDir)//caseparentsFile
       caseFatherDuosFile=TRIM(inputDir)//caseFatherDuosFile
       caseFathersFile=TRIM(inputDir)//caseFathersFile
       caseMotherDuosFile=TRIM(inputDir)//caseMotherDuosFile
       caseMothersFile=TRIM(inputDir)//caseMothersFile
       casesFile=TRIM(inputDir)//casesFile
       conFatherDuosFile=TRIM(inputDir)//conFatherDuosFile
       conMotherDuosFile=TRIM(inputDir)//conMotherDuosFile
       conParentsFile=TRIM(inputDir)//conParentsFile
       consFile=TRIM(inputDir)//consFile
      ENDIF

     
      open (unit=3, file=markersFile, status='old')
      open (unit=4, file=resultsFile, status='unknown')
      open (unit=17, file=summaryFile, status='unknown')
      open (unit=5, file='emimparams.dat', status='old')

      write(4,*) "EMIM version 3.22"
      write(4,*)

      warning=0

      sumall=0
      sumcasetrio=0
      sumcasemaduo=0
      sumcasepaduo=0
      summaofcase=0
      sumpaofcase=0
      sumcase=0
      sumconpar=0
      sumconmaduo=0
      sumconpaduo=0
      sumcon=0   


      char1 = 'snp '
      char2 = 'snpID '
      c01='lnR1 '
      c02='sd_lnR1 '
      c03='cl_lnR1 '
      c04='cu_lnR1 '
      c11='lnR2 '
      c12='sd_lnR2 '
      c13='cl_lnR2 '
      c14='cu_lnR2 '
      c21='lnS1 '
      c22='sd_lnS1 '
      c23='cl_lnS1 '
      c24='cu_lnS1 '
      c31='lnS2 '
      c32='sd_lnS2 '
      c33='cl_lnS2 '
      c34='cu_lnS2 '
      c41='lnIm '
      c42='sd_lnIm '
      c43='cl_lnIm '
      c44='cu_lnIm '
      c51='lnIp '
      c52='sd_lnIp '
      c53='cl_lnIp '
      c54='cu_lnIp '
      c61='lngam11 '
      c62='sd_lngam11 '
      c63='cl_lngam11 '
      c64='cu_lngam11 '
      c71='lngam12 '
      c72='sd_lngam12 '
      c73='cl_lngam12 '
      c74='cu_lngam12 '
      c81='lngam21 '
      c82='sd_lngam21 '
      c83='cl_lngam21 '
      c84='cu_lngam21 '
      c91='lngam22 '
      c92='sd_lngam22 '
      c93='cl_lngam22 '
      c94='cu_lngam22 '
      char3 = 'lnliknull '
      char4 = 'lnlikfull '
      char5 = 'twicediff '
      char6 = 'WARNINGS '
      charfr= 'freq '


c     NEED TO ALTER THIS IF USING SINSHEIMER/LI/PALMER PARAMS


c      write(17,7777) char1, char2, charfr,
c     +c01, c02, c03, c04,
c     +c11, c12, c13, c14,
c     +c21, c22, c23, c24,    
c     +c31, c32, c33, c34,
c     +c41, c42, c43, c44,
c     +c51, c52, c53, c54, 
c     +c61, c62, c63, c64,    
c     +c71, c72, c73, c74,
c     +c81, c82, c83, c84,
c     +c91, c92, c93, c94, 
c     +char3, char4, char5
c
c 7777 format(A4,A6,A5,6(A5,A8,A8,A8),4(A8,A11,A11,A11) 3A10)
c
c
c      write(17,*) "snp snpID lnR1 sd_lnR1 cl_lnR1 cu_lnR1
c     + lnR2 sd_lnR2 cl_lnR2 cu_lnR2
c     + lnS1 sd_lnS1 cl_lnS1 cu_lnS1
c     + lnS2 sd_lnS2 cl_lnS2 cu_lnS2
c     + lnIm sd_lnIm cl_lnIm cu_lnIm
c     + lnIp sd_lnIp cl_lnIp cu_lnIp
c     + lngam11 sd_lngam11 cl_lngam11 cu_lngam11
c     + lngam12 sd_lngam12 cl_lngam12 cu_lngam12
c     + lngam21 sd_lngam21 cl_lngam21 cu_lngam21
c     + lngam22 sd_lngam22 cl_lngam22 cu_lngam22
c     + lnliknull lnlikfull twice(lnlikfull-lnliknull)"

      read(5,*)
      read(5,*) casetrioIND 
      read(5,*) parofcaseIND
      read(5,*) casemaduoIND
      read(5,*) casepaduoIND 
      read(5,*) maofcaseIND
      read(5,*) paofcaseIND
      read(5,*) caseIND
      read(5,*) conparIND 
      read(5,*) conmaduoIND 
      read(5,*) conpaduoIND 
      read(5,*) conIND
      read(5,*) nsnps
    
      read(5,*) 
      read(5,*) fixAIND
      read(5,*) HWEIND
      read(5,*) ALSYMIND
      read(5,*) CPGIND

      read(5,*) INDR1
      read(5,*) INDR2
      read(5,*) INDR1eqR2
      read(5,*) INDR1sqeqR2
      read(5,*) INDS1
      read(5,*) INDS2     
      read(5,*) INDS1eqS2 
      read(5,*) INDS1sqeqS2
      read(5,*) INDIm
      read(5,*) INDIp 
      read(5,*) INDgam11
      read(5,*) INDgam12
      read(5,*) INDgam21
      read(5,*) INDgam22
      read(5,*) INDgameq

c      write(4,*) "INDgam=",INDgam11,INDgam12,INDgam21,INDgam22

      read(5,*)
      read(5,*) INDweinIm
      read(5,*) INDweinIp
      read(5,*) INDsinsgam01
      read(5,*) INDsinsgam21
      read(5,*) INDpalmer
      read(5,*) INDjc


      if ((INDweinIm.eq.1).and.(INDIm.eq.1)) then
         write(4,*) "CANNOT ESTIMATE Im AND WEINBERG (1999b) Im"
         write(4,*) "MUST CHOOSE WHICH PARAMETRIZATION TO USE"
         write(4,*) "ESTIMATING Im ONLY"
         warning=1
         write(4,*)
         INDweinIm=0
         end if


c         if ((INDR1.eq.1).or.(INDR2.eq.1).or.(INDR1eqR2.eq.1)
c     +.or.(INDR1sqeqR2.eq.1).or.(INDgam11.eq.1).or.(INDgam12.eq.1)
c     +.or.(INDgam21.eq.1).or.(INDgam22.eq.1).or.(INDgameq.eq.1)
c     +.or.(INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1).or.
c     +(INDpalmer.eq.1).or.(INDjc.eq.1)) then
c
c     DONT ALLOW THIS EVEN IT IS THEORETICALLY ESTIMABLE
c     DOESNT MAKE SCIENTIFIC SENSE AND I HAVENT PROGRAMMED
c     MAXFUN FUNCTION TO ALLOW IT

      if ((INDweinIm.eq.1).and.(INDIp.eq.1)) then
         write(4,*) "CANNOT ESTIMATE Ip AND WEINBERG (1999b) Im"
         write(4,*) "MUST CHOOSE WHICH PARAMETRIZATION TO USE"
         write(4,*) "ESTIMATING Ip ONLY"
         warning=1
         write(4,*)
         INDweinIm=0
         end if

c         end if


      if ((INDweinIp.eq.1).and.(INDIp.eq.1)) then
         write(4,*) "CANNOT ESTIMATE Ip AND WEINBERG (1999b) Ip"
         write(4,*) "MUST CHOOSE WHICH PARAMETRIZATION TO USE"
         write(4,*) "ESTIMATING Ip ONLY"
         warning=1
         write(4,*)
         INDweinIp=0
         end if


c         if ((INDR1.eq.1).or.(INDR2.eq.1).or.(INDR1eqR2.eq.1)
c     +.or.(INDR1sqeqR2.eq.1).or.(INDgam11.eq.1).or.(INDgam12.eq.1)
c     +.or.(INDgam21.eq.1).or.(INDgam22.eq.1).or.(INDgameq.eq.1)
c     +.or.(INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1).or.
c     +(INDpalmer.eq.1).or.(INDjc.eq.1)) then
c
c     DONT ALLOW THIS EVEN IT IS THEORETICALLY ESTIMABLE
c     DOESNT MAKE SCIENTIFIC SENSE AND I HAVENT PROGRAMMED
c     MAXFUN FUNCTION TO ALLOW IT

      if ((INDweinIp.eq.1).and.(INDIm.eq.1)) then
         write(4,*) "CANNOT ESTIMATE Im AND WEINBERG (1999b) Ip"
         write(4,*) "MUST CHOOSE WHICH PARAMETRIZATION TO USE"
         write(4,*) "ESTIMATING Im ONLY"
         warning=1
         write(4,*)
         INDweinIp=0
         end if

c         end if


c     RE-USE TR(5) and TR(6) FOR WEINBERG PARAMETERIZATION

         if (INDweinIm.eq.1) INDIm=INDweinIm
         if (INDweinIp.eq.1) INDIp=INDweinIp


      if (((INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1)).and.
     +((INDgam11.eq.1).or.(INDgam12.eq.1).or.(INDgam21.eq.1).or.
     +(INDgam22.eq.1).or.(INDgameq.eq.1))) then

         write(4,*) "CANNOT ESTIMATE BOTH AINSWORTH (2010) AND"
         write(4,*) "SINSHEIMER (2003) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING AINSWORTH (2010) gammas ONLY"
         warning=1
         write(4,*)
         INDsinsgam01=0
         INDsinsgam21=0
         end if

      if (((INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1)).and.
     +(INDpalmer.eq.1)) then

         write(4,*) "CANNOT ESTIMATE BOTH PALMER (2006) AND"
         write(4,*) "SINSHEIMER (2003) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING PALMER (2006) ONLY"
         warning=1
         write(4,*)
         INDsinsgam01=0
         INDsinsgam21=0
         end if

      if ((INDpalmer.eq.1).and.
     +((INDgam11.eq.1).or.(INDgam12.eq.1).or.(INDgam21.eq.1).or.
     +(INDgam22.eq.1).or.(INDgameq.eq.1))) then

         write(4,*) "CANNOT ESTIMATE BOTH AINSWORTH (2010) AND"
         write(4,*) "PALMER (2006) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING AINSWORTH (2010) gammas ONLY"
         warning=1
         write(4,*)
         INDpalmer=0
         end if

      if ((INDjc.eq.1).and.     
     +((INDgam11.eq.1).or.(INDgam12.eq.1).or.(INDgam21.eq.1).or.
     +(INDgam22.eq.1).or.(INDgameq.eq.1))) then

         write(4,*) "CANNOT ESTIMATE BOTH AINSWORTH (2010) AND"
         write(4,*) "LI (2009) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING AINSWORTH (2010) gammas ONLY"
         warning=1
         write(4,*)         
         INDjc=0
         end if

      if ((INDjc.eq.1).and.     
     +((INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1))) then

         write(4,*) "CANNOT ESTIMATE BOTH LI(2009) AND"
         write(4,*) "SINSHEIMER (2003) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING LI (2009) ONLY"
         warning=1
         write(4,*)
         INDsinsgam01=0
         INDsinsgam21=0
         end if

      if ((INDjc.eq.1).and.(INDpalmer.eq.1)) then

         write(4,*) "CANNOT ESTIMATE BOTH LI(2009) AND"
         write(4,*) "PALMER (2006) INTERACTIONS (gammas)"
         write(4,*) "ESTIMATING LI (2009) ONLY"
         warning=1
         write(4,*)
         INDpalmer=0
         end if


         
c     RE-USE TR(7) and TR(10) FOR Sinsheimer PARAMETERIZATION

      if (INDsinsgam01.eq.1) INDgam11=1
      if (INDsinsgam21.eq.1) INDgam22=1

      if (INDsinsgam01.eq.1) then
      c61='lnMFG01 '
      c62='sd_lnMFG01 '
      c63='cl_lnMFG01 '
      c64='cu_lnMFG01 '
      c71='lnMFG01 '
      end if

      if (INDsinsgam21.eq.1) then
      c91='lnMFG21 '
      c92='sd_lnMFG21 '
      c93='cl_lnMFG21 '
      c94='cu_lnMFG21 '
      end if

c     RE-USE TR(7) FOR PALMER PARAMETERIZATION

      if (INDpalmer.eq.1) INDgam11=1

      if (INDpalmer.eq.1) then
      c61='lnMFGmu '
      c62='sd_lnMFGmu '
      c63='cl_lnMFGmu '
      c64='cu_lnMFGmu '
      c71='lnMFGmu '
         end if


c     RE-USE TR(10) FOR LI PARAMETERIZATION

      if (INDjc.eq.1) INDgam22=1

      if (INDjc.eq.1) then
      c91='lnJc    '
      c92='sd_lnJc    '
      c93='cl_lnJc    '
      c94='cu_lnJc    '
      end if



      write(17,7777) char1, char2, charfr,
     +c01, c02, c03, c04,
     +c11, c12, c13, c14,
     +c21, c22, c23, c24,    
     +c31, c32, c33, c34,
     +c41, c42, c43, c44,
     +c51, c52, c53, c54, 
     +c61, c62, c63, c64,    
     +c71, c72, c73, c74,
     +c81, c82, c83, c84,
     +c91, c92, c93, c94, 
     +char3, char4, char5, char6

 7777 format(A4,A6,A5,6(A5,A8,A8,A8),4(A8,A11,A11,A11) 3A10, A9)





      if ((fixAIND.eq.1).and.(HWEIND.eq.0)) then
         write(4,*) "IF ALLELE FREQ FIXED, MUST ASSUME HWE"
         warning=1
         write(4,*)
         HWEIND=1
         end if

      if (((casetrioIND.ne.1).and.(casetrioIND.ne.2)).and.
     +((casemaduoIND.ne.1).and.(casemaduoIND.ne.2)).and.
     +((casepaduoIND.ne.1).and.(casepaduoIND.ne.2)).and.
     +(maofcaseIND.ne.1).and.
     +(paofcaseIND.ne.1).and.(caseIND.ne.1).and.
     +(conparIND.ne.1).and.
     +(conmaduoIND.ne.1).and.
     +(conpaduoIND.ne.1).and.
     +(conIND.ne.1)) then
      write(4,*) "NO INPUT DATA FILES SPECIFIED"
      write(4,*)
      goto 888
      end if

      if (((casetrioIND.ne.1).and.(casetrioIND.ne.2)).and.
     +((casemaduoIND.ne.1).and.(casemaduoIND.ne.2)).and.
     +((casepaduoIND.ne.1).and.(casepaduoIND.ne.2))) then
          write(4,*) "NO INPUT DATA FILES SPECIFIED WITH INFORMATION"
          write(4,*) "ON BOTH CASES AND PARENTS (MOTHER AND/OR FATHER)"
          write(4,*)  
          goto 888
         end if



      if ((HWEIND.eq.1).and.(ALSYMIND.eq.1)) then
      write(4,*) "Cannot assume both HWE and PARENTAL ALLELIC SYMMETRY"
      write(4,*) "ASSUMING HWE only"
         warning=1
      write(4,*)
      ALSYMIND=0
      write(4,*)
      end if

c      if ((fixAIND.eq.1).and.(ALSYMIND.eq.1)) then
c      write(4,*) "Cannot fix allele freq and PARENTAL ALLELIC SYMMETRY"
c      write(4,*) "Assuming PARENTAL ALLELIC SYMMETRY only"
c      fixAIND=0
c      write(4,*)
c      end if
c     THIS BIT NEVER GETS ACCESSED AS fixAIND implies HWEIND FROM ABOVE


      if ((CPGIND.eq.1).and.(ALSYMIND.eq.1)) then
      write(4,*) "Cannot assume both CPG and PARENTAL ALLELIC SYMMETRY"
      write(4,*) "ASSUMING  PARENTAL ALLELIC SYMMETRY only"
         warning=1
      write(4,*)
      CPGIND=0
      write(4,*)
      end if

      if ((HWEIND.eq.1).and.(CPGIND.eq.1)) then
      write(4,*) "Cannot assume both HWE and CPG"
      write(4,*) "ASSUMING HWE only"
         warning=1
      write(4,*)
      CPGIND=0
      write(4,*)
      end if




        if ((INDR1.eq.1).or.(INDR2.eq.1).or.(INDR1eqR2.eq.1)
     +.or.(INDR1sqeqR2.eq.1).or.(INDgam11.eq.1).or.(INDgam12.eq.1)
     +.or.(INDgam21.eq.1).or.(INDgam22.eq.1).or.(INDgameq.eq.1)
     +.or.(INDsinsgam01.eq.1).or.(INDsinsgam21.eq.1).or.
     +(INDpalmer.eq.1).or.(INDjc.eq.1)) then

      if ((INDIm.eq.1).and.(INDIp.eq.1)) then
      write(4,*) "Cannot estimate both Im and Ip (unless assume"
      write(4,*) "(no child genotype or interaction effects)"
      write(4,*) "Estimating Im only"
         warning=1
      write(4,*)
      INDIp=0
      end if

      end if


c	COMMENT OUT BELOW FOR DEC 2011 CHECK	
c	WHEN WERE THESE TWO CONDITIONS NEEDED ANYWAY?
c	PROBLEM IF TRYING TO ESTMATE BOTH R2 and gam12, and S2 and gam21?
c 	USE CHECK NOW INSERTED BELOW INSTEAD

c      if ((INDR2.eq.1).and.(INDgam12.eq.1)) then
c      write(4,*) "Cannot estimate both R2 and gamma12"
c      write(4,*) "Estimating R2 only"
c      write(4,*)
c      INDgam12=0
c      end if
c
c      if ((INDS2.eq.1).and.(INDgam21.eq.1)) then
c      write(4,*) "Cannot estimate both S2 and gamma21"
c      write(4,*) "Estimating S2 only"
c      write(4,*)
c      INDgam21=0
c      end if

c	USE THIS CHECK INSTEAD:

	if (((INDR2.eq.1).and.(INDgam12.eq.1)).and.
     +((INDS2.eq.1).and.(INDgam21.eq.1))) then
	write(4,*) "Cannot estimate all 4 of R2, gamma12, S2, gamma21"
	write(4,*) "Estimating R2, S2, gamma12 only"
         warning=1
	INDgam21=0
	end if






c     Extra restrictions if only have duos

c     NEED TO CHECK WHEN NEED TO FIX A 


c     Need to implement these for individual loci if required
c     - and change back afterwards!
c     Rather than implementing different restrictions for individual loci
c     safer to just report that cannot analyse that locus under desired model


      if (casetrioIND.eq.0) then
      

      if ((conparIND.eq.0).and.
     +((casemaduoIND.eq.0).or.(casepaduoIND.eq.0)).and.
     +(conmaduoIND.eq.0).and.(conpaduoIND.eq.0)) then
         HWEIND=1
         write(4,*) "No case/parent trios nor parents of controls"
         if (casemaduoIND.eq.0) write(4,*) "No case/mother duos"
         if (casepaduoIND.eq.0) write(4,*) "No case/father duos"
         write(4,*) "No information on mu's"
         write(4,*) "Must analyse assuming HWE"
         warning=1
         write(4,*)
         end if



      if ((conIND.eq.0).and.(conparIND.eq.0).and.
     +(conmaduoIND.eq.0).and.(conpaduoIND.eq.0)) then

      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "Estimating Im only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "Estimating Ip only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if

      end if

      

      if (casemaduoIND.eq.0) then

c     so only have case pa duos
c     Dont even try to estimate anything about mothers genotype
      
      write(4,*) "No case/parent trios or case/mother duos"
      write(4,*) "Cannot estimate maternal genotype effects (S1,S2)"
      write(4,*) "or mother/child interactions (gammas)"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam12=0
      INDgam21=0
      INDgam22=0
      INDgameq=0
      INDS1=0
      INDS2=0
      INDS1eqS2=0
      INDS1sqeqS2=0

      end if


      if (casepaduoIND.eq.0) then
c      so only have case ma duos


      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "Estimating Im only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "Estimating Ip only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if


      if ((conIND.eq.0).and.(conparIND.eq.0).and.
     +(conmaduoIND.eq.0).and.(conpaduoIND.eq.0)) then


      if ((INDIm.eq.1).and.(INDS1.eq.1).and.(INDS2.eq.1)
     +.and.(fixAIND.ne.1)) then
      write(4,*) "Cannot estimate all three of Im, S1, S2"
      write(4,*) "Estimating Im only"
         warning=1
      write(4,*)
      INDS1=0
      INDS2=0
      end if

      if ((INDIp.eq.1).and.(INDS1.eq.1).and.(INDS2.eq.1)
     +.and.(fixAIND.ne.1)) then
      write(4,*) "Cannot estimate all three of Ip, S1, S2"
      write(4,*) "Estimating Ip only"
         warning=1
      write(4,*)
      INDS1=0
      INDS2=0
      end if


         end if


      end if



c     else must have both case pa duos and  case ma duos
c     but no case/parent trios

      if ((conparIND.eq.0).and.(conmaduoIND.eq.0).and.
     +(conpaduoIND.eq.0).and.(HWEIND.eq.0)) then
      write(4,*) "Insufficient information on mu's"
      write(4,*) "Must analyse assuming HWE"   
            warning=1
      write(4,*)
      HWEIND=1
      end if


      if ((conparIND.eq.0).and.(HWEIND.eq.0)) then

      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "(unless assume HWE)"
      write(4,*) "Estimating Im only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "(unless assume HWE)"
      write(4,*) "Estimating Ip only"
         warning=1
      write(4,*)
      INDgam11=0
      INDgam22=0
      end if

      end if


      end if
c     ends if no case/parent trios



      if ((casetrioIND.eq.1).or.(casetrioIND.eq.2)) then
      open (unit=7, file=caseParentTriosFile, status='old')      
      read(7,*)
      end if

      if (parofcaseIND.eq.1) then
      open (unit=18, file=caseParentsFile, status='old')      
      read(18,*)
      end if

      if ((casemaduoIND.eq.1).or.(casemaduoIND.eq.2)) then
      open (unit=8, file=caseMotherDuosFile, status='old')      
      read(8,*)
      end if

      if ((casepaduoIND.eq.1).or.(casepaduoIND.eq.2)) then
      open (unit=9, file=casefatherduosFile, status='old')      
      read(9,*)
      end if

      if (maofcaseIND.eq.1) then
      open (unit=19, file=caseMothersFile, status='old')      
      read(19,*)
      end if

      if (paofcaseIND.eq.1) then
      open (unit=20, file=caseFathersFile, status='old')      
      read(20,*)
      end if

      if (caseIND.eq.1) then
      open (unit=10, file=casesFile, status='old')      
      read(10,*)
      end if

      if (conparIND.eq.1) then
      open (unit=13, file=conparentsFile, status='old')      
      read(13,*)
      end if

      if (conmaduoIND.eq.1) then
      open (unit=14, file=conmotherduosFile, status='old')      
      read(14,*)
      end if

      if (conpaduoIND.eq.1) then
      open (unit=15, file=confatherduosFile, status='old')      
      read(15,*)
      end if

      if (conIND.eq.1) then
      open (unit=16, file=consFile, status='old')      
      read(16,*)
      end if


      do 887 snp=1, nsnps


         warning=0

         lnlikfull=0.0
         lnliknull=0.0

         do 885 i=1,10
            THETA(i)=0.0
            sd(i)=0.0
 885        continue


         read(3, *, IOSTAT=IOStatus) snpno, fixA2
         if(IOStatus < 0) exit 

         fixA1=1.0-fixA2
         oldsnpno=snpno


         write(*,*) "ANALYSING SNP NUMBER",snp, ", SNPID=", snpno
         write(4,*)"--------------------------------------"
         write(4,*) "ANALYSING SNP NUMBER",snp, ", SNPID=", snpno
         write(4,*)


         if (fixAIND.eq.1) write(4,*) "ALLELE FREQ FIXED AT", fixA2


c       assume HWE to get starting values for mu's

         mew1dash=(1-fixA1)**4
         mew2dash=2*fixA1*((1-fixA1)**3)
         mew3dash=(fixA1**2)*((1-fixA1)**2)
         mew4dash=4*(fixA1**2)*((1-fixA1)**2)
         mew5dash=2*(fixA1**3)*(1-fixA1)
         mew6dash=fixA1**4

         start(12)=mew1dash
         start(13)=mew2dash*0.5
         start(14)=mew3dash
         start(15)=mew4dash*0.25
         start(16)=mew5dash*0.5
         start(17)=mew6dash


c     mew 1-6 are params 12-17
c     For CPG
c     let mew2a, 3a, 5a be params 13, 14, 16
c     and mew2b, 3b, 5b be params 18-20





      if (casetrioIND.eq.1) then

c     Read in data
c     Counts for genotypes of 
c         ma pa child

      read(7,*) snpno, (casetriocount(j), j=1,15)

      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if
      

      totcount=0
      do 20 j=1,15
      totcount=totcount+casetriocount(j)
 20   continue


c     Use these to get better starting estimates of mu's

      if (totcount.gt.0) then
 
      start(12)=casetriocount(1)/totcount
      start(13)=(casetriocount(2)+casetriocount(3)+
     +casetriocount(4)+casetriocount(5))/totcount
      start(14)=(casetriocount(6)+casetriocount(7))/totcount
      start(15)=(casetriocount(8)+casetriocount(9)+
     +casetriocount(10))/totcount
      start(16)=(casetriocount(11)+casetriocount(12)+
     +casetriocount(13)+casetriocount(14))/totcount
      start(17)=casetriocount(15)/totcount


c	if results in bad values, go back to assuming HWE

	if ((start(12).lt.0.0001).or.(start(13).lt.0.0001).or.
     +(start(14).lt.0.0001).or.(start(15).lt.0.0001).or.
     +(start(16).lt.0.0001).or.(start(17).lt.0.0001)) then
         start(12)=mew1dash
         start(13)=mew2dash*0.5
         start(14)=mew3dash
         start(15)=mew4dash*0.25
         start(16)=mew5dash*0.5
         start(17)=mew6dash
	end if
 
        end if

c      write(4,*) start(12), start(13), start(14), start(15), 
c     +start(16), start(17)

         end if


c     Start of code for estimated POO, read in extra cells
c     Cell 16=9a, 17=9b as in EMIM paper

      if (casetrioIND.eq.2) then

c     Read in data
c     Counts for genotypes of 
c         ma pa child

      read(7,*) snpno, (casetriocount(j), j=1,17), cntNoPhTRIO,
     +cntPhTRIO
      
      propNoPhTRIO=cntNoPhTRIO/(cntPhTRIO+cntNoPhTRIO)
      propPhTRIO=1-propNoPhTRIO

      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if
      

      totcount=0
      do 1020 j=1,17
      totcount=totcount+casetriocount(j)
 1020   continue


c     Use these to get better starting estimates of mu's

      if (totcount.gt.0) then
 
      start(12)=casetriocount(1)/totcount
      start(13)=(casetriocount(2)+casetriocount(3)+
     +casetriocount(4)+casetriocount(5))/totcount
      start(14)=(casetriocount(6)+casetriocount(7))/totcount
      start(15)=(casetriocount(8)+casetriocount(9)+
     +casetriocount(16)+casetriocount(17)+
     +casetriocount(10))/totcount
      start(16)=(casetriocount(11)+casetriocount(12)+
     +casetriocount(13)+casetriocount(14))/totcount
      start(17)=casetriocount(15)/totcount

c	if results in bad values, go back to assuming HWE

	if ((start(12).lt.0.0001).or.(start(13).lt.0.0001).or.
     +(start(14).lt.0.0001).or.(start(15).lt.0.0001).or.
     +(start(16).lt.0.0001).or.(start(17).lt.0.0001)) then
         start(12)=mew1dash
         start(13)=mew2dash*0.5
         start(14)=mew3dash
         start(15)=mew4dash*0.25
         start(16)=mew5dash*0.5
         start(17)=mew6dash
	end if
 
        end if

c      write(4,*) start(12), start(13), start(14), start(15), 
c     +start(16), start(17)

         end if

c     End of code for estimated POO, read in extra cells


         if (parofcaseIND.eq.1) then
            read(18,*) snpno, (parofcasecount(j), j=1,9)
            end if
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
        
      sumparofcase=0
      do 920 j=1,9
      sumparofcase=sumparofcase+parofcasecount(j)
 920   continue

            end if


         if (casemaduoIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        ma child
c            write(*,*)  "counting"
      read(8,*) snpno, (casemacount(j), j=1,7)


      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcasemaduo=0
      do 921 j=1,7
c           write(*,*)  "counting", casemacount(j)
      sumcasemaduo=sumcasemaduo+casemacount(j)
 921  continue

      end if

c     Start code for reading in phased mother/case duos
c
	          if (casemaduoIND.eq.2) then
c     Read in data
c     Counts for genotypes of 
c        ma child
c            write(*,*)  "counting"
      read(8,*) snpno, (casemacount(j), j=1,9), cntNoPhCAMD,
     +cntPhCAMD
      
      propNoPhCAMD=cntNoPhCAMD/(cntPhCAMD+cntNoPhCAMD)
      propPhCAMD=1-propNoPhCAMD     

      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcasemaduo=0
      do 1021 j=1,9
c           write(*,*)  "counting", casemacount(j)
      sumcasemaduo=sumcasemaduo+casemacount(j)
 1021  continue

      end if
c     End of code for reading in phased mother/case duos


         if (casepaduoIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        pa child
      read(9,*) snpno, (casepacount(j), j=1,7)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcasepaduo=0
      do 922 j=1,7
      sumcasepaduo=sumcasepaduo+casepacount(j)
 922  continue

      end if

c     Start code for reading in phased father/case duos
c
        if (casepaduoIND.eq.2) then
c     Read in data
c     Counts for genotypes of 
c        pa child
      read(9,*) snpno, (casepacount(j), j=1,9), cntNoPhCAFD, cntPhCAFD
      
	  propNoPhCAFD=cntNoPhCAFD/(cntPhCAFD+cntNoPhCAFD)
	  propPhCAFD=1-propNoPhCAFD

      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcasepaduo=0
      do 1022 j=1,9
      sumcasepaduo=sumcasepaduo+casepacount(j)
 1022 continue

      end if
c     End of code for reading in phased mother/case duos


         if (maofcaseIND.eq.1) then
            read(19,*) snpno, (maofcasecount(j), j=1,3)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      summaofcase=0
      do 923 j=1,3
      summaofcase=summaofcase+maofcasecount(j)
 923  continue

            end if

         if (paofcaseIND.eq.1) then
            read(20,*) snpno, (paofcasecount(j), j=1,3)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumpaofcase=0
      do 924 j=1,3
      sumpaofcase=sumpaofcase+paofcasecount(j)
 924  continue

            end if


         if (caseIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        child
      read(10,*) snpno, (casecount(j), j=1,3)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcase=0
      do 925 j=1,3
      sumcase=sumcase+casecount(j)
 925  continue


      end if

         if (conparIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        mother father
      read(13,*) snpno, (conparcount(j), j=1,9)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if
      end if


c     If have parents of controls, use them to produce better
c     starting estimates for mu's.

      partotcount=0
      do 21 j=1,9
      partotcount=partotcount+conparcount(j)
 21   continue


      if ((conparcount(1).gt.0).and.(conparcount(2).gt.0).and.
     +(conparcount(3).gt.0).and.(conparcount(4).gt.0).and.
     +(conparcount(5).gt.0).and.(conparcount(6).gt.0).and.
     +(conparcount(7).gt.0).and.(conparcount(8).gt.0).and.
     +(conparcount(9).gt.0)) then

      start(12)=conparcount(1)/partotcount
      start(13)=(conparcount(2)+conparcount(4))/(4*partotcount)
      start(14)=(conparcount(3)+conparcount(7))/(2*partotcount)
      start(15)=(conparcount(5))/(4*partotcount)
      start(16)=(conparcount(6)+conparcount(8))/(4*partotcount)
      start(17)=conparcount(9)/partotcount

      end if
 
c      write(4,*) start(12), start(13), start(14), start(15), 
c     +start(16), start(17)


c     Get starting values for (9 mu) CPG

c     mew 1-6 are params 12-17
c     For CPG
c     let mew2a, 3a, 5a be params 13, 14, 16
c     and mew2b, 3b, 5b be params 18-20

      cpgstart(12)=start(12)
      cpgstart(13)=start(13)/2.0
      cpgstart(14)=start(14)/2.0
      cpgstart(15)=start(15)
      cpgstart(16)=start(16)/2.0
      cpgstart(17)=start(17)
      cpgstart(18)=start(13)/2.0
      cpgstart(19)=start(14)/2.0
      cpgstart(20)=start(16)/2.0


         if (conmaduoIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        ma child
      read(14,*) snpno, (conmacount(j), j=1,7)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumconmaduo=0
      do 926 j=1,7
      sumconmaduo=sumconmaduo+conmacount(j)
 926  continue

      end if

         if (conpaduoIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        pa child
      read(15,*) snpno, (conpacount(j), j=1,7)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumconpaduo=0
      do 927 j=1,7
      sumconpaduo=sumconpaduo+conpacount(j)
 927  continue

      end if


         if (conIND.eq.1) then
c     Read in data
c     Counts for genotypes of 
c        child
      read(16,*) snpno, (concount(j), j=1,3)
      if (snpno.ne.oldsnpno) then
        write(4,*) "WARNING: SNP IDS IN INPUT FILES DON'T MATCH"  
         warning=1
c        write(4,*)
         end if

      sumcon=0
      do 928 j=1,3
      sumcon=sumcon+concount(j)
 928  continue

      end if



c     EXTRA CHECKS OF WHETHER NEED MORE RESTRICTIONS IF NOT ENOUGH DATA
C     IN SOME CATAGORIES

      
      sumcasetrio=totcount
      sumconpar=partotcount

      sumall=sumcasetrio+sumcasemaduo+sumcasepaduo+summaofcase
     +sumpaofcase+sumcase+sumconpar+sumconmaduo+sumconpaduo+sumcon


c      write(*,*) "sumall=",sumall
c      write(*,*) "sumcasetrio=",sumcasetrio
c      write(*,*) "sumcasemaduo=",sumcasemaduo
c      write(*,*) "sumcasepaduo=",sumcasepaduo
c      write(*,*) "summaofcase=",summaofcase
c      write(*,*) "sumpaofcase=",sumpaofcase
c      write(*,*) "sumcase=",sumcase
c      write(*,*) "sumconpar=",sumconpar
c      write(*,*) "sumconmaduo=",sumconmaduo
c      write(*,*) "sumconpaduo=",sumconpaduo
c      write(*,*) "sumcon=",sumcon



      if (sumall.eq.0) then
      write(4,*) "NO INPUT DATA PROVIDED FOR THIS LOCUS"
         warning=1
      write(4,*)
      goto 886
      end if

      if ((sumcasetrio.eq.0).and.(sumcasemaduo.eq.0)
     +.and.(sumcasepaduo.eq.0)) then
          write(4,*) "NO INPUT DATA FOR THIS LOCUS WITH INFORMATION"
          write(4,*) "ON BOTH CASES AND PARENTS (MOTHER AND/OR FATHER)"
         warning=1
          write(4,*)  
          goto 886
         end if


      if (sumcasetrio.eq.0) then
      
         
      if ((HWEIND.eq.0).and.(sumconpar.eq.0).and.
     +((sumcasemaduo.eq.0).or.(sumcasepaduo.eq.0)).and.
     +(sumconmaduo.eq.0).and.(sumconpaduo.eq.0)) then
         write(4,*) "No case/parent trios nor parents of controls"
         write(4,*) "for this locus"
         if (sumcasemaduo.eq.0) write(4,*) 
     +"No case/mother duos for this locus"
         if (sumcasepaduo.eq.0) write(4,*) 
     +"No case/father duos for this locus"
         write(4,*) "No information on mu's"
         write(4,*)
         write(4,*) "Can't implement desired analysis"
         write(4,*) "Try re-running assuming HWE"
         warning=1
         write(4,*)
         goto 886
         end if



      if ((sumcon.eq.0).and.(sumconpar.eq.0).and.
     +(sumconmaduo.eq.0).and.(sumconpaduo.eq.0)) then

      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Re-run analysis e.g. estimating Im only"
         warning=1
      write(4,*)         
         goto 886
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Re-run analysis e.g. estimating Ip only"
         warning=1
      write(4,*)
         goto 886
      end if

      end if




      if ((sumcasemaduo.eq.0).and.
     +((INDgam11.eq.1).or.
     +(INDgam12.eq.1).or.
     +(INDgam21.eq.1).or.
     +(INDgam22.eq.1).or.
     +(INDgameq.eq.1).or.
     +(INDS1.eq.1).or.
     +(INDS2.eq.1).or.
     +(INDS1eqS2.eq.1).or.
     +(INDS1sqeqS2.eq.1))) then

c     so only have case pa duos
c     Dont even try to estimate anything about mothers genotype
      
      write(4,*) "No case/parent trios or case/mother duos"
      write(4,*) "for this locus"
      write(4,*) "Cannot estimate maternal genotype effects (S1,S2)"
      write(4,*) "or mother/child interactions (gammas)"
      write(4,*)
      write(4,*) "Try re-runnning with alternative model"
         warning=1
      write(4,*)
        goto 886

      end if



      if (sumcasepaduo.eq.0) then
c      so only have case ma duos


      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Im only"
         warning=1
      write(4,*)
        goto 886
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Ip only"
         warning=1
      write(4,*)
        goto 886
      end if


      if ((sumcon.eq.0).and.(sumconpar.eq.0).and.
     +(sumconmaduo.eq.0).and.(sumconpaduo.eq.0)) then


      if ((INDIm.eq.1).and.(INDS1.eq.1).and.(INDS2.eq.1)
     +.and.(fixAIND.ne.1)) then
      write(4,*) "Cannot estimate all three of Im, S1, S2"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Im only"
         warning=1
      write(4,*)
        goto 886
      end if

      if ((INDIp.eq.1).and.(INDS1.eq.1).and.(INDS2.eq.1)
     +.and.(fixAIND.ne.1)) then
      write(4,*) "Cannot estimate all three of Ip, S1, S2"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Ip only"
         warning=1
      write(4,*)
        goto 886
      end if


         end if


      end if



c     else must have both case pa duos and  case ma duos
c     but no case/parent trios

      if ((sumconpar.eq.0).and.(sumconmaduo.eq.0).and.
     +(sumconpaduo.eq.0).and.(HWEIND.eq.0)) then
      write(4,*) "Insufficient information on mu's"
      write(4,*) "for this locus"
      write(4,*)
      write(4,*) "Try re-running assuming HWE"  
             warning=1
      write(4,*)
        goto 886
      end if


      if ((sumconpar.eq.0).and.(HWEIND.eq.0)) then

      if ((INDIm.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Im, gamma11, gamma22"
      write(4,*) "for this locus (unless assume HWE)"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Im only"
         warning=1
      write(4,*)
        goto 886
      end if

      if ((INDIp.eq.1).and.(INDgam11.eq.1).and.(INDgam22.eq.1)) then
      write(4,*) "Cannot estimate all three of Ip, gamma11, gamma22"
      write(4,*) "for this locus (unless assume HWE)"
      write(4,*)
      write(4,*) "Try re-runnning e.g. estimating Ip only"
         warning=1
      write(4,*)
        goto 886
      end if

      end if


      end if
c     ends if sumcasetrio.eq.0



        write(4,*)



c      write(*,*) start(12), start(13), start(14), start(15), 
c     +start(16), start(17)


C     Set up maxfun structures
C     See manual for details
c     DONT OUTPUT DETAILED OUTPUT FROM MAXFUN
C     AS FILE WILL GET TOO BIG

      IOUT     = 0
      IDET     = 0
      IXVC     = 2
      METHOD   =1

      EPSC1=1.0D-16
      IHIT=1

      NT=20


c      write(4,*) "INDgam=",INDgam11,INDgam12,INDgam21,INDgam22



c      NULL MODEL

c     Evaluate likelihood once as null

c     Dont estimate any params of interest

      do 352 i=1,10
      THIN(i) = 0.0
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 4  
 352  continue

         if (HWEIND.eq.0) then

           write(4,*) 'NOT ASSUMING HWE'

      THIN(11) = fixA1
      THL(11)  =1.0D-16
      THU(11)  = 0.999999
      ISTIN(11)= 4  


c     allele freq is fixed in the sense that mu's 
c     are estimated, not allele freq


      do 362 i=12,17
      THIN(i) = LOG(start(i))
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 1  
 362  continue


         if (ALSYMIND.eq.1) then
      ISTIN(15)= 3
           write(4,*) 'ASSUMING PARENTAL ALLELE SYMMETRY'
            end if

c     If assuming parental allelic symmetry need to make 
c     mewstar4=4mewstar3 i.e. mew4=mew3

      else

      write(4,*) 'ASSUMING HWE'    

         
      THIN(11) = fixA1
      THL(11)  =1.0D-16
      THU(11)  = 0.999999
      ISTIN(11)= 1  

c     by default, above says you should estimate allele freq

    
      if (fixAIND.eq.1) ISTIN(11)= 4 


      do 372 i=12,17
      THIN(i) = LOG(start(i))
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 4  
 372  continue

         end if


         if (CPGIND.eq.0) then

      do 472 i=18,20
      THIN(i) = 0.0
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 4  
 472  continue

         end if


         if (CPGIND.eq.1) then

      do 572 i=12,20
      THIN(i) = cpgstart(i)
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 1 
 572  continue

         end if






      MAXIT= 5000

c     If nothing estimated, need to just evaluate likelihood once
c     Need to pretend we are estimating something (e.g. allele freq)
c     for Maxfun to work?

      if ((fixAIND.eq.1).and.(HWEIND.eq.1)) then
      MAXIT=0
      ISTIN(11)= 1
      end if


c     Initialise AV(i,i) to zero for re-running MAXFUN repeatedly
c     (In case maximisation fails, you are left with old values)

      do 800 i=1,20
         AV(i,i)=0.0
 800     continue

      CALL MAXFUN(FUNCTION, DEPAR1, THETA, F, NFE, LFL)
      write(4,*)
      write(4,*) "RUNNING NULL MODEL"
c      write(*,*) "RUNNING NULL MODEL"

        lnliknull=F

	R1=exp(THETA(1))
	R2=exp(THETA(2))
	S1=exp(THETA(3))
	S2=exp(THETA(4))
	Im=exp(THETA(5))
	Ip=exp(THETA(6))
	gam11=exp(THETA(7))
	gam12=exp(THETA(8))
	gam21=exp(THETA(9))
	gam22=exp(THETA(10))
        A1=THETA(11)

        if (HWEIND.eq.0) then

        mew1=exp(THETA(12))
        mew2=exp(THETA(13))
        mew3=exp(THETA(14))
        mew4=exp(THETA(15))
        mew5=exp(THETA(16))
        mew6=exp(THETA(17))


        if (CPGIND.eq.1) then
        mew2a=exp(THETA(13))
        mew3a=exp(THETA(14))
        mew5a=exp(THETA(16))
        mew2b=exp(THETA(18))
        mew3b=exp(THETA(19))
        mew5b=exp(THETA(20))
               end if


        else

c       if assuming HWE

c           write(4,*) "THETA11=", THETA(11)


         mew1dash=(1-THETA(11))**4
         mew2dash=2*THETA(11)*((1-THETA(11))**3)
         mew3dash=(THETA(11)**2)*((1-THETA(11))**2)
         mew4dash=4*(THETA(11)**2)*((1-THETA(11))**2)
         mew5dash=2*(THETA(11)**3)*(1-THETA(11))
         mew6dash=THETA(11)**4

         mew1=mew1dash
         mew2=mew2dash*0.5
         mew3=mew3dash
         mew4=mew4dash*0.25
         mew5=mew5dash*0.5
         mew6=mew6dash

         end if


        write(4,*)
        write(4,*) "NULL MODEL RESULTS:"

        write(4,*) "PARAMETER ESTIMATES"
        write(4,*)

c        if (HWEIND.eq.0) then
c        write(4,*) "Fixed (not used) allele freq A2 =", 1.0-THETA(11)
c        end if

        if ((HWEIND.eq.1).and.(fixAIND.ne.1)) then
        write(4,*) "Estimated allele freq A2 =", 1.0-THETA(11)
        end if

        if ((HWEIND.eq.1).and.(fixAIND.eq.1)) then
        write(4,*) "Fixed allele freq A2 =", 1.0-THETA(11)
        end if

        if (CPGIND.ne.1) then
        write(4,*) "mu1=", mew1
        write(4,*) "mu2=", mew2
        write(4,*) "mu3=", mew3
        write(4,*) "mu4=", mew4
        write(4,*) "mu5=", mew5
        write(4,*) "mu6=", mew6
        end if

       if (CPGIND.eq.1) then
        write(4,*) "mu1=", mew1
        write(4,*) "mu2a=", mew2a
        write(4,*) "mu2b=", mew2b
        write(4,*) "mu3a=", mew3a
        write(4,*) "mu3b=", mew3b
        write(4,*) "mu4=", mew4
        write(4,*) "mu5a=", mew5a
        write(4,*) "mu5b=", mew5b
        write(4,*) "mu6=", mew6
        end if



        write(4,*)

        write(4,*) "R1=", R1
        write(4,*) "R2=", R2
        write(4,*) "S1=", S1
        write(4,*) "S2=", S2
        write(4,*) "Im=", Im
        write(4,*) "Ip=", Ip
        write(4,*) "gamma11=", gam11
        write(4,*) "gamma12=", gam12
        write(4,*) "gamma21=", gam21
        write(4,*) "gamma22=", gam22
        write(4,*)


c     ALTERNATIVE MODEL


c      write(4,*) "INDgam=",INDgam11,INDgam12,INDgam21,INDgam22


     
c     First assume we are estimating all specified
c     params of interest

      do 350 i=1,10

      THIN(i) = 0.0
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 1  

 350     continue

         
      if (INDR1.eq.0) ISTIN(1)=4
      if (INDR2.eq.0) ISTIN(2)=4
      if (INDS1.eq.0) ISTIN(3)=4
      if (INDS2.eq.0) ISTIN(4)=4
      if (INDIm.eq.0) ISTIN(5)=4
      if (INDIp.eq.0) ISTIN(6)=4
      if (INDgam11.eq.0) ISTIN(7)=4
      if (INDgam12.eq.0) ISTIN(8)=4
      if (INDgam21.eq.0) ISTIN(9)=4
      if (INDgam22.eq.0) ISTIN(10)=4


      if (INDR1eqR2.eq.1)   ISTIN(2)=3
      if (INDR1sqeqR2.eq.1) ISTIN(2)=3
      if (INDS1eqS2.eq.1)   ISTIN(4)=3
      if (INDS1sqeqS2.eq.1) ISTIN(4)=3


      if (INDgameq.eq.1) then
         ISTIN(8)=3
         ISTIN(9)=3
         ISTIN(10)=3
         end if

c     Make sure you are actually estimating something in
c     case of restrictions as specified - if necessary
c     overruling what they said

      if (INDR1eqR2.eq.1)   ISTIN(1)=1
      if (INDR1sqeqR2.eq.1) ISTIN(1)=1
      if (INDS1eqS2.eq.1)   ISTIN(3)=1
      if (INDS1sqeqS2.eq.1) ISTIN(3)=1
      if (INDgameq.eq.1) ISTIN(7)=1


c      write(4,*) "ISTIN(8)=",ISTIN(8)
c      write(4,*) "INDgam=",INDgam11,INDgam12,INDgam21,INDgam22


         if (HWEIND.eq.0) then

c           write(4,*) 'NOT ASSUMING HWE'

      THIN(11) = fixA1
      THL(11)  =1.0D-16
      THU(11)  = 0.999999
      ISTIN(11)= 4  


c     allele freq is fixed in the sense that mu's 
c     are estimated, not allele freq



      do 360 i=12,17
      THIN(i) = LOG(start(i))
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 1  
 360     continue


         if (ALSYMIND.eq.1) then
      ISTIN(15)= 3
c           write(4,*) 'ASSUMING PARENTAL ALLELE SYMMETRY'
            end if

c     If assuming parental allelic symmetry need to make 
c     mewstar4=4mewstar3 i.e. mew4=mew3

      else

c      write(4,*) 'ASSUMING HWE'    

      THIN(11) = fixA1
      THL(11)  =1.0D-16
      THU(11)  = 0.999999
      ISTIN(11)= 1  

c     by default, above says you should estimate allele freq

    
      if (fixAIND.eq.1) ISTIN(11)= 4 


      do 370 i=12,17
      THIN(i) = LOG(start(i))
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 4  
 370     continue

         end if


         if (CPGIND.eq.0) then

      do 470 i=18,20
      THIN(i) = 0.0
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 4  
 470  continue

         end if

         if (CPGIND.eq.1) then

      do 570 i=12,20
      THIN(i) = cpgstart(i)
      THL(i)  = -1000.0
      THU(i)  = 1000.0
      ISTIN(i)= 1 
 570  continue

         end if


      MAXIT= 5000

c     If nothing estimated, need to just evaluate likelihood once
c     Need to pretend we are estimating something (e.g. allele freq)
c     for Maxfun to work?

c      write(4,*) 
c      write(4,*) 'TESTING, fixAIND=',fixAIND
c      write(4,*) 'TESTING, HWEIND=',HWEIND
c      write(4,*) 

c      if (((INDR1+INDR2+INDS1+INDS2+INDIm+INDIp+
c     +INDgam11+INDgam12+INDgam21+INDgam22).eq.0).and.(fixAIND.eq.1)
c     +.and.(HWEIND.eq.1)) then

      if (((ISTIN(1)+ISTIN(2)+ISTIN(3)+ISTIN(4)+ISTIN(5)+ISTIN(6)+
     +ISTIN(7)+ISTIN(8)+ISTIN(9)+ISTIN(10)).eq.0).and.(fixAIND.eq.1)
     +.and.(HWEIND.eq.1)) then

      MAXIT=0
      ISTIN(11)= 1

      end if


c      THL(12)  = -30.0
c      THU(12)  = 30.0

c     Initialise AV(i,i) to zero for re-running MAXFUN repeatedly
c     (In case maximisation fails, you are left with old values)

      do 801 i=1,20
         AV(i,i)=0.0
 801  continue

      CALL MAXFUN(FUNCTION, DEPAR1, THETA, F, NFE, LFL)
      write(4,*)
      write(4,*) "RUNNING ALTERNATIVE MODEL"
c      write(*,*) "RUNNING ALTERNATIVE MODEL"

        if ((ISTIN(1).eq.1).or.(ISTIN(1).eq.3))   IND(1)=1
        if ((ISTIN(2).eq.1).or.(ISTIN(2).eq.3))   IND(2)=1
        if ((ISTIN(3).eq.1).or.(ISTIN(3).eq.3))   IND(3)=1
        if ((ISTIN(4).eq.1).or.(ISTIN(4).eq.3))   IND(4)=1
        if ((ISTIN(5).eq.1).or.(ISTIN(5).eq.3))   IND(5)=1
        if ((ISTIN(6).eq.1).or.(ISTIN(6).eq.3))   IND(6)=1
        if ((ISTIN(7).eq.1).or.(ISTIN(7).eq.3))   IND(7)=1
        if ((ISTIN(8).eq.1).or.(ISTIN(8).eq.3))   IND(8)=1
        if ((ISTIN(9).eq.1).or.(ISTIN(9).eq.3))   IND(9)=1
        if ((ISTIN(10).eq.1).or.(ISTIN(10).eq.3))  IND(10)=1


c     Jeremie found that you actually need to output those with ISTIN=3 
c     after all (MAXFUN does output them)
 

         numestparam=0
         
      totvar=0.00

        do 910 i=1,10
c           write(4,*) "ISTIN",i,"=",ISTIN(i)
c           write(4,*) "IND",i,"=",IND(i)
c           write(4,*) "AV",i,"=",AV(i,i)
              if (IND(i).eq.1) then
                 totvar=totvar+AV(i,i)
                 end if
                  if (ISTIN(i).eq.1) numestparam=numestparam+1
910       continue


c      If CIs not estimable, re-run with mew1 and mew6 
c      bounds altered to hopefully get MAXFUN to fix at bound

c       write(4,*) "ESTIMATED TOTAL VARIANCE=", totvar

      if (totvar.le.1.0D-16) then


      do 371 i=1,20
      THL(i)  = -25.0
      THU(i)  = 25.0
 371  continue


  
c      THL(12)  = -25.0
c      THU(12)  = 25.0   
c      THL(17)  = -25.0
c      THU(17)  = 25.0   


c     Initialise AV(i,i) to zero for re-running MAXFUN repeatedly
c     (In case maximisation fails, you are left with old values)

      do 802 i=1,20
         AV(i,i)=0.0
 802  continue


      CALL MAXFUN(FUNCTION, DEPAR1, THETA, F, NFE, LFL)
      write(4,*) "RE-RUNNING ALTERNATIVE MODEL"

      end if


c      write(4,*)
c        do 911 i=1,10
c           write(4,*) "IND",i,"=",IND(i)
c           write(4,*) "AV",i,"=",AV(i,i)
c 911    continue

        if ((ISTIN(1).eq.1).or.(ISTIN(1).eq.3))   IND(1)=1
        if ((ISTIN(2).eq.1).or.(ISTIN(2).eq.3))   IND(2)=1
        if ((ISTIN(3).eq.1).or.(ISTIN(3).eq.3))   IND(3)=1
        if ((ISTIN(4).eq.1).or.(ISTIN(4).eq.3))   IND(4)=1
        if ((ISTIN(5).eq.1).or.(ISTIN(5).eq.3))   IND(5)=1
        if ((ISTIN(6).eq.1).or.(ISTIN(6).eq.3))   IND(6)=1
        if ((ISTIN(7).eq.1).or.(ISTIN(7).eq.3))   IND(7)=1
        if ((ISTIN(8).eq.1).or.(ISTIN(8).eq.3))   IND(8)=1
        if ((ISTIN(9).eq.1).or.(ISTIN(9).eq.3))   IND(9)=1
        if ((ISTIN(10).eq.1).or.(ISTIN(10).eq.3))  IND(10)=1


c     Jeremie found that you actually need to output those with ISTIN=3 
c     after all (MAXFUN does output them)

	R1=exp(THETA(1))
	R2=exp(THETA(2))
	S1=exp(THETA(3))
	S2=exp(THETA(4))
	Im=exp(THETA(5))
	Ip=exp(THETA(6))
	gam11=exp(THETA(7))
	gam12=exp(THETA(8))
	gam21=exp(THETA(9))
	gam22=exp(THETA(10))
        A1=THETA(11)

        if (HWEIND.eq.0) then

        mew1=exp(THETA(12))
        mew2=exp(THETA(13))
        mew3=exp(THETA(14))
        mew4=exp(THETA(15))
        mew5=exp(THETA(16))
        mew6=exp(THETA(17))


        if (CPGIND.eq.1) then
        mew2a=exp(THETA(13))
        mew3a=exp(THETA(14))
        mew5a=exp(THETA(16))
        mew2b=exp(THETA(18))
        mew3b=exp(THETA(19))
        mew5b=exp(THETA(20))
               end if

        else

c       if assuming HWE

c           write(4,*) "THETA11=", THETA(11)


         mew1dash=(1-THETA(11))**4
         mew2dash=2*THETA(11)*((1-THETA(11))**3)
         mew3dash=(THETA(11)**2)*((1-THETA(11))**2)
         mew4dash=4*(THETA(11)**2)*((1-THETA(11))**2)
         mew5dash=2*(THETA(11)**3)*(1-THETA(11))
         mew6dash=THETA(11)**4

         mew1=mew1dash
         mew2=mew2dash*0.5
         mew3=mew3dash
         mew4=mew4dash*0.25
         mew5=mew5dash*0.5
         mew6=mew6dash

         end if


        write(4,*)
        write(4,*) "ALTERNATIVE MODEL RESULTS:"

        write(4,*) "PARAMETER ESTIMATES"
        write(4,*)

c        if (HWEIND.eq.0) then
c        write(4,*) "Fixed (not used) allele freq A2 =", 1.0-THETA(11)
c        end if


        if ((HWEIND.eq.1).and.(fixAIND.ne.1)) then
        write(4,*) "Estimated allele freq A2 =", 1.0-THETA(11)
        end if

        if ((HWEIND.eq.1).and.(fixAIND.eq.1)) then
        write(4,*) "Fixed allele freq A2 =", 1.0-THETA(11)
        end if



        if (CPGIND.ne.1) then
        write(4,*) "mu1=", mew1
        write(4,*) "mu2=", mew2
        write(4,*) "mu3=", mew3
        write(4,*) "mu4=", mew4
        write(4,*) "mu5=", mew5
        write(4,*) "mu6=", mew6
        end if

       if (CPGIND.eq.1) then
        write(4,*) "mu1=", mew1
        write(4,*) "mu2a=", mew2a
        write(4,*) "mu2b=", mew2b
        write(4,*) "mu3a=", mew3a
        write(4,*) "mu3b=", mew3b
        write(4,*) "mu4=", mew4
        write(4,*) "mu5a=", mew5a
        write(4,*) "mu5b=", mew5b
        write(4,*) "mu6=", mew6
        end if

        write(4,*)







        write(4,*) "R1=", R1
        write(4,*) "R2=", R2
        write(4,*) "S1=", S1
        write(4,*) "S2=", S2
        if (INDweinIm.ne.1) write(4,*) "Im=", Im
        if (INDweinIm.eq.1) write(4,*) "Im=", Im,
     +"[Weinberg (1999b) parameterization]"
        if (INDweinIp.ne.1) write(4,*) "Ip=", Ip
        if (INDweinIp.eq.1) write(4,*) "Ip=", Ip, 
     +"[Weinberg (1999b) and Li (2009) parameterization]"
        if ((INDsinsgam01.ne.1).and.(INDpalmer.ne.1)) write(4,*) 
     +"gamma11=", gam11
        if (INDsinsgam01.eq.1) write(4,*) "gamma01=", gam11, 
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDpalmer.eq.1) write(4,*) "mu=", gam11,
     +"[Palmer (2006) match parameter]"
        write(4,*) "gamma12=", gam12
        write(4,*) "gamma21=", gam21
        if ((INDsinsgam21.ne.1).and.(INDjc.ne.1)) write(4,*) 
     +"gamma22=", gam22
        if (INDsinsgam21.eq.1) write(4,*) "gamma21=", gam22,
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDjc.eq.1) write(4,*) "Jc=", gam22,
     +"[Li (2009) and Parami (2008) conflict parameter]"
        write(4,*)

      
        j=0
        do 110 i=1,10

           sd(i)=0

c           write(4,*) "AV",i,"=",AV(i,i)

           if (IND(i).eq.1) then
              j=j+1
              sd(i)=sqrt(AV(j,j))
              end if
 110          continue
           
        
        write(4,*) "ln R1=", THETA(1), "SE=",sd(1)
        write(4,*) "ln R2=", THETA(2), "SE=",sd(2)
        write(4,*) "ln S1=", THETA(3), "SE=",sd(3)
        write(4,*) "ln S2=", THETA(4), "SE=",sd(4)
        if (INDweinIm.ne.1) write(4,*) "ln Im=", THETA(5), "SE=",sd(5)
        if (INDweinIm.eq.1) write(4,*) "ln Im=", THETA(5), "SE=",sd(5),
     +"[Weinberg (1999b) parameterization]"
        if (INDweinIp.ne.1) write(4,*) "ln Ip=", THETA(6), "SE=",sd(6)
        if (INDweinIp.eq.1) write(4,*) "ln Ip=", THETA(6), "SE=",sd(6),
     +"[Weinberg (1999b) and Li (2009) parameterization]"
        if ((INDsinsgam01.ne.1).and.(INDpalmer.ne.1))        
     +write(4,*) "ln gamma11=", THETA(7), "SE=",sd(7)
        if (INDsinsgam01.eq.1) write(4,*) "ln gamma01=", THETA(7), 
     +"SE=",sd(7), 
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDpalmer.eq.1) write(4,*) "ln mu=", THETA(7), "SE=",sd(7),
     +"[Palmer (2006) match parameter]"
        write(4,*) "ln gamma12=", THETA(8), "SE=",sd(8)
        write(4,*) "ln gamma21=", THETA(9), "SE=",sd(9)
        if ((INDsinsgam21.ne.1).and.(INDjc.ne.1))
     +write(4,*) "ln gamma22=", THETA(10), "SE=",sd(10)
       if (INDsinsgam21.eq.1) write(4,*) "ln gamma21=", THETA(10), 
     +"SE=",sd(10), 
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDjc.eq.1) write(4,*) "ln Jc=", THETA(10), 
     +"SE=",sd(10), 
     +"[Li (2009) and Parami (2008) conflict parameter]"

        write(4,*)


        do 810 i=1,10           
           if ((IND(i).eq.1).and.(sd(i).le.1.0D-16)) then
        write(4,*) "WARNING: ESTIMATED PARAMETERS WITH SE=",sd(i)  
        write(4,*) "SUGGESTS INSUFFICIENT DATA AVAILABLE" 
        write(4,*) "TO ESTIMATE THESE PARAMETERS"
        write(4,*) "(OTHER PARAMETER ESTIMATES SHOULD BE OK)"
         warning=1
              write(4,*)
              goto 812
              end if
 810          continue
 812          continue

              if ((sd(1).le.1.0D-16).and.(sd(2).le.1.0D-16).and.
     +(sd(3).le.1.0D-16).and.(sd(4).le.1.0D-16).and.
     +(sd(5).le.1.0D-16).and.(sd(6).le.1.0D-16).and.
     +(sd(7).le.1.0D-16).and.(sd(8).le.1.0D-16).and.
     +(sd(9).le.1.0D-16).and.(sd(10).le.1.0D-16)) then
        write(4,*) "IF NO PARAMETERS HAVE SE>0, TRY RE-RUNNING EMIM"
        write(4,*) "WITH REDUCED NUMBER OF PARAMETERS"
         warning=1
        write(4,*) 
        end if


        write(4,*)
        write(4,*) "95% CI for parameters of interest"
        write(4,*)

        do 811 i=1,10           
           if ((IND(i).eq.1).and.(sd(i).le.1.0D-16)) then
        write(4,*) "WARNING: PARAMETERS ABOVE WITH SE=",sd(i) 
        write(4,*) "SUGGESTS INSUFFICIENT DATA AVAILABLE" 
        write(4,*) "TO ESTIMATE THESE PARAMETERS"
         warning=1
              write(4,*)
              goto 813
              end if
 811       continue
 813       continue



        if (IND(1).eq.1) write(4,*) "R1=", 
     +exp(THETA(1)-1.96*sd(1)), 
     +exp(THETA(1)+1.96*sd(1))
        if (IND(2).eq.1) write(4,*) "R2=", 
     +exp(THETA(2)-1.96*sd(2)), 
     +exp(THETA(2)+1.96*sd(2))
        if (IND(3).eq.1) write(4,*) "S1=", 
     +exp(THETA(3)-1.96*sd(3)), 
     +exp(THETA(3)+1.96*sd(3))
        if (IND(4).eq.1) write(4,*) "S2=", 
     +exp(THETA(4)-1.96*sd(4)), 
     +exp(THETA(4)+1.96*sd(4))

        if ((IND(5).eq.1).and.(INDweinIm.ne.1)) 
     +write(4,*) "Im=", 
     +exp(THETA(5)-1.96*sd(5)), 
     +exp(THETA(5)+1.96*sd(5))
        if ((IND(5).eq.1).and.(INDweinIm.eq.1)) 
     +write(4,*) "Im=", 
     +exp(THETA(5)-1.96*sd(5)), 
     +exp(THETA(5)+1.96*sd(5)),
     +"[Weinberg (1999b) parameterization]"

        if ((IND(6).eq.1).and.(INDweinIp.ne.1)) 
     +write(4,*) "Ip=", 
     +exp(THETA(6)-1.96*sd(6)), 
     +exp(THETA(6)+1.96*sd(6))
        if ((IND(6).eq.1).and.(INDweinIp.eq.1)) 
     +write(4,*) "Ip=", 
     +exp(THETA(6)-1.96*sd(6)), 
     +exp(THETA(6)+1.96*sd(6)),
     +"[Weinberg (1999b) and Li (2009) parameterization]"

       if (IND(7).eq.1) then
        if ((INDsinsgam01.ne.1).and.(INDpalmer.ne.1))  
     +write(4,*) "gamma11=", 
     +exp(THETA(7)-1.96*sd(7)), 
     +exp(THETA(7)+1.96*sd(7))
        if (INDsinsgam01.eq.1)  write(4,*) 
     +"gamma01=", 
     +exp(THETA(7)-1.96*sd(7)), 
     +exp(THETA(7)+1.96*sd(7)),
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDpalmer.eq.1)  write(4,*) 
     +"mu=", 
     +exp(THETA(7)-1.96*sd(7)), 
     +exp(THETA(7)+1.96*sd(7)),
     +"[Palmer (2006) match parameter]"
        end if
     
       if (IND(8).eq.1) write(4,*) "gamma12=", 
     +exp(THETA(8)-1.96*sd(8)), 
     +exp(THETA(8)+1.96*sd(8))
       if (IND(9).eq.1) write(4,*) "gamma21=", 
     +exp(THETA(9)-1.96*sd(9)), 
     +exp(THETA(9)+1.96*sd(9))

       if (IND(10).eq.1) then
        if ((INDsinsgam21.ne.1).and.(INDjc.ne.1))  
     +write(4,*) "gamma22=", 
     +exp(THETA(10)-1.96*sd(10)), 
     +exp(THETA(10)+1.96*sd(10))
       if (INDsinsgam21.eq.1)  write(4,*) 
     +"gamma21=", 
     +exp(THETA(10)-1.96*sd(10)), 
     +exp(THETA(10)+1.96*sd(10)),
     +"[Sinsheimer (2003) MFG parameterization]"
       if (INDjc.eq.1)  write(4,*) 
     +"Jc=", 
     +exp(THETA(10)-1.96*sd(10)), 
     +exp(THETA(10)+1.96*sd(10)),
     +"[Li (2009) and Parami (2008) conflict parameter]"
          end if

        write(4,*)  
        write(4,*)  


         if (IND(1).eq.1) write(4,*) "ln R1=", 
     +THETA(1)-1.96*sd(1), 
     +THETA(1)+1.96*sd(1)
         if (IND(2).eq.1) write(4,*) "ln R2=", 
     +THETA(2)-1.96*sd(2), 
     +THETA(2)+1.96*sd(2)
         if (IND(3).eq.1) write(4,*) "ln S1=", 
     +THETA(3)-1.96*sd(3), 
     +THETA(3)+1.96*sd(3)
         if (IND(4).eq.1) write(4,*) "ln S2=", 
     +THETA(4)-1.96*sd(4), 
     +THETA(4)+1.96*sd(4)

        if ((IND(5).eq.1).and.(INDweinIm.ne.1)) 
     +write(4,*) "ln Im=", 
     +THETA(5)-1.96*sd(5), 
     +THETA(5)+1.96*sd(5)
        if ((IND(5).eq.1).and.(INDweinIm.eq.1)) 
     +write(4,*) "ln Im=", 
     +THETA(5)-1.96*sd(5), 
     +THETA(5)+1.96*sd(5),
     +"[Weinberg (1999b) parameterization]"

      if ((IND(6).eq.1).and.(INDweinIp.ne.1)) 
     +write(4,*) "ln Ip=", 
     +THETA(6)-1.96*sd(6), 
     +THETA(6)+1.96*sd(6)
        if ((IND(6).eq.1).and.(INDweinIp.eq.1)) 
     +write(4,*) "ln Ip=",
     +THETA(6)-1.96*sd(6), 
     +THETA(6)+1.96*sd(6),
     +"[Weinberg (1999b) and Li (2009) parameterization]"

        if (IND(7).eq.1) then
        if ((INDsinsgam01.ne.1).and.(INDpalmer.ne.1))
     +write(4,*) "ln gamma11=", 
     +THETA(7)-1.96*sd(7), 
     +THETA(7)+1.96*sd(7)
        if (INDsinsgam01.eq.1)  write(4,*) 
     +"ln gamma01=", 
     +THETA(7)-1.96*sd(7), 
     +THETA(7)+1.96*sd(7),
     +"[Sinsheimer (2003) MFG parameterization]"
        if (INDpalmer.eq.1)  write(4,*) 
     +"ln mu=", 
     +THETA(7)-1.96*sd(7), 
     +THETA(7)+1.96*sd(7),
     +"[Palmer (2006) match parameter]"
        end if


        if (IND(8).eq.1) write(4,*) "ln gamma12=", 
     +THETA(8)-1.96*sd(8), 
     +THETA(8)+1.96*sd(8)
        if (IND(9).eq.1) write(4,*) "ln gamma21=", 
     +THETA(9)-1.96*sd(9), 
     +THETA(9)+1.96*sd(9)

        if (IND(10).eq.1) then
        if ((INDsinsgam21.ne.1).and.(INDjc.ne.1)) 
     +write(4,*) "ln gamma22=", 
     +THETA(10)-1.96*sd(10), 
     +THETA(10)+1.96*sd(10)
       if (INDsinsgam21.eq.1)  write(4,*) 
     +"ln gamma21=", 
     +THETA(10)-1.96*sd(10), 
     +THETA(10)+1.96*sd(10),
     +"[Sinsheimer (2003) MFG parameterization]"
       if (INDjc.eq.1)  write(4,*) 
     +"ln Jc=", 
     +THETA(10)-1.96*sd(10), 
     +THETA(10)+1.96*sd(10),
     +"[Li (2009) and Parami (2008) conflict parameter]"
          end if

c        write(4,*) "SDs are", (sqrt(AV(i,i)), i=1,10)




        lnlikfull=F

        write(4,*)  
        write(4,*) "COMPARISON OF MODEL RESULTS:"
        write(4,*)  
        write(4,*) "Null ln likelihood=", lnliknull
c        write(4,*)
        write(4,*) "Maximized ln likelihood=", lnlikfull 
c        write(4,*)
        write(4,*) "Twice difference in ln likelihoods=", 
     +2*(lnlikfull-lnliknull)
        write(4,*)
        write(4,*) "Compare to chi-squared on",numestparam, "df"
        write(4,*)
        write(4,*) "Estimated parameters contributing to these df:"
c        write(4,*)

          if (ISTIN(1).eq.1) write(4,*) "R1"
          if (ISTIN(2).eq.1) write(4,*) "R2"
          if (ISTIN(3).eq.1) write(4,*) "S1"
          if (ISTIN(4).eq.1) write(4,*) "S2"
          if (ISTIN(5).eq.1) then
         if (INDweinIm.ne.1) write(4,*) "Im"
         if (INDweinIm.eq.1) write(4,*) 
     +"Im [Weinberg 1999b parameterization]"
          end if
          if (ISTIN(6).eq.1) then
          if (INDweinIp.ne.1) write(4,*) "Ip"
          if (INDweinIp.eq.1) write(4,*) 
     +"Ip [Weinberg 1999b; Li 2009 parameterization]"
          end if
          if (ISTIN(7).eq.1) then
          if ((INDsinsgam01.ne.1).and.(INDpalmer.ne.1)) 
     +write(4,*) "gamma11"
          if (INDsinsgam01.eq.1)  write(4,*)  
     +"gamma01 [Sinsheimer 2003 parameterization]"
          if (INDpalmer.eq.1)  write(4,*) 
     +"mu [Palmer 2006 parameterization]"
             end if
          if (ISTIN(8).eq.1) write(4,*) "gamma12"
          if (ISTIN(9).eq.1) write(4,*) "gamma21"
          if (ISTIN(10).eq.1) then        
          if ((INDsinsgam21.ne.1).and.(INDjc.ne.1))    
     +write(4,*) "gamma22"
       if (INDsinsgam21.eq.1)  
     +write(4,*) "gamma21 [Sinsheimer 2003 parameterization]"
      if (INDjc.eq.1)  write(4,*) 
     +"Jc =exp(i_c) [Li 2009; Parami 2008 parameterization]"
      end if

          write(4,*)


c      write(4,*)
c      write(4,*) "Var/cov matrix of ln RRs" 
c      write(4,*)
c
c      do 200 k=1,17
c      write(4,990) (AV(k,l), l=1,17)
c 200  continue

 990  format(17f10.4)

 886  continue

      write(17,890) snp, snpno, 1.0-THETA(11), 
     +(THETA(i), sd(i), 
     +THETA(i)-1.96*sd(i), THETA(i)+1.96*sd(i), i=1,10),
     +lnliknull,lnlikfull, 2*(lnlikfull-lnliknull), warning
 890  format(i12,f20.3, f10.5, 40f15.5, 3f15.5, i5)



 887  continue
c     ENDS DO LOOP FOR DIFFERENT SNPS


 888  continue



c     close input/output files

      close(3)
      close(4)
      close(5)
      close(7)
      close(8)
      close(9)
      close(10)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)


      stop
      END PROGRAM EMIM




	SUBROUTINE FUNCTION(TR, FTR, NFE, LEX)
      use shared_var
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C The function to be evaluated:
C  Parameters:  TR  - trial parameters
C               FTR - function value with parameters tr
C               NFE - number of times the function has been evaluated
C               LEX - error code
c


c     LOCAL VARIABLES

      
      double precision R1, R2, S1, S2, Im, Ip, gam11, gam12, 
     +gam21, gam22, 
     +A1, mew1, mew2, mew3, mew4, mew5, mew6, tot, celltot(17), 
     +parcell(9), duocell(9), cell(17), totcount,
     +multinomcell(17), casecell(3), concell(3),
     +mew1dash, mew2dash, mew3dash, mew4dash, mew5dash, mew6dash,
     +mew2a, mew2b, mew3a, mew3b, mew5a, mew5b

      DIMENSION TR(*)
      EXTERNAL DEPAR1

C 	Check the parameters are valid
      	CALL DEPAR1(TR, LEX)
      	IF (LEX .GT. 0) RETURN

C Calculate function

c       FUNCTION MAXIMISED IS SUPPOSED TO BE LN LIK, IF ONE WISHES 
c       TO INTERPRET
c	INVERSE OF HESSIAN AS VAR/COV MARTRIX



	FTR=0


	R1=exp(TR(1))
	R2=exp(TR(2))
	S1=exp(TR(3))
	S2=exp(TR(4))
	Im=exp(TR(5))
	Ip=exp(TR(6))
	gam11=exp(TR(7))
	gam12=exp(TR(8))
	gam21=exp(TR(9))
	gam22=exp(TR(10))
	A1=TR(11)



        if (HWEIND.eq.0) then

        mew1=exp(TR(12))
        mew2=exp(TR(13))
        mew3=exp(TR(14))
        mew4=exp(TR(15))
        mew5=exp(TR(16))
        mew6=exp(TR(17))

        if (CPGIND.eq.1) then
        mew2a=exp(TR(13))
        mew3a=exp(TR(14))
        mew5a=exp(TR(16))
        mew2b=exp(TR(18))
        mew3b=exp(TR(19))
        mew5b=exp(TR(20))
               end if

        else

c       if assuming HWE

         mew1dash=(1-TR(11))**4
         mew2dash=2*TR(11)*((1-TR(11))**3)
         mew3dash=(TR(11)**2)*((1-TR(11))**2)
         mew4dash=4*(TR(11)**2)*((1-TR(11))**2)
         mew5dash=2*(TR(11)**3)*(1-TR(11))
         mew6dash=TR(11)**4

         mew1=mew1dash
         mew2=mew2dash*0.5
         mew3=mew3dash
         mew4=mew4dash*0.25
         mew5=mew5dash*0.5
         mew6=mew6dash

         end if


c        write(4,*) 'params= ', (exp(TR(i)), i=1,17)

c     Likelihood is basically a multinomial, but with 15 cell probs
c     parameterized in terms of parameters as above

c     Also add in terms for case/mothers, case/fathers, cases,
c     controls/parents, controls/mothers, control/fathers, controls
        


         if (CPGIND.eq.0) then
            mew2a=mew2
            mew2b=mew2
            mew3a=mew3
            mew3b=mew3
            mew5a=mew5
            mew5b=mew5
            end if


        celltot(1)=R2*S2*Im*Ip*gam22*mew1

        celltot(2)=R2*S2*Im*Ip*gam22*mew2a
        celltot(3)=R1*S2*Im*gam21*mew2a
        celltot(4)=R2*S1*Im*Ip*gam12*mew2b
        celltot(5)=R1*S1*Ip*gam11*mew2b

        celltot(6)=R1*S2*Im*gam21*mew3a
        celltot(7)=R1*Ip*mew3b

        celltot(8)=R2*S1*Im*Ip*gam12*mew4
        celltot(9)=R1*S1*(Im+Ip)*gam11*mew4
        celltot(10)=S1*mew4
	
c       cells will be added or not in appro props later
		celltot(16)=R1*S1*Ip*gam11*mew4
		celltot(17)=R1*S1*Im*gam11*mew4
         
        celltot(11)=R1*S1*Im*gam11*mew5a
        celltot(12)=S1*mew5a
        celltot(13)=R1*Ip*mew5b
        celltot(14)=mew5b

        celltot(15)=mew6

 
        if ((INDweinIm.eq.1).or.((INDweinIp.eq.1).or.
     +(INDsinsgam01.eq.1).or.INDsinsgam21.eq.1).or.
     +(INDpalmer.eq.1).or.(INDjc.eq.1)) then

c     THIS IN PRINCIPAL ALLOWS BOTH WEINBERG Im AND WEINBERG Ip
c     TO BE ESTIMATED (PROVIDED Rs AND gammas not estimated)

        celltot(1)=mew1*R2*S2
        if (INDpalmer.eq.1) celltot(1)=mew1*R2*S2*gam11
        if ((INDIm.eq.1).and.(INDweinIm.ne.1)) celltot(1)=celltot(1)*Im
        if ((INDIp.eq.1).and.(INDweinIp.ne.1)) celltot(1)=celltot(1)*Ip

c        write(4,*) "INDIp, INDweinIp", INDIp, INDweinIp
c        write(4,*) "mew1, R2, S2, Ip=",mew1, R2, S2,  Ip
c        write(4,*) "celltot(1)=",celltot(1)

        celltot(2)=mew2a*R2*S2
        if (INDpalmer.eq.1) celltot(2)=mew2a*R2*S2*gam11
        if ((INDIm.eq.1).and.(INDweinIm.ne.1)) celltot(2)=celltot(2)*Im
        if ((INDIp.eq.1).and.(INDweinIp.ne.1)) celltot(2)=celltot(2)*Ip

        celltot(3)=mew2a*R1*S2*Im
        if (INDsinsgam21.eq.1) celltot(3)=mew2a*R1*S2*Im*gam22
        if (INDjc.eq.1) celltot(3)=mew2a*R1*S2*Im*gam22

        celltot(4)=mew2b*R2*S1
        if (INDpalmer.eq.1) celltot(4)=mew2b*R2*S1*gam11
        if (INDjc.eq.1) celltot(4)=mew2b*R2*S1*gam22
        if ((INDIm.eq.1).and.(INDweinIm.ne.1)) celltot(4)=celltot(4)*Im
        if ((INDIp.eq.1).and.(INDweinIp.ne.1)) celltot(4)=celltot(4)*Ip

        celltot(5)=mew2b*R1*S1*Ip
        if (INDpalmer.eq.1) celltot(5)=mew2b*R1*S1*Ip*gam11

        celltot(6)=mew3a*R1*S2*Im
        if (INDsinsgam21.eq.1) celltot(6)=mew3a*R1*S2*Im*gam22
        if (INDjc.eq.1) celltot(6)=mew3a*R1*S2*Im*gam22

        celltot(7)=mew3b*R1*Ip
        if (INDsinsgam01.eq.1) celltot(7)=mew3b*R1*Ip*gam11
        if (INDjc.eq.1) celltot(7)=mew3b*R1*Ip*gam22

        celltot(8)=mew4*R2*S1
        if (INDpalmer.eq.1) celltot(8)=mew4*R2*S1*gam11
        if (INDjc.eq.1) celltot(8)=mew4*R2*S1*gam22
        if ((INDIm.eq.1).and.(INDweinIm.ne.1)) celltot(8)=celltot(8)*Im
        if ((INDIp.eq.1).and.(INDweinIp.ne.1)) celltot(8)=celltot(8)*Ip

        celltot(9)=mew4*R1*S1*(Im+Ip)
        if (INDpalmer.eq.1) celltot(9)=mew4*R1*S1*(Im+Ip)*gam11

		 celltot(16)=mew4*R1*S1*Ip
         if (INDpalmer.eq.1) celltot(16)=mew4*R1*S1*Ip*gam11

		 celltot(17)=mew4*R1*S1*Im
         if (INDpalmer.eq.1) celltot(17)=mew4*R1*S1*Im*gam11
		 
        celltot(10)=mew4*S1
        if (INDpalmer.eq.1) celltot(10)=mew4*S1*gam11
        if (INDjc.eq.1) celltot(10)=mew4*S1*gam22

        celltot(11)=mew5a*R1*S1*Im
        if (INDpalmer.eq.1) celltot(11)=mew5a*R1*S1*Im*gam11
        celltot(12)=mew5a*S1
        if (INDpalmer.eq.1) celltot(12)=mew5a*S1*gam11
        if (INDjc.eq.1) celltot(12)=mew5a*S1*gam22

        celltot(13)=mew5b*R1*Ip
        if (INDsinsgam01.eq.1) celltot(13)=mew5b*R1*Ip*gam11
        if (INDjc.eq.1) celltot(13)=mew5b*R1*Ip*gam22
        celltot(14)=mew5b
        if (INDpalmer.eq.1) celltot(14)=mew5b*gam11

        celltot(15)=mew6
        if (INDpalmer.eq.1) celltot(15)=mew6*gam11

        end if
       


        tot=celltot(1)+celltot(2)+celltot(3)+
     +celltot(4)+celltot(5)+celltot(6)+
     +celltot(7)+celltot(8)+celltot(9)+
     +celltot(10)+celltot(11)+celltot(12)+
     +celltot(13)+celltot(14)+celltot(15)

c     already added thro cell 9
c	    if (casetrioIND.eq.2) then
c		 tot=tot+celltot(16)+celltot(17)
c		 end if


        if ((casetrioIND.eq.1).or.(casetrioIND.eq.2)) then



c        totcount=casetriocount(1)+casetriocount(2)+
c     +casetriocount(3)+
c     +casetriocount(4)+casetriocount(5)+casetriocount(6)+
c     +casetriocount(7)+casetriocount(8)+casetriocount(9)+
c     +casetriocount(10)+casetriocount(11)+casetriocount(12)+
c     +casetriocount(13)+casetriocount(14)+casetriocount(15)
c
c        do 500 i=1,15
c          cell(i)=celltot(i)/tot
c          multinomcell(i)=casetriocount(i)/totcount
c 500      continue
c
c          write(4,995) (cell(i), i=1,15)
c          write(4,995) (multinomcell(i), i=1,15)
c          write(4,*)
c
c 995      format(15f8.5)
c





	FTR=FTR+
     +casetriocount(1)*LOG(celltot(1)/tot)+
     +casetriocount(2)*LOG(celltot(2)/tot)+
     +casetriocount(3)*LOG(celltot(3)/tot)+
     +casetriocount(4)*LOG(celltot(4)/tot)+
     +casetriocount(5)*LOG(celltot(5)/tot)+
     +casetriocount(6)*LOG(celltot(6)/tot)+
     +casetriocount(7)*LOG(celltot(7)/tot)+
     +casetriocount(8)*LOG(celltot(8)/tot)+    
     +casetriocount(10)*LOG(celltot(10)/tot)+
     +casetriocount(11)*LOG(celltot(11)/tot)+
     +casetriocount(12)*LOG(celltot(12)/tot)+
     +casetriocount(13)*LOG(celltot(13)/tot)+
     +casetriocount(14)*LOG(celltot(14)/tot)+
     +casetriocount(15)*LOG(celltot(15)/tot)

		 
      if (casetrioIND.eq.1) then
	   FTR=FTR+
     +casetriocount(9)*LOG(celltot(9)/tot)      
	  else if (casetrioIND.eq.2) then

	   if (propNoPhTRIO>0) then
	    FTR=FTR+casetriocount(9)*LOG(propNoPhTRIO*celltot(9)/tot)
	   end if

	   if (propPhTRIO>0) then 
        FTR=FTR+casetriocount(16)*LOG(propPhTRIO*celltot(16)/tot)+
     +casetriocount(17)*LOG(propPhTRIO*celltot(17)/tot)
	   end if
 
	  end if

        end if


        if (casemaduoIND.eq.1) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(2)
          duocell(2)=celltot(3)+celltot(6)
          duocell(3)=celltot(4)+celltot(8)
          duocell(4)=celltot(5)+celltot(9)+celltot(11)
          duocell(5)=celltot(10)+celltot(12)
          duocell(6)=celltot(7)+celltot(13)
          duocell(7)=celltot(14)+celltot(15)

              FTR=FTR+casemacount(1)*LOG(duocell(1)/tot)+
     +casemacount(2)*LOG(duocell(2)/tot)+
     +casemacount(3)*LOG(duocell(3)/tot)+
     +casemacount(4)*LOG(duocell(4)/tot)+
     +casemacount(5)*LOG(duocell(5)/tot)+
     +casemacount(6)*LOG(duocell(6)/tot)+
     +casemacount(7)*LOG(duocell(7)/tot)
          
           end if


      if (casemaduoIND.eq.2) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(2)
          duocell(2)=celltot(3)+celltot(6)
          duocell(3)=celltot(4)+celltot(8)
          duocell(4)=propNoPhCAMD*(celltot(5)+celltot(9)+celltot(11))
c          duocell(4)=celltot(5)+celltot(9)+celltot(11)
          duocell(5)=celltot(10)+celltot(12)
          duocell(6)=celltot(7)+celltot(13)
          duocell(7)=celltot(14)+celltot(15)
          duocell(8)=propPhCAMD*(celltot(16) + celltot(5))
          duocell(9)=propPhCAMD*(celltot(17) + celltot(11))

           FTR=FTR+casemacount(1)*LOG(duocell(1)/tot)+
     +casemacount(2)*LOG(duocell(2)/tot)+
     +casemacount(3)*LOG(duocell(3)/tot)+
     +casemacount(5)*LOG(duocell(5)/tot)+
     +casemacount(6)*LOG(duocell(6)/tot)+
     +casemacount(7)*LOG(duocell(7)/tot)  

c          FTR=FTR+casemacount(8)*LOG(duocell(4)/tot)
c          FTR=FTR+casemacount(9)*LOG(duocell(4)/tot)

c		  write(*, *) 'counts'
c		  write(*, *) casemacount(1), casemacount(2)
c		  write(*, *) casemacount(3), casemacount(4)
c		  write(*, *) casemacount(5), casemacount(6)
c		  write(*, *) casemacount(7), casemacount(8)
c		  write(*, *) casemacount(9)

c      write(*, *) duocell(4)

       if (propNoPhCAMD>0) then 
        FTR=FTR+casemacount(4)*LOG(duocell(4)/tot)    
      end if

      if (propPhCAMD>0) then 
        FTR=FTR+casemacount(8)*LOG(duocell(8)/tot)+
     +casemacount(9)*LOG(duocell(9)/tot)
       end if

c	   if(FTR*0.ne.0) then
c	     write(*, *) duocell(1), duocell(2), duocell(3), duocell(4), duocell(5)
c		 write(*, *) duocell(6), duocell(7), duocell(8), duocell(9)
c	     write(*, *) FTR
c	   end if

c	     write(*, *) duocell(1), duocell(2), duocell(3), duocell(4), duocell(5)
c		 write(*, *) duocell(6), duocell(7), duocell(8), duocell(9)

      end if



        if (casepaduoIND.eq.1) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(4)
          duocell(2)=celltot(5)+celltot(7)
          duocell(3)=celltot(2)+celltot(8)
          duocell(4)=celltot(3)+celltot(9)+celltot(13)
          duocell(5)=celltot(10)+celltot(14)
          duocell(6)=celltot(6)+celltot(11)
          duocell(7)=celltot(12)+celltot(15)

              FTR=FTR+casepacount(1)*LOG(duocell(1)/tot)+
     +casepacount(2)*LOG(duocell(2)/tot)+
     +casepacount(3)*LOG(duocell(3)/tot)+
     +casepacount(4)*LOG(duocell(4)/tot)+
     +casepacount(5)*LOG(duocell(5)/tot)+
     +casepacount(6)*LOG(duocell(6)/tot)+
     +casepacount(7)*LOG(duocell(7)/tot)
          
           end if


      if (casepaduoIND.eq.2) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(4)
          duocell(2)=celltot(5)+celltot(7)
          duocell(3)=celltot(2)+celltot(8)
          duocell(4)=propNoPhCAFD*(celltot(3)+celltot(9)+celltot(13))
          duocell(5)=celltot(10)+celltot(14)
          duocell(6)=celltot(6)+celltot(11)
          duocell(7)=celltot(12)+celltot(15)
		  duocell(8)=propPhCAFD*(celltot(16)+celltot(13))
		  duocell(9)=propPhCAFD*(celltot(17)+celltot(3))

              FTR=FTR+casepacount(1)*LOG(duocell(1)/tot)+
     +casepacount(2)*LOG(duocell(2)/tot)+
     +casepacount(3)*LOG(duocell(3)/tot)+
     +casepacount(5)*LOG(duocell(5)/tot)+
     +casepacount(6)*LOG(duocell(6)/tot)+
     +casepacount(7)*LOG(duocell(7)/tot)
      
	   if (propNoPhCAFD>0) then 
         FTR=FTR+casepacount(4)*LOG(duocell(4)/tot)
	   end if  
	     
	   if (propPhCAFD>0) then 
        FTR=FTR+casepacount(8)*LOG(duocell(8)/tot)+
     +casepacount(9)*LOG(duocell(9)/tot)
	   end if

      end if


        if (caseIND.eq.1) then

          casecell(1)=celltot(1)+celltot(2)+celltot(4)+celltot(8)
          casecell(2)=celltot(3)+celltot(5)+celltot(6)+celltot(7)+
     +celltot(9)+celltot(11)+celltot(13)
          casecell(3)=celltot(10)+celltot(12)+celltot(14)+celltot(15)

          FTR=FTR+casecount(1)*LOG(casecell(1)/tot)+
     +casecount(2)*LOG(casecell(2)/tot)+
     +casecount(3)*LOG(casecell(3)/tot)

           end if


        if (maofcaseIND.eq.1) then

          casecell(1)=celltot(1)+celltot(2)+celltot(3)+celltot(6)
          casecell(2)=celltot(4)+celltot(8)+celltot(5)+celltot(9)+
     +celltot(11)+celltot(10)+celltot(12)
          casecell(3)=celltot(7)+celltot(13)+celltot(14)+celltot(15)

          FTR=FTR+maofcasecount(1)*LOG(casecell(1)/tot)+
     +maofcasecount(2)*LOG(casecell(2)/tot)+
     +maofcasecount(3)*LOG(casecell(3)/tot)

           end if


        if (paofcaseIND.eq.1) then

          casecell(1)=celltot(1)+celltot(4)+celltot(5)+celltot(7)
          casecell(2)=celltot(2)+celltot(8)+celltot(3)+celltot(9)+
     +celltot(13)+celltot(10)+celltot(14)
          casecell(3)=celltot(6)+celltot(11)+celltot(12)+celltot(15)

          FTR=FTR+paofcasecount(1)*LOG(casecell(1)/tot)+
     +paofcasecount(2)*LOG(casecell(2)/tot)+
     +paofcasecount(3)*LOG(casecell(3)/tot)

           end if



        if (parofcaseIND.eq.1) then

          parcell(1)=celltot(1)
          parcell(2)=celltot(2)+celltot(3)
          parcell(3)=celltot(6)
          parcell(4)=celltot(4)+celltot(5)
          parcell(5)=celltot(8)+celltot(9)+celltot(10)
          parcell(6)=celltot(11)+celltot(12)
          parcell(7)=celltot(7)
          parcell(8)=celltot(13)+celltot(14)
          parcell(9)=celltot(15)

              FTR=FTR+parofcasecount(1)*LOG(parcell(1)/tot)+
     +parofcasecount(2)*LOG(parcell(2)/tot)+
     +parofcasecount(3)*LOG(parcell(3)/tot)+
     +parofcasecount(4)*LOG(parcell(4)/tot)+
     +parofcasecount(5)*LOG(parcell(5)/tot)+
     +parofcasecount(6)*LOG(parcell(6)/tot)+
     +parofcasecount(7)*LOG(parcell(7)/tot)+
     +parofcasecount(8)*LOG(parcell(8)/tot)+
     +parofcasecount(9)*LOG(parcell(9)/tot)         

           end if


c     RESET CELL PROBS FOR CONTROLS



        celltot(1)=mew1

        celltot(2)=mew2a
        celltot(3)=mew2a
        celltot(4)=mew2b
        celltot(5)=mew2b

        celltot(6)=mew3a
        celltot(7)=mew3b

        celltot(8)=mew4
        celltot(9)=2*mew4		
        celltot(10)=mew4

        celltot(11)=mew5a
        celltot(12)=mew5a
        celltot(13)=mew5b
        celltot(14)=mew5b

        celltot(15)=mew6

        tot=celltot(1)+celltot(2)+celltot(3)+
     +celltot(4)+celltot(5)+celltot(6)+
     +celltot(7)+celltot(8)+celltot(9)+
     +celltot(10)+celltot(11)+celltot(12)+
     +celltot(13)+celltot(14)+celltot(15)
	  
       if (conparIND.eq.1) then
                          
c              do as proper multinomial in terms of parcells

          parcell(1)=celltot(1)
          parcell(2)=celltot(2)+celltot(3)
          parcell(3)=celltot(6)
          parcell(4)=celltot(4)+celltot(5)
          parcell(5)=celltot(8)+celltot(9)+celltot(10)
          parcell(6)=celltot(11)+celltot(12)
          parcell(7)=celltot(7)
          parcell(8)=celltot(13)+celltot(14)
          parcell(9)=celltot(15)


              FTR=FTR+conparcount(1)*LOG(parcell(1)/tot)+
     +conparcount(2)*LOG(parcell(2)/tot)+
     +conparcount(3)*LOG(parcell(3)/tot)+
     +conparcount(4)*LOG(parcell(4)/tot)+
     +conparcount(5)*LOG(parcell(5)/tot)+
     +conparcount(6)*LOG(parcell(6)/tot)+
     +conparcount(7)*LOG(parcell(7)/tot)+
     +conparcount(8)*LOG(parcell(8)/tot)+
     +conparcount(9)*LOG(parcell(9)/tot)         

           end if


        if (conmaduoIND.eq.1) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(2)
          duocell(2)=celltot(3)+celltot(6)
          duocell(3)=celltot(4)+celltot(8)
          duocell(4)=celltot(5)+celltot(9)+celltot(11)
          duocell(5)=celltot(10)+celltot(12)
          duocell(6)=celltot(7)+celltot(13)
          duocell(7)=celltot(14)+celltot(15)

              FTR=FTR+conmacount(1)*LOG(duocell(1)/tot)+
     +conmacount(2)*LOG(duocell(2)/tot)+
     +conmacount(3)*LOG(duocell(3)/tot)+
     +conmacount(4)*LOG(duocell(4)/tot)+
     +conmacount(5)*LOG(duocell(5)/tot)+
     +conmacount(6)*LOG(duocell(6)/tot)+
     +conmacount(7)*LOG(duocell(7)/tot)
          
           end if


        if (conpaduoIND.eq.1) then
                         
c              do as proper multinomial in terms of duocells

          duocell(1)=celltot(1)+celltot(4)
          duocell(2)=celltot(5)+celltot(7)
          duocell(3)=celltot(2)+celltot(8)
          duocell(4)=celltot(3)+celltot(9)+celltot(13)
          duocell(5)=celltot(10)+celltot(14)
          duocell(6)=celltot(6)+celltot(11)
          duocell(7)=celltot(12)+celltot(15)

              FTR=FTR+conpacount(1)*LOG(duocell(1)/tot)+
     +conpacount(2)*LOG(duocell(2)/tot)+
     +conpacount(3)*LOG(duocell(3)/tot)+
     +conpacount(4)*LOG(duocell(4)/tot)+
     +conpacount(5)*LOG(duocell(5)/tot)+
     +conpacount(6)*LOG(duocell(6)/tot)+
     +conpacount(7)*LOG(duocell(7)/tot)
          
           end if



        if (conIND.eq.1) then

          concell(1)=celltot(1)+celltot(2)+celltot(4)+celltot(8)
          concell(2)=celltot(3)+celltot(5)+celltot(6)+celltot(7)+
     +celltot(9)+celltot(11)+celltot(13)
          concell(3)=celltot(10)+celltot(12)+celltot(14)+celltot(15)

          FTR=FTR+concount(1)*LOG(concell(1)/tot)+
     +concount(2)*LOG(concell(2)/tot)+
     +concount(3)*LOG(concell(3)/tot)

           end if

c           write(4,*) "FTR=", FTR



C Add 1 to the times we have evaluated
      NFE = NFE + 1
c       print*,'nfe=',nfe
C Set no error
      LEX = 0

C And return
      RETURN
      END






	SUBROUTINE HESSIN(DF,HESS,OBS,PAR,F,ITER,KASE,MODEL,NCASE
     :,NOBS,NPAR,DIFFER)
C
C     THIS SUBROUTINE PERMITS RECOMPUTATION OF THE OBJECTIVE FUNCTION
C     F AND ITS DIFFERENTIAL DF.  IF DIFFER(2) IS TRUE, THEN PROVIDE
C     AN APPROXIMATION TO THE SECOND DIFFERENTIAL, I.E, HESSIAN OF F.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION DF(NPAR),HESS(NPAR,NPAR),OBS(NOBS)
     :,PAR(NPAR)
      LOGICAL DIFFER(2)
C
C
      END








	SUBROUTINE DEPAR1(TR,LEX)
      use shared_var
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C The function to check the parameters are in bounds and to calculate
C dependent parameters.
C
      PARAMETER (NP=45)

	COMMON /MAXF1/ THIN(NP),THL(NP),THU(NP),STPIN(NP),EPSD,YOTA,EPST,
     $               EPSC1,EPSC2,EPSC3,ISTIN(NP),NT,MAXIT,METHOD,IXVC,
     $               IHIT
      COMMON /MKUND/ KP2,KK
      DIMENSION TR(*)

c       Local variables

      double precision cell(17), R1, R2, S1, S2, 
     +Im, Ip, gam11, gam12, gam21, gam22, 
     +A1, mew1, mew2, mew3, mew4, mew5, mew6, tot, celltot(17),
     +mew1dash, mew2dash, mew3dash, mew4dash, mew5dash, mew6dash,
     +mew2a, mew2b, mew3a, mew3b, mew5a, mew5b
C
C Calculate the dependent parameters
c

        if (HWEIND.eq.0) then
        if (ALSYMIND.eq.1) then
           TR(15)=TR(14)
           end if
           end if

           if (INDR1eqR2.eq.1) then
          TR(2)=TR(1)    
              end if

           if (INDR1sqeqR2.eq.1) then
          TR(2)=2*TR(1)    
              end if

           if (INDS1eqS2.eq.1) then
          TR(4)=TR(3)    
              end if

           if (INDS1sqeqS2.eq.1) then
          TR(4)=2*TR(3)    
              end if

              if (INDgameq.eq.1) then
                 TR(8)=TR(7)
                 TR(9)=TR(7)
                 TR(10)=TR(7)
                 end if


C     Check that all parameters are in bounds 


	do 20 i=1,NT

	   if ((TR(i).lt.-100.D0).or.(TR(i).gt.1000.0D0)) then
         LEX = 1
         RETURN
       END IF

 20	continue


c     Adjust lower bound for ln mews and (indeed ln other params)
c     to -100 since Vals data suggested -500 could cause problems


        if ((TR(11).lt.0.D0).or.(TR(11).gt.1.0D0)) then
         LEX = 1
         RETURN
       END IF


	R1=exp(TR(1))
	R2=exp(TR(2))
	S1=exp(TR(3))
	S2=exp(TR(4))
	Im=exp(TR(5))
	Ip=exp(TR(6))
	gam11=exp(TR(7))
	gam12=exp(TR(8))
	gam21=exp(TR(9))
	gam22=exp(TR(10))
	A1=TR(11)


        if (HWEIND.eq.0) then

        mew1=exp(TR(12))
        mew2=exp(TR(13))
        mew3=exp(TR(14))
        mew4=exp(TR(15))
        mew5=exp(TR(16))
        mew6=exp(TR(17))

        if (CPGIND.eq.1) then
        mew2a=exp(TR(13))
        mew3a=exp(TR(14))
        mew5a=exp(TR(16))
        mew2b=exp(TR(18))
        mew3b=exp(TR(19))
        mew5b=exp(TR(20))
               end if

        else

c       if assuming HWE

         mew1dash=(1-TR(11))**4
         mew2dash=2*TR(11)*((1-TR(11))**3)
         mew3dash=(TR(11)**2)*((1-TR(11))**2)
         mew4dash=4*(TR(11)**2)*((1-TR(11))**2)
         mew5dash=2*(TR(11)**3)*(1-TR(11))
         mew6dash=TR(11)**4

         mew1=mew1dash
         mew2=mew2dash*0.5
         mew3=mew3dash
         mew4=mew4dash*0.25
         mew5=mew5dash*0.5
         mew6=mew6dash

         end if




         if (CPGIND.eq.0) then
            mew2a=mew2
            mew2b=mew2
            mew3a=mew3
            mew3b=mew3
            mew5a=mew5
            mew5b=mew5
            end if

        celltot(1)=R2*S2*Im*Ip*gam22*mew1

        celltot(2)=R2*S2*Im*Ip*gam22*mew2a
        celltot(3)=R1*S2*Im*gam21*mew2a
        celltot(4)=R2*S1*Im*Ip*gam12*mew2b
        celltot(5)=R1*S1*Ip*gam11*mew2b

        celltot(6)=R1*S2*Im*gam21*mew3a
        celltot(7)=R1*Ip*mew3b

        celltot(8)=R2*S1*Im*Ip*gam12*mew4
        celltot(9)=R1*S1*(Im+Ip)*gam11*mew4
        celltot(10)=S1*mew4
	
		celltot(16)=R1*S1*Ip*gam11*mew4
		celltot(17)=R1*S1*Im*gam11*mew4

        celltot(11)=R1*S1*Im*gam11*mew5a
        celltot(12)=S1*mew5a
        celltot(13)=R1*Ip*mew5b
        celltot(14)=mew5b

        celltot(15)=mew6

        tot=celltot(1)+celltot(2)+celltot(3)+
     +celltot(4)+celltot(5)+celltot(6)+
     +celltot(7)+celltot(8)+celltot(9)+
     +celltot(10)+celltot(11)+celltot(12)+
     +celltot(13)+celltot(14)+celltot(15)


c     Model is overparameterized in sense that get 
c     same likelihood if multiply all mu's by a const

c     Could fix so that mu's chosen to make celltots
c     all equal to probs that add up to 1, however gives
c     problems

c     Better just leave it - seems to work OK (!)

      
    
             
			 do 40 i=1,17

         cell(i)=celltot(i)/tot

c       write(4,*) i, cell(i)

           if ((cell(i).lt.0.D0).or.(cell(i).gt.1.0D0)
     +.and.((casetrioIND.eq.1).or.(i.ne.9))) then
         LEX = 1
         RETURN
       END IF

 40         continue
      
	
      LEX = 0

      RETURN
      END
c
