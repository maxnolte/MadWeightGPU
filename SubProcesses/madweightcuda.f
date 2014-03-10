c----------MadWeight Matrix Element Function in CUDA Fortran ----------
c--                 -------------- Edited by: Max Nolte ------------ 
c --       only for p p > H > mu+ mu- generated with heftold-full

      module madweightcuda
      use dhelascuda

c     CONSTANT VALUES
      integer, constant :: max_particles, max_branches, max_configs, max_channel
      integer, constant :: nexternal, num_inv, nincoming, npoints_d, ndim_d
      integer, constant :: inv_matching_type_part_d(3:12), NHEL(4,16)
      logical, constant :: mirrorprocs_d(2), lpp_d(2)
      double precision, constant :: pi, s_d, px_visible_d, py_visible_d, q2fact_d(2)
      double precision, constant :: mm_d,mh_d,wh_d,mc_d,t_ISR(2)
      double complex, constant :: gc_8_d,GC_42_d,GC_40_d
      double precision, constant ::  missPhi_EXP_d, missPT_EXP_d
      double precision, constant :: c_point_d(1:12,3,2), pexp_d(0:3,4)
      double precision, constant :: tf_lepton_E_7_d, tf_lepton_E_8_d,
     &     tf_lepton_E_9_d, tf_lepton_E_10_d, tf_lepton_E_11_d, tf_lepton_E_12_d
      double precision, constant :: mvir2_d(-11:24),  pmass_d(12)     
      double precision, constant :: Al, XVpow(0:96), XV(0:96), TV(0:20)
      real*4, constant :: UPD(15360) 
      integer, constant :: Nx, Nt, NfMx

      contains
      attributes(global) subroutine fct_cuda(xarray,wgtarray,resfctarray)
      implicit none

      double precision, dimension(npoints_d,ndim_d) :: xarray
      double precision, dimension(npoints_d) :: wgtarray,resfctarray

      double precision :: X1,X2,PSWGT,JAC
      integer id

c----------local variables -----------------
      double precision twgt
c     ----- global variables ----------------
c     get_PS_point
      integer n_var,  vis(2)
      double precision Emiss_weight
c     generate_visible
      double precision jac_temp,Emax,sqrts,jac_visible
      integer i,j,k,nu
      double precision Etot,pztot,misspx,misspy
      double precision, momenta_d(0:3,-11:24)
      double precision pexp(0:3,4)
c     class a
      integer p1,p2
      double precision pboost(0:3), CMS_mom(0:3,12)
      double precision Ptot(0:3),PtotCMS(0:3)
      double precision measureLAB, measureCMS
      double precision normp1,normp2,det,sqrts
      double precision angles(2,2),px(2),py(2),jac_loc,c1,c2,ss,xx
      integer MG,k
      double precision c3(0:3)

      double precision misspx_reco, misspy_reco
      double precision test2
      double precision pxISR, pyISR, gen_var(4,3)

      double precision xbk(2)
      double precision pp(0:3,4)
      double precision p1_l(0:3,4), xdum
      DOUBLE PRECISION PTEMP(0:3)
      integer imode, IMIRROR, iproc
      double precision ib_d(2), dsigproc_l, dsig_l,dsig1_l

      INTEGER LP
      DOUBLE PRECISION G1,G2
      DOUBLE PRECISION DSIGUU  
      DOUBLE PRECISION PD(0:2)

c     ********************* VARIABLES SMATRIX1 ************************************
c     
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=16)
      REAL*8 ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER NTRY
      REAL*8 T
      INTEGER IHEL,IDEN
      INTEGER JC(4)
      LOGICAL GOODHEL(NCOMB)

      iden = 256 
      
      pxISR=t_ISR(1)
      pyISR=t_ISR(2)
      pexp=pexp_d
      

c     CUDA number of thread
      id = threadidx%x + (blockidx%x-1)*blockdim%x
      if ( id <= npoints_d) then
c     -------------  call get_PS_point(x) --------- start
         jac=1d0/((2d0*pi)**(3*(nexternal-2)-4))
         
         n_var=0
c     -------------  call generate_visible(x,n_var) ------ start
         sqrts = dsqrt(s_d)
         Emax = sqrts
         jac_visible=1d0
         misspx=0d0
         misspy=0d0
         Etot=0d0
         pztot=0d0
         n_var=0                ! n_var labels the variables of integration
c     !   Rest cut out because num_vis(config_pos)=0

c     -------------  call generate_visible(x,n_var) ------ end
         
         if (missPT_EXP_d.ge.0d0) then ! using reconstructed missing pT to apply the boost correction
c     --------------- call generate_miss_parton --- start

            misspx_reco=missPT_EXP_d*dcos(missPhi_EXP_d)
            misspy_reco=missPT_EXP_d*dsin(missPhi_EXP_d)
            Emiss_weight=1d0
c     call smear_missing_reco(misspx_reco,misspy_reco,weight)
            pxISR =-misspx_reco -px_visible_d
            pyISR =-misspy_reco -py_visible_d
            misspx=misspx-pxISR
            misspy=misspy-pyISR
            jac=jac*Emiss_weight
         else
            pxISR=0d0
            pyISR=0d0
         endif
c     ************************************************************************
c     ------------ block a start
c     subroutine class_h(x,n_var,p1 :=3 ,p2 = 4)
c     ************************************************************************
         
         p1=3
         p2=4
         jac_loc=1d0
         test2=1d0
         sqrts=dsqrt(s_d)
         
         vis(1)=p1
         vis(2)=p2
         do i=1,2
            do j=1,3
c     
c     if width is zero, just take the exp. component (TF=delta function)
               if(c_point_d(vis(i),j,2).lt.1d-6) then
                  gen_var(i,j)=c_point_d(vis(i),j,1) 
c     
c     if width is positive, generate the component
               elseif(c_point_d(vis(i),j,2).gt.0d0) then

                  n_var=n_var+1 ! update the component of random variable

                  c1=c_point_d(vis(i),j,1)
                  c2=c_point_d(vis(i),j,2)
                  xx=xarray(id,n_var)
                  call get_component_cuda(c1,c2,
     &                 xx, gen_var(i,j),jac_temp,j,sqrts-Etot)

                  jac_loc=jac_loc*jac_temp
                  test2=test2*jac_temp
                  
               endif
            enddo
            
c-----------------------------------------------------------------
c     Now theta,phi and |p| of particle i (MG label) are defined.
c     define the momentum in a Lorentz fashion,
c     and record the result in momenta(#,i)
c------------------------------------------------------------------
            c1=pmass_d(vis(i))
            
            call four_momentum_cuda(gen_var(i,1),gen_var(i,2),gen_var(i,3),
     &           c1, momenta_d(0,vis(i)))

c----------------------------------------------
c     update missing transverse mome ntum
c----------------------------------------------
            misspx=misspx-momenta_d(1,vis(i))
            misspy=misspy-momenta_d(2,vis(i))
c----------------------------------------------
c     update Etot and pztot for visible particles
c----------------------------------------------
            Etot=Etot+momenta_d(0,vis(i))
            pztot=pztot+momenta_d(3,vis(i))
            Emax=sqrts-Etot
c----------------------------------------------
c     update jacobian
c----------------------------------------------
            
            
            jac_loc=jac_loc*gen_var(i,3)*gen_var(i,3)*dsin(gen_var(i,1))/(2d0*momenta_d(0,vis(i)))
c     write(89,*) 'momenta',i,c1
c     write(89,*) 'jac_loc after',i,jac_loc  
         enddo


c     Apply the boost correction
c     

c     First evaluated the total momentum in the LAB frame
         do j=0,3
            Ptot(j)=0d0
            do k=3,nexternal
               Ptot(j)=Ptot(j)+momenta_d(j,k)
            enddo
            pboost(j)=Ptot(j)
         enddo
         
c     Then calculate the momenta in the CMS frame
         pboost(1)=-pboost(1)
         pboost(2)=-pboost(2)
         pboost(3)=0d0
         do j=3,nexternal
c     write(*,*) "p",j,momenta(0,j), momenta(1,j),momenta(2,j),momenta(3,j)
            
            call boostx_cuda(momenta_d(0,j),pboost,CMS_mom(0,j))
         enddo
         call boostx_cuda(Ptot,pboost,PtotCMS)

c     Evaluate the initial momenta in the CMS frame
         x1=(PtotCMS(0)+PtotCMS(3))/sqrts
         x2=(PtotCMS(0)-PtotCMS(3))/sqrts

         CMS_mom(0,1)=sqrts*x1/2d0
         CMS_mom(1,1)=0d0
         CMS_mom(2,1)=0d0
         CMS_mom(3,1)=sqrts*x1/2d0
         CMS_mom(0,2)=sqrts*x2/2d0
         CMS_mom(1,2)=0d0
         CMS_mom(2,2)=0d0
         CMS_mom(3,2)=-sqrts*x2/2d0

c     Evaluate the initial momenta in the LAB frame
         pboost(1)=Ptot(1)
         pboost(2)=Ptot(2)
         
         call boostx_cuda(CMS_mom(0,1),pboost,momenta_d(0,1))
         
         call boostx_cuda(CMS_mom(0,2),pboost,momenta_d(0,2))


         measureLAB=1d0
         do j=3,nexternal-num_inv
            MG=inv_matching_type_part_d(j)
            measureLAB=measureLAB*dsqrt(momenta_d(1,MG)**2+momenta_d(2,MG)**2)
         enddo

         measureCMS=1d0
         do j=3,nexternal-num_inv
            MG=inv_matching_type_part_d(j)
            measureCMS=measureCMS*dsqrt(CMS_mom(1,MG)**2+CMS_mom(2,MG)**2)
         enddo

         jac=jac*measureCMS/measureLAB
c     
c     flux factor
c     
         jac_loc=jac_loc/(2d0*S_d*x1*x2) ! flux 
         jac=jac*jac_loc



c     ************************************************************************
c     ----------- block a end
c*************************************************************************

c     --------------- call generate_miss_parton --- end


c     call main_code(x,n_var):
c     call class_a(x,n_var,3,4)
         
c     -------------  call get_PS_point(x) --------- end

         


         if (jac.gt.0d0) then
            
            resfctarray(id)=jac
            xbk(1)=X1
            xbk(2)=X2
c$$$  resfct=resfct*dsig(momenta(0,1),wgt)
C     BEGIN CODE
C     ----------
            DSIG_l=0D0
c     change this one later
            call change_shape_cuda(momenta_d(0,1),pp)
            
            imode =1            ! Dummy variable
C     Select among the subprocesses based on PDF weight
            DO IPROC=1,2
               DO IMIRROR=1,2
                  IF(IMIRROR.EQ.1.OR.MIRRORPROCS_d(IPROC))THEN
                     
c     ------------ dsigproc start ------------------------------------------------------------------------------- end
                     DO I=1,NEXTERNAL
                        DO J=0,3
                           P1_l(J,I)=PP(J,I)
                        ENDDO
                     ENDDO
                     IB_d(1)=1
                     IB_d(2)=2
                     IF(IMIRROR.EQ.2)THEN
                        DO I=0,3
                           PTEMP(I)=P1_l(I,1)
                           P1_l(I,1)=P1_l(I,2)
                           P1_l(I,2)=PTEMP(I)
                        ENDDO
C     Flip x values 
                        XDUM=XBK(1)
                        XBK(1)=XBK(2)
                        XBK(2)=XDUM
                     ENDIF

                     DSIGPROC_l=0D0

                     IF(IPROC.EQ.2)THEN
c     **************    DSIGPROC_l=DSIG2(P1_l,WGT,IMODE) ! c c~ > mu+ mu-

c     ******************* START DSIG1 ***************************************
                        DSIG1_l=0D0
                        IF (ABS(LPP_d(IB_d(1))).GE.1) THEN
                           LP=SIGN(1,LPP_d(IB_d(1)))
                           G1=PDG2PDF_cuda(ABS(LPP_d(IB_d(1))),4*LP,XBK(IB_d(1)),DSQRT(Q2FACT_d(1)))
                           
                        ENDIF
                        IF (ABS(LPP_d(IB_d(2))).GE.1) THEN
                           LP=SIGN(1,LPP_d(IB_d(2)))
                           G2=PDG2PDF_cuda(ABS(LPP_d(IB_d(2))),-4*LP,XBK(IB_d(2)),DSQRT(Q2FACT_d(2)))
!G2=PDG2PDF_cuda(1,0,0.0277d0,91.188d0)
                           
                        ENDIF
                        PD(0) = 0D0
!IPROC = 0
!IPROC=IPROC+1 ! g g > mu+ mu-
!PD(IPROC)=G1*G2 
                        
                        PD(0)=PD(0)+DABS(G1*G2)
c     ******************* START SMATRIX1 ***************************************

c     CALL SMATRIX1(PP,DSIGUU)
                        NTRY=NTRY+1
                        DO IHEL=1,NEXTERNAL
                           JC(IHEL) = +1
                        ENDDO
                        DSIGUU = 0D0 
                        DO IHEL=1,NCOMB
                           IF (GOODHEL(IHEL) .LT. 2) THEN
                              T=MATRIX2_cuda(Pp ,NHEL(1,IHEL),JC(1))
                              DSIGUU=DSIGUU+T
!DSIGUU=dsiguu+1d0
                              
                              IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                                 GOODHEL(IHEL)=.TRUE.
                              ENDIF
                           ENDIF
                        ENDDO
                        DSIGUU=DSIGUU/36d0
                        
c     ******************* END SMATRIX1 *******************************************
                        IF (DSIGUU.LT.1D199) THEN
                           DSIG1_l=PD(0)*DSIGUU
                        ELSE
c     WRITE(*,*) 'Error in matrix element'
                           DSIGUU=0D0
                           DSIG1_l=0D0
                        ENDIF
                        dsigproc_l=dsig1_l
c     ******************* END DISG2 *******************************************
                     ENDIF

                     IF(IPROC.EQ.1)THEN
c     **************    DSIGPROC_l=DSIG1(P1_l,WGT,IMODE) ! g g > mu+ mu-

c     ******************* START DSIG1 ***************************************
                        DSIG1_l=0D0
                        IF (ABS(LPP_d(IB_d(1))).GE.1) THEN
                           LP=SIGN(1,LPP_d(IB_d(1)))
                           G1=PDG2PDF_cuda(ABS(LPP_d(IB_d(1))),0*LP,XBK(IB_d(1)),DSQRT(Q2FACT_d(1)))
                           
                        ENDIF
                        IF (ABS(LPP_d(IB_d(2))).GE.1) THEN
                           LP=SIGN(1,LPP_d(IB_d(2)))
                           G2=PDG2PDF_cuda(ABS(LPP_d(IB_d(2))),0*LP,XBK(IB_d(2)),DSQRT(Q2FACT_d(2)))
!G2=PDG2PDF_cuda(1,0,0.0277d0,91.188d0)
                           
                        ENDIF
                        PD(0) = 0D0
                        IPROC = 0
                        IPROC=IPROC+1 ! g g > mu+ mu-
                        PD(IPROC)=G1*G2 
                        
                        PD(0)=PD(0)+DABS(PD(IPROC))
c     ******************* START SMATRIX1 ***************************************

c     CALL SMATRIX1(PP,DSIGUU)
                        NTRY=NTRY+1
                        DO IHEL=1,NEXTERNAL
                           JC(IHEL) = +1
                        ENDDO
                        DSIGUU = 0D0 
                        DO IHEL=1,NCOMB
                           IF (GOODHEL(IHEL) .LT. 2) THEN
                              T=MATRIX1_cuda(Pp ,NHEL(1,IHEL),JC(1))
                              DSIGUU=DSIGUU+T
!DSIGUU=dsiguu+1d0
                              
                              IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
                                 GOODHEL(IHEL)=.TRUE.
                              ENDIF
                           ENDIF
                        ENDDO
                        DSIGUU=DSIGUU/DBLE(IDEN)
                        
c     ******************* END SMATRIX1 *******************************************
                        IF (DSIGUU.LT.1D199) THEN
                           DSIG1_l=PD(0)*DSIGUU
                        ELSE
c     WRITE(*,*) 'Error in matrix element'
                           DSIGUU=0D0
                           DSIG1_l=0D0
                        ENDIF
                        dsigproc_l=dsig1_l
c     ******************* END DISG1 *******************************************
                     ENDIF
c     INCLUDE AT LATER TIME *******************
c     *****    IF(IPROC.EQ.2) DSIGPROC_l=DSIG2(P1_l,WGT,IMODE) ! c c~ > mu+ mu-
c     -------------- dsigproc end --------------------------------------------------------------------------------- start
                     DSIG_l=DSIG_l+DSIGPROC_l !(pp,IPROC,IMIRROR)
                     IF(IMIRROR.EQ.2)THEN
C     Need to flip back x values
                        XDUM=XBK(1)
                        XBK(1)=XBK(2)
                        XBK(2)=XDUM
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            
            resfctarray(id)=resfctarray(id)*dsig_l
c     -------- CRUCIAL ----------------------------------- !!!!!!!!!!!!!!!!!!!!!!!

c$$$  call transfer_fct(momenta(0,1),TWGT)
            twgt=1d0
            call tf_E_lepton_cuda(pexp(0,3),  momenta_d(0,3),twgt)
            call tf_THETA_lepton_cuda(pexp(0,3),  momenta_d(0,3),twgt)
            call tf_PHI_lepton_cuda(pexp(0,3),  momenta_d(0,3),twgt)

            call tf_E_lepton_cuda(pexp(0,4),  momenta_d(0,4),twgt)
            call tf_THETA_lepton_cuda(pexp(0,4),  momenta_d(0,4),twgt)
            call tf_PHI_lepton_cuda(pexp(0,4),  momenta_d(0,4),twgt)

            resfctarray(id)=resfctarray(id)*twgt*wgtarray(id)
         endif
         
         
      endif
      end subroutine fct_cuda

      attributes(device) subroutine boostx_cuda(p,q,pboost)
      implicit none
      double precision p(0:3),q(0:3),pboost(0:3),pq,qq,m,lf

      double precision rZero
      parameter( rZero = 0.0d0 )
      qq = q(1)**2+q(2)**2+q(3)**2
      if ( qq.ne.rZero ) then
         pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m = sqrt(max(q(0)**2-qq,1d-99))
         lf = ((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
      else
         pboost(0) = p(0)
         pboost(1) = p(1)
         pboost(2) = p(2)
         pboost(3) = p(3)
      endif
c     
      return
      end subroutine boostx_cuda

      attributes(device) subroutine four_momentum_cuda(theta,phi,rho,m,p)
      implicit none
      double precision theta,phi,rho,m,p(0:3)

      P(1)=rho*dsin(theta)*dcos(phi)
      P(2)=rho*dsin(theta)*dsin(phi)
      P(3)=rho*dcos(theta)
      P(0)=dsqrt(rho**2+m**2)

      return
      end subroutine four_momentum_cuda

      attributes(device)  subroutine get_component_cuda(c_point,gam,x,gen_point,jac,var_num,Emax)
      implicit none
      double precision c_point,gam,x,Emax
      double precision gen_point,jac
      integer, intent(in) :: var_num
      double precision t1, t2
c     
c     local
c     
      double precision point_min, point_max
c     
c     Parameter
c     
      double precision pi,zero
      parameter (pi=3.141592654d0,zero=0d0)

      t1 = 2d0
      t2 = 3d0

c     var_num=1 means that we generate a theta
      if (var_num.eq.1) then

         point_max=c_point +5d0*gam
         if (point_max.gt.pi) point_max=pi
         point_min=c_point -5d0*gam
         if (point_min.lt.0d0) point_min=0d0

         gen_point=(point_max-point_min)*x+point_min

c     var_num=2 means that we generate a phi (note that phi is a cyclic variable) 
      elseif(var_num.eq.2) then

         if(gam.lt.(2d0*pi/10)) then
            point_max=c_point +5d0*gam
            point_min=c_point -5d0*gam
         else
            point_max=2*pi
            point_min=0d0
         endif
         t1=((point_max-point_min)*x+point_min)
         t2=2d0*pi
         gen_point=dble(t1-floor(t1/t2)*t2)
c     gen_point=dble(mod(((point_max-point_min)*x+point_min),2d0*pi))
         if(gen_point.lt.zero) then
            gen_point=gen_point+2d0*pi ! this is true since phi is cyclic
         endif



c     var_num=3 means that we generate a rho
      elseif (var_num.eq.3) then
         point_max=dble(min(c_point +5d0*gam,Emax))
         point_min=dble(max(c_point -5d0*gam,0.d0))
         if (point_max.le.point_min) then
            jac=-1d0
            return
         endif
         gen_point=(point_max-point_min)*x+point_min

      endif

c**************************************************************
c     III.  compute the jacobian                              *
c**************************************************************

      jac=(point_max-point_min)

      return
      end subroutine get_component_cuda

      attributes(device) subroutine tf_E_lepton_cuda(pexp,p,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'


      prov1=(tf_lepton_E_7_d+tf_lepton_E_8_d*dsqrt(p(0))+tf_lepton_E_9_d*p(0))
      prov2=(tf_lepton_E_10_d+tf_lepton_E_11_d*dsqrt(p(0))+tf_lepton_E_12_d*p(0))

      tf=(exp(-(p(0)-pexp(0)-prov1)**2/2d0/prov2**2)) !first gaussian
      tf=tf*((1d0/dsqrt(2d0*pi))/(prov2)) !normalisation



      weight=weight*tf

      return
      end subroutine tf_E_lepton_cuda

      attributes(device) subroutine tf_THETA_lepton_cuda(pexp,p,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

      tf=1d0


      weight=weight*tf

      return
      end subroutine tf_THETA_lepton_cuda

      attributes(device) subroutine tf_PHI_lepton_cuda(pexp,p,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

      tf=1d0


      weight=weight*tf

      return
      end  subroutine tf_PHI_lepton_cuda

      attributes(device) subroutine change_shape_cuda(p1,p2)
      implicit none

      double precision p1(0:3,4)
      double precision p2(0:3,4)

      p2=p1

      return
      end  subroutine change_shape_cuda


      attributes(device) double precision function pdg2pdf_cuda(ih,ipdg,x,xmu)
c***************************************************************************
c     Based on pdf.f, wrapper for calling the pdf of MCFM
c***************************************************************************
      implicit none
c     
c     Arguments
c     
      DOUBLE  PRECISION x,xmu
      INTEGER IH,ipdg

      double precision Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf,Ctq5L
      integer mode,Irt,i,j
      double precision xlast(2),xmulast(2),pdflast(-7:7,2),q2max
      double precision epa_electron,epa_proton,Ctq6Pdf
      integer ipart,ireuse,iporg

! ******************** REMOVE DATA!!! ****************************** check save!!
!save xlast,xmulast,pdflast
!data xlast/2*0d0/
!data pdflast/30*0d0/

      ipart=ipdg
      if(iabs(ipart).eq.21) ipart=0
      if(iabs(ipart).eq.22) ipart=7
      iporg=ipart

c     This will be called for any PDG code, but we only support up to 7
      if(iabs(ipart).gt.7)then
         pdg2pdf_cuda=0d0
         return 
      endif

      ireuse = 0
      do i=1,2
c     Check if result can be reused since any of last two calls
         if (x.eq.xlast(i).and.xmu.eq.xmulast(i)) then
            ireuse = i
         endif
      enddo

c     If both x non-zero and not ireuse, then zero x
      if (ireuse.eq.0.and.xlast(1).ne.0d0.and.xlast(2).ne.0d0)then
         do i=1,2
            xlast(i)=0d0
            xmulast(i)=0d0
            do j=-7,7
               pdflast(j,i)=0d0
            enddo
         enddo
         ireuse=1
      else if(ireuse.eq.0.and.xlast(1).ne.0d0)then
         ireuse=2
      else if(ireuse.eq.0)then
         ireuse=1
      endif


      xlast(ireuse)=x
      xmulast(ireuse)=xmu

c$$$  if(iabs(ipart).eq.7.and.ih.gt.1) then
c$$$  q2max=xmu*xmu
c$$$  if(ih.eq.3) then       !from the electron
c$$$  pdg2pdf=epa_electron(x,q2max)
c$$$  elseif(ih .eq. 2) then !from a proton without breaking
c$$$  pdg2pdf=epa_proton(x,q2max)
c$$$  endif 
c$$$  pdflast(iporg,ireuse)=pdg2pdf
c$$$  return
c$$$  endif

      if(iabs(ipart).ge.1.and.iabs(ipart).le.2)
     $     ipart=sign(3-iabs(ipart),ipart)

      Ctq6Pdf = PartonX6_cuda (ipart, x, xmu)
      if(Ctq6Pdf.lt.0.D0)  Ctq6Pdf = 0.D0
      pdg2pdf_cuda=Ctq6Pdf
      
c     pdflast(iporg,ireuse)=pdg2pdf
      return
      end function pdg2pdf_cuda

      
      attributes(device) SUBROUTINE POLINT6_cuda (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Adapted from "Numerical Recipes"
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
         DIFT=ABS(X-XA(I))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         ENDIF
         C(I)=YA(I)
         D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) go to 12 ! instead of STOP, possible bug
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END subroutine polint6_cuda


         attributes(device) function partonx6_cuda (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      Parameter (MXX = 96, MXQ = 20, MXF = 5)
      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))


      double precision fvec(4), fij(4)
      double precision xpow,onep
      integer nqvec
      !Data OneP / 1.00001 
      onep=1.00001d0
      !Data xpow / 0.3d0 /       !**** choice of interpolation variable
      xpow=0.3d0
      !Data nqvec / 4 /
      nqvec=4
   
             
      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Go to 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        ff=0d0
        
        go to 17
      ElseIf (JLx .Eq. 0) Then    
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
           
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        ff=0d0
        
        go to 17
      Endif
C     ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C     This is the variable to be interpolated in
      ss = x**xpow
      
      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
         svec1 = xvpow(jx)
         svec2 = xvpow(jx+1)
         svec3 = xvpow(jx+2)
         svec4 = xvpow(jx+3)

         s12 = svec1 - svec2
         s13 = svec1 - svec3
         
         s23 = svec2 - svec3
         s24 = svec2 - svec4
         s34 = svec3 - svec4

         sy2 = ss - svec2
         sy3 = ss - svec3

c     constants needed for interpolating in s at fixed t lattice points...
         const1 = s13/s23
         const2 = s12/s23
         const3 = s34/s23
         const4 = s24/s23
         s1213 = s12 + s13
         s2434 = s24 + s34
         sdet = s12*s34 - s1213*s2434
         tmp = sy2*sy3/sdet
         const5 = (s34*sy2-s2434*sy3)*tmp/s12
         const6 = (s1213*sy2-s12*sy3)*tmp/s34
      EndIf

c     --------------Now find lower end of interval containing Q, i.e.,
c     get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C     keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C     JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
         Jq = Nt - 3

      Endif
C     This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c     store the lattice points in t...
         tvec1 = Tv(jq)
         tvec2 = Tv(jq+1)
         tvec3 = Tv(jq+2)
         tvec4 = Tv(jq+3)

         t12 = tvec1 - tvec2
         t13 = tvec1 - tvec3
         t23 = tvec2 - tvec3
         t24 = tvec2 - tvec4
         t34 = tvec3 - tvec4

         ty2 = tt - tvec2
         ty3 = tt - tvec3

         tmp1 = t12 + t13
         tmp2 = t24 + t34

         tdet = t12*t34 - tmp1*tmp2

      EndIf


c     get the pdf function values at the lattice points...

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec
         J1  = jtmp + it*(NX+1)
         

         If (Jx .Eq. 0) Then
C     For the first 4 x points, interpolate x^2*f(x,Q)
C     This applies to the two lowest bins JLx = 0, 1
C     We can not put the JLx.eq.1 bin into the "interrior" section
C     (as we do for q), since Upd(J1) is undefined.
            fij(1) = 0
            fij(2) = dble(Upd(J1+1)) * XV(1)**2
            fij(3) = dble(Upd(J1+2)) * XV(2)**2
            fij(4) = dble(Upd(J1+3)) * XV(3)**2
C     
C     Use Polint6 which allows x to be anywhere w.r.t. the grid

            Call Polint6_cuda (XVpow(0), Fij(1), 4, ss, Fx, Dfx)
            
            If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
            
C     Pdf is undefined for x.eq.0
         ElseIf  (JLx .Eq. Nx-1) Then
C     This is the highest x bin:

            Call Polint6_cuda (XVpow(Nx-3), dble(upd(j1)), 4, ss, Fx, Dfx)

            Fvec(it) = Fx

         Else
C     for all interior points, use Jon's in-line function
C     This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
            sf2 = dble(upd(j1+1))
            sf3 = dble(upd(j1+2))

            g1 =  sf2*const1 - sf3*const2
            g4 = -sf2*const3 + sf3*const4
            
            Fvec(it) = (const5*(dble(upd(j1))-g1)
     &           + const6*(dble(upd(j1+3))-g4)
     &           + sf2*sy3 - sf3*sy2) / s23
            

         Endif

      enddo
C     We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C     1st Q-bin, as well as extrapolation to lower Q
         Call Polint6_cuda (TV(0), Fvec(1), 4, tt, ff, Dfq)
         

      ElseIf (JLq .GE. Nt-1) Then
C     Last Q-bin, as well as extrapolation to higher Q
         Call Polint6_cuda (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
         
      Else
C     Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C     which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C     the full range QV(0:Nt)  (in contrast to XV)
         tf2 = fvec(2)
         tf3 = fvec(3)

         g1 = ( tf2*t13 - tf3*t12) / t23
         g4 = (-tf2*t34 + tf3*t24) / t23

         h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &        +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)
         
         ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
         
      EndIf
      
 17   PartonX6_cuda = ff
      Return
c     ********************
      End function partonx6_cuda



      attributes(device) REAL*8 FUNCTION MATRIX1_cuda(P,NHEL,IC)
C     
C     Generated by MadGraph 5 v. 2.0.0.beta1, 2012-11-08
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g g > h > mu+ mu- WEIGHTED=4 HIG=1 HIW=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=5, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      
      
C     BEGIN CODE
C     ----------
      CALL VXXXXX_cuda(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX_cuda(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX_cuda(P(0,3),MM_d,NHEL(3),-1*IC(3),W(1,3))
      CALL OXXXXX_cuda(P(0,4),MM_d,NHEL(4),+1*IC(4),W(1,4))
      CALL VVS3_3_cuda(W(1,1),W(1,2),GC_8_d,MH_d,WH_d,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL FFS1_0_cuda(W(1,3),W(1,4),W(1,5),GC_42_d,AMP(1))
      JAMP(1)=+2D0*(+AMP(1))

      
      
      ZTEMP = (0.D0,0.D0)
      
      ZTEMP = ZTEMP + 2d0*JAMP(1)
      
      MATRIX1_cuda = ZTEMP*DCONJG(JAMP(1)) !/DENOM(I)
      
      END FUNCTION MATRIX1_cuda

      attributes(device) REAL*8 FUNCTION MATRIX2_cuda(P,NHEL,IC)
C     
C     Generated by MadGraph 5 v. 2.0.0.beta1, 2012-11-08
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g g > h > mu+ mu- WEIGHTED=4 HIG=1 HIW=1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=5, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1

C     ----------
      CALL IXXXXX_cuda(P(0,1),MC_d,NHEL(1),+1*IC(1),W(1,1))
      CALL OXXXXX_cuda(P(0,2),MC_d,NHEL(2),-1*IC(2),W(1,2))
      CALL IXXXXX_cuda(P(0,3),MM_d,NHEL(3),-1*IC(3),W(1,3))
      CALL OXXXXX_cuda(P(0,4),MM_d,NHEL(4),+1*IC(4),W(1,4))
      CALL FFS1_3_cuda(W(1,1),W(1,2),GC_40_d,MH_d,WH_d,W(1,5))
C     Amplitude(s) for diagram number 1
      CALL FFS1_0_cuda(W(1,3),W(1,4),W(1,5),GC_42_d,AMP(1))
      JAMP(1)=+AMP(1)

      
      
      ZTEMP = (0.D0,0.D0)
      
      ZTEMP = ZTEMP + 3d0*JAMP(1)
      
      MATRIX2_cuda = ZTEMP*DCONJG(JAMP(1)) !/DENOM(I)
      
      END FUNCTION MATRIX2_cuda

      end module madweightcuda

      
