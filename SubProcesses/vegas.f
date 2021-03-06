
C****************************************************************************
C     THIS IS THE SOUBROUTINE VEGAS for CUDA Fortran
C                   EDITED BY MAX NOLTE
C****************************************************************************
C     NCALL IS THE NUMBER OF CALLS TO VEGAS.
C     NPRN >  0 VEGAS PRINTS THE RESULTS OF EACH ITERATION.
C     NPRN 0 VEGAS PRINTS NOTHING.
C     NPRN < 0 VEGAS PRINTS ALL.
C     XL(I) IS LOWER INTEGRATION LIMIT ON I TH AXIS.
C     XU(I) IS UPPER INTEGRATION LIMIT ON I THE AXIS.
     
      subroutine vegas(fxn,avgi,sd,chi2a)
      use madweightcuda
      use cudafor 
 
c
c     Routine performs n-dimensional Monte Carlo integration
c     Originally written by P. Lepage
c     Modified by Max Nolte for CUDA Fortran
c
c
      implicit double precision (a-h,o-z)
     
      implicit integer*4 (i-n)
      double precision time(20)
c ----- new arrays for cuda ------------
      double precision, dimension(:), allocatable :: wgtarray
      double precision, dimension(:), allocatable :: farray
      double precision, dimension(:,:), allocatable :: xarray
      integer, dimension(:,:), allocatable :: iaarray
      double precision, device, dimension(:), allocatable :: wgtarray_d
      double precision, device, dimension(:), allocatable :: farray_d
      double precision, device, dimension(:,:), allocatable :: xarray_d
      type(dim3) :: grid, block
      double precision, target, device :: at(1)
c ----- end new arrays for cuda ------
      common/bveg1/xl(20),xu(20),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,20),si,si2,swgt,schi,ndo,it
      common/bveg3/alph,ndmx,mds
      common/bveg4/calls,ti,tsi
      COMMON/SEED/NUM,NUM2
      dimension d(50,20),di(50,20),xin(50),r(50),dx(20),dt(20),
     1     x(20),kg(20),ia(20)
      dimension RAND(20)
      DATA NPRN/0/,NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/1/
      DATA XL/20*0.d0/,XU/20*1.d0/
      integer flagTime, index, npoints, allocatestatus
      common /to_seed/ iseed
c     channel position (phasespace.inc) !!! TO DO
      integer config_pos,perm_pos
      common /to_config/config_pos,perm_pos
      double precision momenta(0:3,-11:24)
      double precision mvir2(-11:24)    
      common /to_diagram_kin/ momenta, mvir2
      double precision  missPhi_EXP, missPT_EXP
      common /to_missEXP/  missPhi_EXP, missPT_EXP
      double precision pxISR, pyISR
      common /to_ISR/  pxISR, pyISR
      double precision              S,X1,X2,PSWGT,JAC
      common /PHASESPACE/ S,X1,X2,PSWGT,JAC
      double precision px_visible,py_visible
      common /to_pTrec_visible/px_visible,py_visible
      double precision c_point(1:12,3,2)
      common/ph_sp_init/c_point
      double precision pmass(12)    
      common / to_mass/pmass
      integer nexternal2, num_inv2
      COMMON/to_num_inv/nexternal2, num_inv2
      integer matching_type_part(3:12) !modif/link between our order by type for permutation
      integer inv_matching_type_part(3:12)
      common/madgraph_order_type/matching_type_part, inv_matching_type_part
      double precision pexp(0:3,4)
      common/to_pexp/pexp
      double precision  tf_lepton_E_7
      double precision  tf_lepton_E_8
      double precision  tf_lepton_E_9
      double precision  tf_lepton_E_10
      double precision  tf_lepton_E_11
      double precision  tf_lepton_E_12
      double precision  tf_jet_E_1
      double precision  tf_jet_E_2
      double precision  tf_jet_E_3
      double precision  tf_jet_E_4
      double precision  tf_jet_E_5
      double precision  tf_jet_E_6
      Common/to_TF_param/tf_lepton_E_7,tf_lepton_E_8,
     &tf_lepton_E_9,tf_lepton_E_10,tf_lepton_E_11,tf_lepton_E_12,
     &tf_jet_E_1,tf_jet_E_2,tf_jet_E_3,tf_jet_E_4,tf_jet_E_5,tf_jet_E_6
      LOGICAL MIRRORPROCS(2)
      INCLUDE 'mirrorprocs.inc'
      integer  lpp(2)
      double precision    ebeam(2), xbk(2),q2fact(2)
      common/to_collider/ ebeam   , xbk   ,q2fact,   lpp

      DOUBLE PRECISION ME,MC,MB,MM,MH,MT,MW,MP,MTA,MZ
      COMMON/MASSES/ ME,MC,MB,MM,MH,MT,MW,MP,MTA,MZ
      DOUBLE PRECISION WTAU,WH,WW,WT,WH1,WZ
      COMMON/WIDTHS/ WTAU,WH,WW,WT,WH1,WZ
      DOUBLE COMPLEX GC_8, GC_40, GC_42
      COMMON/COUPLINGS/ GC_8, GC_40, GC_42
      INTEGER NHEL_h(4,16)
    
c     for parton6x (cteq)

      real*4 UPD_hr(15360)

      Common
     > / CtqPar1 / Al_h, XV_h(0:96), TV_h(0:20), UPD_h(15360)
     > / CtqPar2 / Nx_h, Nt_h, NfMx_h
      double precision xvpow_h(0:96)
      double precision xpow 
      
      call cpu_time(time(1))
      write(89,*) 'Point 1 (start Vegas.f) ',time(1)
      nop = cudaThreadSynchronize()
      call cpu_time(time(2))
      write(89,*) 'Point 2 (1st cuda call) ',time(2)
      t_ISR(1)=pxISR
      t_ISR(2)=pyISR
   
      px_visible_d=px_visible
      py_visible_d=py_visible

c ---------- CUDA variables TEXTURE --------------- end

      nhel_h=reshape((/ -1,-1,-1,-1,-1,-1,-1, 1,-1,-1, 1,-1,-1,-1, 1, 1,-1, 1,-1,-1,-1, 1,-1,
     & 1,-1, 1, 1,-1,-1, 1, 1, 1,1,-1,-1,-1,1,-1,-1, 1,1,-1, 1,-1,1,-1, 1, 
     & 1,1, 1,-1,-1,1, 1,-1, 1,1, 1, 1,-1,1, 1, 1, 1 /), shape(nhel))
      nhel=nhel_h

c --------- CUDA variables CONSTANT -------------- start
      pi=3.141592653589793d0
      ndim_d = ndim
      s_d=s
c     phasespace.inc
      max_particles=12
      max_branches=max_particles-1
      max_configs=10
      max_channel=1000
c     nexternal.inc
      nexternal=4
      pexp_d=pexp
      nincoming=2
      missPhi_EXP_d=missPhi_EXP
      missPT_EXP_d=missPT_EXP
      mvir2_d=mvir2
      c_point_d=c_point
      pmass_d=pmass
      num_inv=num_inv2
      inv_matching_type_part_d=inv_matching_type_part
      mirrorprocs_d=mirrorprocs
      lpp_d=lpp
      q2fact_d=q2fact

      mm_d=mm
      mc_d=mc
      mh_d=mh
      wh_d=wh
      gc_8_d=gc_8
      GC_42_d=GC_42
      GC_40_d=GC_40


      tf_lepton_E_7_d=tf_lepton_E_7
      tf_lepton_E_8_d=tf_lepton_E_8
      tf_lepton_E_9_d=tf_lepton_E_9
      tf_lepton_E_10_d=tf_lepton_E_10
      tf_lepton_E_11_d=tf_lepton_E_11
      tf_lepton_E_12_d=tf_lepton_E_12


      do i=1,15360
         upd_hr(i)=real(UPD_h(i))
      enddo  
      upd=upd_hr
   
      Al=Al_h
      XV=XV_h

      xvpow_h(0) = 0D0
      xpow=0.3d0
      do i = 1, nx_h
        xvpow_h(i) = xv_h(i)**xpow
      enddo
      xvpow=xvpow_h

      
      TV=TV_h
      Nx=Nx_h
      Nt=Nt_h
      NfMx=NfMx_h
     
      nop = cudaThreadSynchronize()
      call cpu_time(time(3))
      write(89,*) 'Point 3 (constants transferred) ',time(3)
      write(89,*) 'constant transfer duration ',time(3)-time(2)
     

c ---------- CUDA variables CONSTANT --------------- end
      

      NUM=12345
      NUM2=67890
      flagTime=0 
      
c
      ndo=1
      do 1 j=1,ndim
 1       xi(1,j)=one
c     
      entry vegas1(fxn,avgi,sd,chi2a)
c     initialises  cumulative  variables but not grid
      it=0
      si=0.
      si2=si
      swgt=si
      schi=si
c
      entry vegas2(fxn,avgi,sd,chi2a)
c     no initialisation
      nd=ndmx
      ng=1
      if(mds.eq.0)go to 2
      ng=(ncall/2.)**(1./ndim)
      mds=1
      if((2*ng-ndmx).lt.0)go to 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
 2    k=ng**ndim
      npg=ncall/k
      if(npg.lt.2)npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one/calls
      do 3 j=1,ndim
         dx(j)=xu(j)-xl(j)
 3       xjac=xjac*dx(j)
c
c    rebin preserving bin density
c
      if(nd.eq.ndo)go to 8
      rc=ndo/xnd 
      do 7 J=1,ndim 
         k=0
         xn=0.
         dr=xn
         i=k
 4       k=k+1
         dr=dr+one
         xo=xn
         xn=xi(k,j)
 5       if(rc.gt.dr)go to 4
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr
         if(i.lt.ndm)go to 5
         do 6 i=1,ndm
 6          xi(i,j)=xin(i)
 7       xi(nd,j)=one
      ndo=nd
c
 8    if(nprn.ne.0)write(6,200)ndim,calls,it,itmx,acc
     1     ,mds,nd,(xl(j),xu(j),j=1,ndim)
c
      entry vegas3(fxn,avgi,sd,chi2a)
c     change seed
      NUM=NUM+iseed

c ----- Allocate cuda arrays ----------------------------------------------
c
            npoints=ng**ndim
            npoints=npg*npoints
c$$$            write(89,*) 'ncall ',ncall
c$$$            write(89,*) 'calls ',calls
c$$$            write(89,*) 'npoints ',npoints
            npoints_d = npoints


      allocate(wgtarray(npoints))
      allocate(farray(npoints))
      allocate(xarray(npoints,ndim))
      allocate(wgtarray_d(npoints))
      allocate(farray_d(npoints))
      allocate(xarray_d(npoints,ndim))
      allocate(iaarray(npoints,ndim))
      grid = dim3(ceiling(real(npoints)/256),1,1)
      block = dim3(256,1,1)

      nop = cudaThreadSynchronize()
      call cpu_time(time(4))
      write(89,*) 'Point 4 (start vegas algorithm) ',time(4)
      
c     main integration loop ----- Start of integration ------- i -----
 9    it=it+1
      ti=0.
      tsi=ti
      do 10 j=1,ndim
         kg(j)=1
         do 10 i=1,nd
            d(i,j)=ti
 10         di(i,j)=ti


       flagTime = flagTime+1
 

      index=0
 2511   do 2501 knpg=1,npg
           index = index+1
      call randa(ndim,rand)
      wgt=xjac
      do 15 j=1,ndim
         xn=(kg(j)-rand(j))*dxg+one
         ia(j)=xn
         iaarray(index,j)=ia(j)
         if(ia(j).gt.1)go to 13
         xo=xi(ia(j),j)
         rc=(xn-ia(j))*xo
         go to 14
 13      xO=xi(ia(j),j)-xi(ia(j)-1,j)
         rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14      x(j)=xl(j)+rc*dx(j)
          
         xarray(index,j)=x(j)
 15      wgt=wgt*xo*xnd
         wgtarray(index)=wgt
         f=wgt
       
 2501    farray(index)=f

          k=ndim
 2519                kg(k)=mod(kg(k),ng)+1
                   if(kg(k).ne.1)go to 2511
      k=k-1
      if(k.gt.0)go to 2519

      nop = cudaThreadSynchronize()
      call cpu_time(time(5))
      write(89,*) 'Point 5 (preperation done) ',time(5)

c     Move data go GPU
         wgtarray_d=wgtarray
         farray_d=farray
         xarray_d=xarray
c     Integrand function call
         nop = cudaThreadSynchronize()
         call cpu_time(time(6))
         write(89,*) 'Data transfer time ',time(6)-time(5)
         nop = cudaThreadSynchronize()
         call cpu_time(time(7))
         call fct_cuda<<<grid,block>>>(xarray_d,wgtarray_d,farray_d)
         nop = cudaThreadSynchronize()
         call cpu_time(time(8))
         write(89,*) 'Kernel execution time',time(8)-time(7)
         
c     Retrieve results from GPU
         farray=farray_d
         nop = cudaThreadSynchronize()
         call cpu_time(time(9))
         write(89,*) 'Results transfer time',time(9)-time(8)
         write(89,*) 'Point 6 (start vegas rest)',time(9)
          do 2510 j=1,ndim
 2510        kg(j)=1

c ---- 2nd block start ----
         index = 0
 11      fb=0.
         f2b=fb
          do 16 knpg=1,npg
             index = index+1

c TO DO transfer to GPU! next two lines
             f2=farray(index)*farray(index)
             fb=fb+farray(index)
             f2b=f2b+f2
             do 16 j=1,ndim
                di(iaarray(index,j),j)=di(iaarray(index,j),j)+farray(index)
 16             if(mds.ge.0)d(iaarray(index,j),J)=d(iaarray(index,j),J)+f2


 888            FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
                f2b= sqrt(f2b*      NPG)
                f2b=(f2b-fb)*(f2b+fb)
 1661           FORMAT(1X,'F2B',G14.6,'NPG',  I10)
                ti=ti+fb
                tsi=tsi+f2b
 33             FORMAT(1X,'TSI',G14.6,'F2B',G14.6)

              
                if(mds.ge.0)go to 18
                do 17 j=1,ndim
 17                d(iaarray(index,j),j)=d(iaarray(index,j),j)+f2b
 18                k=ndim
 19                kg(k)=mod(kg(k),ng)+1
                   if(kg(k).ne.1)go to 11
                   k=k-1
                   if(k.gt.0)go to 19
c
               
c final results for this iteration
c


      tsi=tsi*dv2g
      ti2=ti*ti
 88   format(1x,'tsi',g14.6)
      wgt=ti2/tsi
      si=si+ti*wgt
      si2=si2+ti2
      swgt=swgt+wgt
      schi=schi+ti2*wgt
 995  FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
      avgi=si/swgt
      sd=swgt*it/si2
      chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
      sd=dsqrt(one/sd)

    

c
      if(nprn.eq.0)go to 21
      tsi=dsqrt(tsi)
      write(6,201)it,ti,tsi,avgi,sd,chi2a
      if(nprn.ge.0)go to 21
      do 20 j=1,ndim
 20      write(6,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)
c
c ------  refine grid  ------------> then next iteration
c
 21   do 23 j=1,ndim
         xo=d(1,j)
         xn=d(2,j)
         d(1,j)=(xo+xn)/2.
         dt(j)=d(1,j)
         do 22 i=2,ndm
            d(i,j)=xo+xn
            xo=xn
            xn=d(i+1,j)
            d(i,j)=(d(i,j)+xn)/3.
 22         dt(j)=dt(j)+d(i,j)
         d(nd,j)=(xn+xo)/2.
 23      dt(j)=dt(j)+d(nd,j)
c
      do 28 j=1,ndim
         rc=0.
         do 24 i=1,nd
            r(i)=0.
            if(d(i,j).le.0.)go to 24
            xo=dt(j)/d(i,j)
            r(i)=((xo-one)/xo/dlog(xo))**alph
 24         rc=rc+r(i)
      rc=rc/xnd
      k=0
      xn=0.
      dr=xn
      i=k
 25   k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
 26   if(rc.gt.dr)go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
      if(i.lt.ndm)go to 26
      do 27 i=1,ndm
 27      xi(i,j)=xin(i)
 28   xi(nd,j)=one
c ------ check if accuracy reached or max dimensions -> next iteration (line 9)

   


      if(it.lt.itmx.and.acc*dabs(avgi).lt.sd)go to 9

        
      deallocate(wgtarray)
      deallocate(farray)
      deallocate(xarray)
      deallocate(wgtarray_d)
      deallocate(farray_d)
      deallocate(xarray_d)
      deallocate(iaarray)

 200  format(1X,'0input parameters for vegas:  ndim=',i3,
     1     '   ncall=',f8.0/28x,'  it=',i5,'    itmx=',i5/28x,
     2     '  acc=',g9.3/28x,'  mds=',i3,'     nd=',i4/28x,
     3     '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
 201  format(///' integration by vegas' / '0iteration no.',i5,
     1     ':  integral=',g14.8/21x,'std dev =',g10.4 /
     2     ' accumulated results:   integral=',g14.8/
     3     24x,'std dev =',g10.4 / 24x,'chi**2 per it''n =',g10.4)
 202  format(1X,'0data for axis',i2,/,' ',6x,'x',7x,'  delt i ',
     1     2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  '
     2     ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/,
     3     (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
      return
      entry vegas4(fxn,avgi,sd,chi2a)
      avgi=si/swgt
      sd=swgt*it/si2
      chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
      sd=dsqrt(one/sd)
      if(nprn.ne.0) write(6,201)it,0.d0,0.d0,avgi,sd,chi2a
      return
      end

*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION random2(iseed1,iseed2)
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION random2
      INTEGER ISEED1,iseed2,K,IZ
      K=ISEED1/53668
      ISEED1=40014*(ISEED1-K*53668)-K*12211
      IF (ISEED1.LT.0) ISEED1=ISEED1+2147483563
      K=ISEED2/52774
      ISEED2=40692*(ISEED2-K*52774)-K*3791
      IF (ISEED2.LT.0) ISEED2=ISEED2+2147483399
      IZ=ISEED1-ISEED2
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANDOM2=DBLE(IZ)/2147483589
      end

      FUNCTION RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
c      RANDOM = SEED*MINV
      RANDOM = SEED
      random = random/m
      END

      subroutine randa(n,rand)
      implicit double precision (a-h,o-z)
      COMMON/SEED/NUM,NUM2
      dimension rand(20)
      do 1 i=1,n
      rand(i)=random2(NUM,num2)
1     continue
      return
      end
C
C
    
