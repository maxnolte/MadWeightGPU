      PROGRAM DRIVER
c     use cudafor
C**************************************************************************
C     THIS IS THE DRIVER FOR THE WHOLE CALCULATION
C     
C**************************************************************************
C**************************************************************************
C     
C     Structure of the driver
C     
C     1) initialisation
C     2) first loop
C     2.1) loop on permutation
C     2.2) loop on configuration
C     3) second loop
C     3.1) loop on permutation
C     3.2) loop on configuration
C     3.3) write result
C     
C**************************************************************************
      IMPLICIT NONE
C     
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'data.inc'
      include 'madweight_param.inc'
C     
C     PARAMETERS
C     
      character*1 integrator
      parameter (integrator='v') ! v= vegas, m=mint
C     
C     LOCAL (Variable of integration)
C     
      DOUBLE PRECISION SD,CHI,CROSS
      INTEGER I,convergence_status
      integer pos,channel_pos,ll
      double precision temp_err, temp_val
      double precision order_value(max_channel)
      double precision order_error(max_channel)
      double precision normalize_perm
      real time,time1,time2,time3,time4
c     
c     variable for permutation, channel
c     
      integer perm_id(nexternal-2) !permutation of 1,2,...,nexternal-2
      integer num_per           !total number of permutations
      double precision weight,weight_error,final
      integer matching_type_part(3:max_particles) !modif/link between our order by type for permutation
      integer inv_matching_type_part(3:max_particles)
      common/madgraph_order_type/matching_type_part, inv_matching_type_part
c     
c     variable for selecting perm in loop 1
c     
      double precision minimal_value1,minimal_value2
      double precision best_precision, chan_val, chan_err
      double precision check_value
      integer loop_index,counter
c     integer best_config
c     
c     logical
c     
      logical debug

C     
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C     
      DOUBLE PRECISION Xl(20),XU(20),ACC
      INTEGER NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN
      Double precision calls,ti,tsi
      COMMON/BVEG4/calls,ti,tsi
      Double precision ALPH
      integer NDMX,MDS
      COMMON/BVEG3/ALPH,NDMX,MDS

C     
C     Global  number of dimensions
C     
      integer Ndimens
      common /to_dimension/ Ndimens
c     
c     external
c     
      double precision fct,fct_mint
      external  fct,fct_mint

      integer          iseed
      common /to_seed/ iseed
c     
c     variables for mint-integrator
c     
      integer ndimmax
      parameter (ndimmax=20)

      double precision xgrid(50,ndimmax),xint,ymax(50,ndimmax)
      double precision mint_ans,mint_err
      integer ncalls0,nitmax,imode

      double precision fun,random,fun_vegas
      external fun,random,fun_vegas

      integer ifold(ndimmax)
      common/cifold/ifold

      double precision store
      common /to_storage/store

      double precision permutation_weight
      common /to_permutation_weight/permutation_weight


C**************************************************************************
C     step 1: initialisation     
C**************************************************************************
      
    
      call global_init

      OPEN(UNIT=89,FILE='./times.out',STATUS='UNKNOWN')
      call cpu_time(time2)
      write(89,*) 't before initialisation',time2
       
      store=0d0
      if(use_perm) then
         call get_num_per(num_per) !calculate the number of permutations.
      else
         num_per=1
      endif

      if (num_per*nb_sol_config.gt.max_channel) then
         write(*,*) "warning: too many channels => increase max_channel"
         stop
      endif

c     VEGAS parameters for the first loop
      ACC = final_prec
      NPRN = 1
      NCALL = nevents
      ITMX = max_it_step1
      ncalls0=ncall  
      nitmax=ITMX  
      imode=0

c     initialization of the variables in the loops
      check_value=0d0
      loop_index=0
      counter=0
      do perm_pos=1,num_per
         do ll=1,nb_sol_config
            counter=counter+1
            order_value(counter)=1d5
         enddo
      enddo
C     
C     initialization of the permutation
C     
 1    normalize_perm=0d0

      loop_index=loop_index+1
      counter=0
      temp_err=0d0
      temp_val=0d0
     
      call cpu_time(time3)
      write(89,*) 't start VEGAS loop',time3

      do perm_pos=1,num_per
         call get_perm(perm_pos, perm_id)
         call assign_perm(perm_id)
         chan_val=0d0
         chan_err=0d0

         permutation_weight=1d0
         call initialize
         normalize_perm =normalize_perm+permutation_weight

         do ll=1,nb_sol_config
            counter=counter+1
            config_pos=ll
            write(*,*) "Current channel of integration: ",config_pos
            write(*,*) "Current parton-jet assignement: ",perm_pos
            write(*,*) "weight of this assignement:     ",permutation_weight

            iseed = iseed + 1   ! avoid to have the same seed
            if (.not. NWA) then
               NDIM=Ndimens
            else
               NDIM=Ndimens-num_propa(config_pos)
            endif
            if(ISR.eq.3) NDIM=NDIM+2
c     
            if(order_value(counter).gt.check_value) then
               if (integrator.eq.'v') then
                  ITMX=4

                  call cpu_time(time1)
                  write(89,*) 't before VEGAS',time1

                  CALL VEGAS(fct,CROSS,SD,CHI)
                  
                  
                  if (loop_index.eq.1) ITMX=max_it_step1 
                  if (loop_index.eq.2) ITMX=max_it_step2
                  CALL VEGAS1(fct,CROSS,SD,CHI)
                  call cpu_time(time)
                  write(89,*) 't end VEGAS1',time
                  write(89,*) 'total VEGAS duration',time-time1
               
               else
                  write(*,*) "problem: unknown integrator "
                  stop
               endif
               if (CROSS.lt.1d99) then
                  temp_val=temp_val+cross
                  chan_val=chan_val+cross
                  temp_err=temp_err+SD
                  chan_err=chan_err+SD
                  order_value(counter) = cross
                  order_error(counter) = sd
                  if (histo) call histo_combine_iteration(counter)
               else
                  order_value(counter)=0d0
               endif
            endif
         enddo
         write(32,*) perm_pos,chan_val, chan_err
      enddo
      
      call cpu_time(time4)
      write(89,*) 't endVEGAS loop',time4

      temp_val=temp_val/dble(normalize_perm)
      temp_err=temp_err/dble(normalize_perm)
      check_value=temp_val/(nb_sol_config) * min_prec_cut1
      if (loop_index.eq.1) then
         NCALL = nevents
         ITMX = max_it_step2
         goto 1
      endif

      write(*,*) "the weight is",temp_val,"+/-",temp_err

      OPEN(UNIT=21,FILE='./weights.out',STATUS='UNKNOWN')
      write(21,*) temp_val,'  ',temp_err
      close(21)

      write(23,*) ' permutation channel   value      error'
      counter=0
      do perm_pos=1,num_per
         call get_perm(perm_pos, perm_id)
         write(23,*) "======================================"
         write(23,*) '1   2', (2+perm_id(i-2), i=3,8)

         do ll=1,nb_sol_config
            counter=counter+1
            write(23,*) perm_pos,'\t',ll,'\t',
     &           order_value(counter),
     &           '\t', order_error(counter)
         enddo
      enddo
      write(23,*) "======================================"
      write(23,*)'Weight: ',temp_val,'+-', temp_err
      write(23,*) "======================================"

      call cpu_time(time)
      write(89,*) 't  end driver',time
      write(89,*) 'duration',time-time2
      close(89)
C**************************************************************************
C     write histogram file (not activated in the standard mode)           *
C**************************************************************************
      if (histo)   call histo_final(num_per*nb_sol_config)  

      close(21)
      OPEN(UNIT=21,FILE='./stop',STATUS='UNKNOWN')
      close(21)

      END


C**************************************************************************************
        subroutine fct(x,wgt,resfct)
        implicit none

        include 'phasespace.inc'
        include 'nexternal.inc'
        include 'run.inc'
        include 'coupl.inc'
c
c       this is the function which is called by the integrator

c
c       parameter
c
        double precision pi
        parameter (pi=3.141592653589793d0)
c
c       arguments
c
        double precision, intent(in) :: x(20),wgt
        double precision, intent(out) :: resfct
c
c       local
c
c        integer i,j ! debug mode
        double precision twgt
c
c       global
c
        double precision              S,X1,X2,PSWGT,JAC
        common /PHASESPACE/ S,X1,X2,PSWGT,JAC
        double precision momenta(0:3,-max_branches:2*max_particles)  ! momenta of external/intermediate legs     (MG order)
        double precision mvir2(-max_branches:2*max_particles)        ! squared invariant masses of intermediate particles (MG order)
        common /to_diagram_kin/ momenta, mvir2


        logical histo
        common /to_histo/histo
c
c       external
c
        double precision dsig
        external dsig
        double precision alphas
        external alphas
        double precision time
        include 'data.inc'
   

c$$$         call get_PS_point(x)


         if (jac.gt.0d0) then
        
c$$$c          here we evaluate the scales if running 
c$$$           if(.not.fixed_ren_scale) then
c$$$             call set_ren_scale(momenta(0,1),scale)
c$$$             if(scale.gt.0) 
c$$$            G = SQRT(4d0*PI*ALPHAS(scale))
c$$$             call UPDATE_AS_PARAM()
c$$$           endif
c$$$           if(.not.fixed_fac_scale) then
c$$$             call set_fac_scale(momenta(0,1),q2fact)
c$$$           endif
c$$$ 
           resfct=jac
c$$$           xbk(1)=X1
c$$$           xbk(2)=X2
c$$$           resfct=resfct*dsig(momenta(0,1),wgt)
c$$$           call cpu_time(time)
c$$$                write(89,*) 'start PS',time
c$$$           call transfer_fct(momenta(0,1),TWGT)
c$$$          call cpu_time(time)
c$$$                write(89,*) 'start PS',time
c$$$           resfct=resfct*twgt
         endif
         end
