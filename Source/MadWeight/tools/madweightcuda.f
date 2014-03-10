c module for function

      module madweightcuda
      contains
        attributes(global) subroutine fct_cuda(xarray,wgtarray,resfctarray, npoints)
        implicit none

c$$$        include 'phasespace.inc'
c$$$        include 'nexternal.inc'
c$$$        include 'run.inc'
c$$$        include 'coupl.inc'

        double precision pi
        parameter (pi=3.141592653589793d0)

c       Arguments
c
        double precision, dimension(:,:) :: xarray
        double precision, dimension(:) :: wgtarray,resfctarray
        integer, value :: npoints
        integer :: i
c$$$c
c$$$c       local
c$$$c
c$$$c        integer i,j ! debug mode
c$$$        double precision twgt
c$$$c
c$$$c       global
c$$$c
c$$$        double precision              S,X1,X2,PSWGT,JAC
c$$$c       common /PHASESPACE/ S,X1,X2,PSWGT,JAC
c$$$c       double precision momenta(0:3,-max_branches:2*max_particles)  ! momenta of external/intermediate legs     (MG order)
c$$$c       double precision mvir2(-max_branches:2*max_particles)        ! squared invariant masses of intermediate particles (MG order)
c$$$c       common /to_diagram_kin/ momenta, mvir2
c$$$
c$$$
c$$$c       logical histo
c$$$c       common /to_histo/histo
c$$$c
c$$$c       external
c$$$c
c$$$        double precision dsig
c$$$c       external dsig
c$$$        double precision alphas
c$$$c       external alphas
c$$$        double precision time
c$$$c       include 'data.inc'

c      CUDA number of thread
         i = threadidx%x + (blockidx%x-1)*blockdim%x
         if ( i <= npoints) resfctarray(i) = 2.d0
c$$$    
c$$$c$$$         call get_PS_point(x)
c$$$
c$$$         resfctarray(i)=1.d0
c$$$         if (jac.gt.0d0) then
c$$$        
c$$$c$$$c          here we evaluate the scales if running 
c$$$c$$$           if(.not.fixed_ren_scale) then
c$$$c$$$             call set_ren_scale(momenta(0,1),scale)
c$$$c$$$             if(scale.gt.0) 
c$$$c$$$            G = SQRT(4d0*PI*ALPHAS(scale))
c$$$c$$$             call UPDATE_AS_PARAM()
c$$$c$$$           endif
c$$$c$$$           if(.not.fixed_fac_scale) then
c$$$c$$$             call set_fac_scale(momenta(0,1),q2fact)
c$$$c$$$           endif
c$$$c$$$ 
c$$$           resfctarray(i)=jac
c$$$c$$$           xbk(1)=X1
c$$$c$$$           xbk(2)=X2
c$$$c$$$           resfct=resfct*dsig(momenta(0,1),wgt)
c$$$c$$$           call cpu_time(time)
c$$$c$$$                write(89,*) 'start PS',time
c$$$c$$$           call transfer_fct(momenta(0,1),TWGT)
c$$$c$$$          call cpu_time(time)
c$$$c$$$                write(89,*) 'start PS',time
c$$$c$$$           resfct=resfct*twgt
c$$$           resfctarray(i)=2.d0
c$$$         endif
c$$$         endif
         end subroutine fct_cuda
      end module madweightcuda

  
