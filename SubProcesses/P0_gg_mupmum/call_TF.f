
C+-----------------------------------------------------------------------+
C|                  TRANSFER FUNCTION FOR MADWEIGHT                      |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-TF_BUILDER        |
C+-----------------------------------------------------------------------+     


C+-----------------------------------------------------------------------+
C|     Subroutine: GET_CENTRAL_POINT                                     |
C|                                                                       |
C|     purpose: Define the central point of integration c_point(I,J,1)   |
C|              and the width associated c_point(I,J,2)                  |
C|                 I: MadGraph Position                                  |
C|                 J: 1,2,3 for THETA,PHI,RHO variable                   |
C+-----------------------------------------------------------------------+
      subroutine get_central_point
      implicit none
c                                                                        
c     parameter                                                          
c                                                                        
      include 'run.inc'
      include 'nexternal.inc'
      include 'TF_param.inc'

c                                                                        
c     global                                                             
c                                                                        
      double precision c_point(1:max_particles,3,2)
      common/ph_sp_init/c_point
c                                                                        
      double precision p_exp(0:3,nexternal)
      common /to_pexp/p_exp
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c                                                                        
c     local                                                              
c                                                                       
      integer i
      integer n_lhco


        external width_E_lepton, width_THETA_lepton, width_PHI_lepton
        double precision width_E_lepton, width_THETA_lepton, width_PHI_lepton

C+-----------------------------------------------------------------------+
C|     Start the definition                                              |
C+-----------------------------------------------------------------------+


        c_point(1,1,2)=-1d0
        c_point(1,2,2)=-1d0
        c_point(1,3,2)=-1d0

        c_point(2,1,2)=-1d0
        c_point(2,2,2)=-1d0
        c_point(2,3,2)=-1d0

        c_point(3,1,1)=theta(p_exp(0,3))
        c_point(3,2,1)=phi(p_exp(0,3))
        c_point(3,3,1)=rho(p_exp(0,3))
        c_point(3,1,2)=width_THETA_lepton(p_exp(0,3),tag_lhco(3))
        c_point(3,2,2)=width_PHI_lepton(p_exp(0,3),tag_lhco(3))
        c_point(3,3,2)=width_E_lepton(p_exp(0,3),tag_lhco(3))

        c_point(4,1,1)=theta(p_exp(0,4))
        c_point(4,2,1)=phi(p_exp(0,4))
        c_point(4,3,1)=rho(p_exp(0,4))
        c_point(4,1,2)=width_THETA_lepton(p_exp(0,4),tag_lhco(4))
        c_point(4,2,2)=width_PHI_lepton(p_exp(0,4),tag_lhco(4))
        c_point(4,3,2)=width_E_lepton(p_exp(0,4),tag_lhco(4))


        return
        end

C+-----------------------------------------------------------------------+
C|     Subroutine: transfer_fct(P,weight)                                |
C|                                                                       |
C|     purpose: returns the weight by the coeficient comming from the    |
C|              transfer_fct                                             |
C+-----------------------------------------------------------------------+
            subroutine transfer_fct(P,weight)
      implicit none
c                                                 
      integer    maxexternal
      parameter (maxexternal=15)
c                               
c     ARGUMENTS                 
c                               
      DOUBLE PRECISION P(0:3,maxexternal)
      DOUBLE PRECISION weight
c                                        
c     INCLUDE and COMMON                 
c                                        
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c                            
c                            
      double precision pexp(0:3,nexternal)
      integer tag_lhco(3:nexternal)
      common/to_pexp/pexp
      common/lhco_order/tag_lhco
C                               
C     SPECIAL FCT              
C                               
      DOUBLE PRECISION R2,DOT,ET,ETA,DJ,SumDot,PT
c                
c     LOCAL                                      
c                
      integer i,k
      integer n_lhco


        weight=1d0
        n_lhco=tag_lhco(3)
        call tf_E_lepton(pexp(0,3),p(0,3),n_lhco,weight)
        call tf_THETA_lepton(pexp(0,3),p(0,3),n_lhco,weight)
        call tf_PHI_lepton(pexp(0,3),p(0,3),n_lhco,weight)

        n_lhco=tag_lhco(4)
        call tf_E_lepton(pexp(0,4),p(0,4),n_lhco,weight)
        call tf_THETA_lepton(pexp(0,4),p(0,4),n_lhco,weight)
        call tf_PHI_lepton(pexp(0,4),p(0,4),n_lhco,weight)


        call check_nan(weight)
        return
        end

C+-----------------------------------------------------------------------+
C|     Subroutine: tf_E_for_part(MG_num)                                 |
C|                                                                       |
C|     purpose: returns the value of the transfert function (in energy)  |
C|              for the particles associated to the MG_number            |
C+-----------------------------------------------------------------------+
      double precision function  tf_E_for_part(MG_num)
      implicit none
c
c     ARGUMENTS
c
      INTEGER MG_num
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer n_lhco
c


        if(MG_num.eq.1) then
        tf_E_for_part=1d0
        return
        endif
        if(MG_num.eq.2) then
        tf_E_for_part=1d0
        return
        endif
        if(MG_num.eq.3) then
        tf_E_for_part=1d0
        n_lhco=tag_lhco(3)
        call tf_E_lepton(pexp(0,3),momenta(0,3),n_lhco,tf_E_for_part)

        return
        endif
        if(MG_num.eq.4) then
        tf_E_for_part=1d0
        n_lhco=tag_lhco(4)
        call tf_E_lepton(pexp(0,4),momenta(0,4),n_lhco,tf_E_for_part)

        return
        endif

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_E_for_XX()                                    |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_E_for_1()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_E_for_1=1d0

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_E_for_XX()                                    |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_E_for_2()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_E_for_2=1d0

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_E_for_XX()                                    |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+




        double precision function tf_E_for_3()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_E_for_3=1d0
        n_lhco=tag_lhco(3)
        call tf_E_lepton(pexp(0,3),momenta(0,3),n_lhco,tf_E_for_3)

        return
        end


C+-----------------------------------------------------------------------+
C|          Subroutine: tf_E_for_XX()                                    |
C|                                                                       |
C|         purpose: returns the value of the transfer function (in energy)
C|                                                                       |
C+-----------------------------------------------------------------------+



        double precision function tf_E_for_4()

      implicit none
c
c     INCLUDE and COMMON
c
      include 'phasespace.inc'
      include 'nexternal.inc'
      include 'run.inc'
c
      double precision pexp(0:3,nexternal)
      common/to_pexp/pexp
      double precision momenta(0:3,-max_branches:2*max_particles)
      double precision mvir2(-max_branches:2*max_particles)
      common/to_diagram_kin/momenta,mvir2
      integer tag_lhco(3:nexternal)
      common/lhco_order/tag_lhco
c
c     LOCAL
c
      integer i,k,n_lhco


        tf_E_for_4=1d0
        n_lhco=tag_lhco(4)
        call tf_E_lepton(pexp(0,4),momenta(0,4),n_lhco,tf_E_for_4)

        return
        end



