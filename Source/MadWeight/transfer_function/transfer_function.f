
C+-----------------------------------------------------------------------+
C|                  TRANSFER FUNCTION FOR MADWEIGHT                      |
C|                                                                       |
C|     Author: Pierre Artoisenet (UCL-CP3)                               |
C|             Olivier Mattelaer (UCL-CP3)                               |
C+-----------------------------------------------------------------------+
C|     This file is generated automaticly by MADWEIGHT-TF_BUILDER        |
C+-----------------------------------------------------------------------+     


C+-----------------------------------------------------------------------+
C|    Transfer function for tf_E_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_E_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'


        prov1=(tf_lepton_E_7+tf_lepton_E_8*dsqrt(p(0))+tf_lepton_E_9*p(0))
        prov2=(tf_lepton_E_10+tf_lepton_E_11*dsqrt(p(0))+tf_lepton_E_12*p(0))

        tf=(exp(-(p(0)-pexp(0)-prov1)**2/2d0/prov2**2))                !first gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2))            !normalisation



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_E_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_E_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'


        width=(tf_lepton_E_10+tf_lepton_E_11*dsqrt(pexp(0))+tf_lepton_E_12*pexp(0))




      width_E_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_THETA_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_lepton
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_lepton(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_lepton
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_lepton(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_PHI_lepton= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_E_jet
C+-----------------------------------------------------------------------+
      subroutine tf_E_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'


        prov1=(tf_jet_E_1+tf_jet_E_2*dsqrt(p(0))+tf_jet_E_3*p(0))
        prov2=(tf_jet_E_4+tf_jet_E_5*dsqrt(p(0))+tf_jet_E_6*p(0))

        tf=(exp(-(p(0)-pexp(0)-prov1)**2/2d0/prov2**2))                !first gaussian
        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2))            !normalisation 	



      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_E_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_E_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'


        width=(tf_jet_E_4+tf_jet_E_5*dsqrt(pexp(0))+tf_jet_E_6*pexp(0))




      width_E_jet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_THETA_jet
C+-----------------------------------------------------------------------+
      subroutine tf_THETA_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_THETA_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_THETA_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_THETA_jet= width

      return
      end



C+-----------------------------------------------------------------------+
C|    Transfer function for tf_PHI_jet
C+-----------------------------------------------------------------------+
      subroutine tf_PHI_jet(pexp,p,n_lhco,weight)
      implicit none

      double precision tf
      double precision pexp(0:3)
      double precision p(0:3)
      integer n_lhco
      double precision weight
      double precision pi
      parameter (pi=3.141592654d0)
      include 'TF_param.inc'

        tf=1d0


      weight=weight*tf

      return
      end

C+-----------------------------------------------------------------------+
C|    Definition of the WIDTH associated to tf_PHI_jet
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION FUNCTION width_PHI_jet(pexp,n_lhco)
      implicit none

       	  double precision width
      double precision pexp(0:3)
      integer n_lhco

      double precision pi
      parameter (pi=3.141592654d0)

      include 'TF_param.inc'

        width=0d0


      width_PHI_jet= width

      return
      end



