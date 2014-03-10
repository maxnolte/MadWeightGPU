      module dhelascuda
      contains
        attributes(device) subroutine ixxxxx_cuda(p, fmass, nhel, nsf ,fi)

      implicit none
      double complex fi(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,ip,im,nh

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )


      fi(1) = dcmplx(p(0),p(3))*nsf*-1
      fi(2) = dcmplx(p(1),p(2))*nsf*-1

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then

            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses sdf
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            ip = (1+nh)/2
            im = (1-nh)/2

            fi(3) = ip     * sqm(ip)
            fi(4) = im*nsf * sqm(ip)
            fi(5) = ip*nsf * sqm(im)
            fi(6) = im     * sqm(im)

         else

            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , p(2) )/dsqrt(rTwo*pp*pp3)
            endif

            fi(3) = sfomeg(1)*chi(im)
            fi(4) = sfomeg(1)*chi(ip)
            fi(5) = sfomeg(2)*chi(im)
            fi(6) = sfomeg(2)*chi(ip)

         endif

      else

         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fi(3) = dcmplx( rZero )
            fi(4) = dcmplx( rZero )
            fi(5) = chi(1)
            fi(6) = chi(2)
         else
            fi(3) = chi(2)
            fi(4) = chi(1)
            fi(5) = dcmplx( rZero )
            fi(6) = dcmplx( rZero )
         endif
      endif
c
      return
      end subroutine ixxxxx_cuda


       attributes(device) subroutine vxxxxx_cuda(p,vmass,nhel,nsv , vc)
c
c This subroutine computes a VECTOR wavefunction.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer nhel = -1, 0, 1: helicity      of vector boson
c                                (0 is forbidden if vmass=0.0)
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c
      implicit none
      double complex vc(6)
      double precision p(0:3),vmass,hel,hel0,pt,pt2,pp,pzpt,emp,sqh
      integer nhel,nsv,nsvahl

      double precision rZero, rHalf, rOne, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0 )


      sqh = dsqrt(rHalf)
      hel = dble(nhel)
      nsvahl = nsv*dabs(hel)
      pt2 = p(1)**2+p(2)**2
      pp = min(p(0),dsqrt(pt2+p(3)**2))
      pt = min(pp,dsqrt(pt2))

      vc(1) = dcmplx(p(0),p(3))*nsv
      vc(2) = dcmplx(p(1),p(2))*nsv


      if ( vmass.ne.rZero ) then

         hel0 = rOne-dabs(hel)

         if ( pp.eq.rZero ) then

            vc(3) = dcmplx( rZero )
            vc(4) = dcmplx(-hel*sqh )
            vc(5) = dcmplx( rZero , nsvahl*sqh )
            vc(6) = dcmplx( hel0 )

         else

            emp = p(0)/(vmass*pp)
            vc(3) = dcmplx( hel0*pp/vmass )
            vc(6) = dcmplx( hel0*p(3)*emp+hel*pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = p(3)/(pp*pt)*sqh*hel
               vc(4) = dcmplx( hel0*p(1)*emp-p(1)*pzpt ,
     &                         -nsvahl*p(2)/pt*sqh       )
               vc(5) = dcmplx( hel0*p(2)*emp-p(2)*pzpt ,
     &                          nsvahl*p(1)/pt*sqh       )
            else
               vc(4) = dcmplx( -hel*sqh )
               vc(5) = dcmplx( rZero , nsvahl*sign(sqh,p(3)) )
            endif

         endif

      else

         pp = p(0)
         pt = sqrt(p(1)**2+p(2)**2)
         vc(3) = dcmplx( rZero )
         vc(6) = dcmplx( hel*pt/pp*sqh )
         if ( pt.ne.rZero ) then
            pzpt = p(3)/(pp*pt)*sqh*hel
            vc(4) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
            vc(5) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
         else
            vc(4) = dcmplx( -hel*sqh )
            vc(5) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
         endif

      endif
c
      return
      end subroutine vxxxxx_cuda

       attributes(device)    subroutine oxxxxx_cuda(p,fmass,nhel,nsf , fo)
c
c This subroutine computes a fermion wavefunction with the flowing-OUT
c fermion number.
c
c input:
c       real    p(0:3)         : four-momentum of fermion
c       real    fmass          : mass          of fermion
c       integer nhel = -1 or 1 : helicity      of fermion
c       integer nsf  = -1 or 1 : +1 for particle, -1 for anti-particle
c
c output:
c       complex fo(6)          : fermion wavefunction               <fo|
c
      implicit none
      double complex fo(6),chi(2)
      double precision p(0:3),sf(2),sfomeg(2),omega(2),fmass,
     &     pp,pp3,sqp0p3,sqm(0:1)
      integer nhel,nsf,nh,ip,im

      double precision rZero, rHalf, rTwo
      parameter( rZero = 0.0d0, rHalf = 0.5d0, rTwo = 2.0d0 )

      fo(1) = dcmplx(p(0),p(3))*nsf
      fo(2) = dcmplx(p(1),p(2))*nsf

      nh = nhel*nsf

      if ( fmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then

            sqm(0) = dsqrt(abs(fmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),fmass) ! possibility of negative fermion masses
            ip = -((1+nh)/2)
            im =  (1-nh)/2

            fo(3) = im     * sqm(im)
            fo(4) = ip*nsf * sqm(im)
            fo(5) = im*nsf * sqm(-ip)
            fo(6) = ip     * sqm(-ip)

         else

            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsf+(1-nsf)*nh)*rHalf
            sf(2) = dble(1+nsf-(1-nsf)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = fmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif

            fo(3) = sfomeg(2)*chi(im)
            fo(4) = sfomeg(2)*chi(ip)
            fo(5) = sfomeg(1)*chi(im)
            fo(6) = sfomeg(1)*chi(ip)

         endif

      else

         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsf
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fo(3) = chi(1)
            fo(4) = chi(2)
            fo(5) = dcmplx( rZero )
            fo(6) = dcmplx( rZero )
         else
            fo(3) = dcmplx( rZero )
            fo(4) = dcmplx( rZero )
            fo(5) = chi(2)
            fo(6) = chi(1)
         endif

      endif
c
      return
      end subroutine oxxxxx_cuda

      attributes(device) SUBROUTINE FFS1_0_cuda(F1, F2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 S3(*)
      COMPLEX*16 TMP0
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      TMP0 = (F2(3)*F1(3)+F2(4)*F1(4)+F2(5)*F1(5)+F2(6)*F1(6))
      VERTEX = COUP*-CI * TMP0*S3(3)
      END subroutine FFS1_0_cuda

      attributes(device) SUBROUTINE FFS1_3_cuda(F1, F2, COUP, M3, W3,S3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 DENOM
      COMPLEX*16 S3(3)
      REAL*8 W3
      COMPLEX*16 TMP0
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 COUP
      S3(1) = +F1(1)+F2(1)
      S3(2) = +F1(2)+F2(2)
      P3(0) = -DBLE(S3(1))
      P3(1) = -DBLE(S3(2))
      P3(2) = -DIMAG(S3(2))
      P3(3) = -DIMAG(S3(1))
      TMP0 = (F2(3)*F1(3)+F2(4)*F1(4)+F2(5)*F1(5)+F2(6)*F1(6))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 
     $ -CI* W3))
      S3(3)= DENOM*CI * TMP0
      END SUBROUTINE FFS1_3_cuda

      attributes(device) SUBROUTINE VVS3_3_cuda(V1, V2, COUP, M3, W3,S3)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 TMP2
      COMPLEX*16 S3(3)
      COMPLEX*16 TMP1
      REAL*8 P1(0:3)
      REAL*8 W3
      REAL*8 P2(0:3)
      REAL*8 P3(0:3)
      REAL*8 M3
      COMPLEX*16 DENOM
      COMPLEX*16 TMP4
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP3
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      S3(1) = +V1(1)+V2(1)
      S3(2) = +V1(2)+V2(2)
      P3(0) = -DBLE(S3(1))
      P3(1) = -DBLE(S3(2))
      P3(2) = -DIMAG(S3(2))
      P3(3) = -DIMAG(S3(1))
      TMP4 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP1 = (V2(3)*P1(0)-V2(4)*P1(1)-V2(5)*P1(2)-V2(6)*P1(3))
      TMP3 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP2 = (V1(3)*P2(0)-V1(4)*P2(1)-V1(5)*P2(2)-V1(6)*P2(3))
      DENOM = COUP/(P3(0)**2-P3(1)**2-P3(2)**2-P3(3)**2 - M3 * (M3 
     $ -CI* W3))
      S3(3)= DENOM*(-CI*(TMP3*TMP4)+CI*(TMP1*TMP2))
      END SUBROUTINE VVS3_3_cuda

      end module
