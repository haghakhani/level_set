C*******************************************************************
C* Copyright (C) 2003 University at Buffalo
C*
C* This software can be redistributed free of charge.  See COPYING
C* file in the top distribution directory for more details.
C*
C* This software is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C*
C* Author: 
C* Description: 
C*
C*******************************************************************
C* $Id: eigen.f 128 2007-06-07 19:51:52Z dkumar $ 
C*

C***********************************************************************
      subroutine eigen(Uvec, eigenvxmax, eigenvymax, evalue, tiny,
     1     kactxy, gravity, v_solid, v_fluid, eps, flowtype)
C***********************************************************************

      include 'rnr.h'
      integer flowtype
      double precision eigenvxmax, eigenvymax, v_solid(2), v_fluid(2)
      double precision evalue, tiny, Uvec(6), kactxy(2), gravity(3)
      double precision eps
      double precision sound_speed

      if (Uvec(1) .lt.0.d0 ) then
c     iverson and denlinger
         if(kactxy(1) .lt. 0.d0) then
            kactxy(1) = -kactxy(1)
         endif

         if (Uvec(2).gt.0.0001) then
           sound_speed = dsqrt(Uvec(2)*kactxy(1)*gravity(3))
!        x-direction
           eigenvxmax=dabs(v_solid(1)+sound_speed)

!        y-direction
           eigenvymax=dabs(v_solid(2)+sound_speed)

         else
           sound_speed = 0.d0
           eigenvxmax = tiny
           eigenvymax = tiny  
         endif

      else
         eigenvxmax=tiny
         eigenvymax=tiny
      endif
      evalue=dmax1(eigenvxmax,eigenvymax)

c      if ( evalue .lt. 0. ) then
c         write (*,*) "ERROR: Eigen-values -ve"
c         write (*,*) sound_speed, v_solid, v_fluid
c      endif

      return
      end
