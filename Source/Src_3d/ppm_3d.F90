module ppm_module

  implicit none

  private

  public ppm

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                 u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                 Ip,Im, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc, &
                 force_type_in)

    use meth_params_module, only : ppm_type

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision     s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision     u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision  cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)
    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision dx,dy,dz,dt
    integer          k3d,kc
   
    integer, intent(in), optional :: force_type_in

    integer :: ppm_type_to_use

    ppm_type_to_use = ppm_type
    if (present(force_type_in)) ppm_type_to_use = force_type_in

    if (ppm_type_to_use == 1) then

        call ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    else if (ppm_type_to_use == 2) then

        call ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine ppm_type1(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use bl_constants_module
  
    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision     s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision     u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision  cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn( f_l1: f_h1, f_l2: f_h2, f_l3: f_h3)

    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision dx,dy,dz,dt
    integer          k3d,kc

    ! local
    integer i,j

    double precision dsl, dsr, dsc, dtdx, dtdy, dtdz
    double precision sigma, s6, dsvlm, dsvlp, dsvlz, ssp, ssm

    ! s_{\ib,+}, s_{\ib,-}
    double precision :: sp(ilo1-1:ihi1+1)
    double precision :: sm(ilo1-1:ihi1+1)

    ! \delta s_{\ib}^{vL}
    double precision :: dsvlx(ilo1-2:ihi1+2)
    double precision, allocatable :: dsvly(:,:)

    ! s_{i+\half}^{H.O.}
    double precision :: sedgex(ilo1-1:ihi1+2)
    double precision, allocatable :: sedgey(:,:)

    allocate(dsvly (ilo1-1:ihi1+1, ilo2-2:ihi2+2))
    allocate(sedgey(ilo1-1:ihi1+1, ilo2-1:ihi2+2))

    if (ppm_type .ne. 1) &
         call bl_error("Should have ppm_type = 1 in ppm_type1")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call bl_error("Need more ghost cells on array in ppm_type1")
    end if

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at x-edges

    do j=ilo2-1,ihi2+1

       ! compute van Leer slopes in x-direction

       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(dsc,dsl,dsr)
       do i=ilo1-2,ihi1+2
          dsc = FOURTH * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl =          (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr =          (s(i+1,j,k3d) - s(i  ,j,k3d))
          dsvlx(i) = (ONE+sign(ONE,dsl*dsr))*sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do

       ! interpolate s to x-edges
       !DIR$ SIMD vectorlength(BL_SIMD_LEN)
       do i=ilo1-1,ihi1+2
          sedgex(i) = HALF*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - SIXTH*(dsvlx(i)-dsvlx(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedgex(i) = max(sedgex(i),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedgex(i) = min(sedgex(i),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do

       ! flatten the parabola BEFORE doing the other                     
       ! monotonization -- this is the method that Flash does       
       if (ppm_flatten_before_integrals == 1) then
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+2
             sedgex(i) = flatn(i,j,k3d)*sedgex(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       end if

       ! modify using quadratic limiters -- note this version of the limiting comes
       ! from Colella and Sekora (2008), not the original PPM paper.
       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(ssm,ssp)
       do i=ilo1-1,ihi1+1
          ssm = s(i,j,k3d) - sedgex(i)
          ssp = sedgex(i+1) - s(i,j,k3d)

          sp(i) = sedgex(i+1)
          sm(i) = sedgex(i)

          if (ssm*ssp .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(ssp) .ge. TWO*abs(ssm)) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(ssm) .ge. TWO*abs(ssp)) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       ! flatten the parabola AFTER doing the monotonization --
       ! this is the method that Miller & Colella do
       if (ppm_flatten_before_integrals == 2) then
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+1
             sm(i) = flatn(i,j,k3d)*sm(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i) = flatn(i,j,k3d)*sp(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       endif

       ! compute x-component of Ip and Im

       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(s6,sigma)
       do i=ilo1-1,ihi1+1
          s6 = SIX*s(i,j,k3d) - THREE*(sm(i)+sp(i))

          ! Ip/m is the integral under the parabola for the extent
          ! that a wave can travel over a timestep
          !
          ! Ip integrates to the right edge of a cell
          ! Im integrates to the left edge of a cell

          ! u-c wave
          sigma = (u(i,j,k3d,1) - cspd(i,j,k3d))*dtdx
          if (sigma < ZERO) then
             Ip(i,j,kc,1,1) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,1,1) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,1,1) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,1,1) = sm(i)              
          end if

          ! u wave             
          sigma = u(i,j,k3d,1)*dtdx
          if (sigma < ZERO) then
             Ip(i,j,kc,1,2) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,1,2) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,1,2) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,1,2) = sm(i)              
          end if

          ! u+c wave
          sigma = (u(i,j,k3d,1) + cspd(i,j,k3d))*dtdx
          if (sigma < ZERO) then
             Ip(i,j,kc,1,3) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,1,3) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,1,3) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,1,3) = sm(i) 
          end if
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    do j=ilo2-2,ihi2+2
       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(dsc,dsl,dsr)
       do i=ilo1-1,ihi1+1
          dsc = FOURTH * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl =          (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr =          (s(i,j+1,k3d) - s(i,j  ,k3d))
          dsvly(i,j) = (ONE+sign(ONE,dsl*dsr))*sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       !DIR$ SIMD vectorlength(BL_SIMD_LEN)
       do i=ilo1-1,ihi1+1
          sedgey(i,j) = HALF*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - SIXTH*(dsvly(i,j)-dsvly(i,j-1))
          ! make sure sedgey lies in between adjacent cell-centered values
          sedgey(i,j) = max(sedgey(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedgey(i,j) = min(sedgey(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    ! flatten the parabola BEFORE doing the other                     
    ! monotonization -- this is the method that Flash does       
    if (ppm_flatten_before_integrals == 1) then
       do j=ilo2-1,ihi2+2
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+1
             sedgey(i,j) = flatn(i,j,k3d)*sedgey(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       end do
    end if
    
    ! modify using quadratic limiters
    do j=ilo2-1,ihi2+1
       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(ssm,ssp)
       do i=ilo1-1,ihi1+1
          ssm = s(i,j,k3d) - sedgey(i,j)
          ssp = sedgey(i,j+1) - s(i,j,k3d)

          sp(i) = sedgey(i,j+1)
          sm(i) = sedgey(i,j)

          if (ssm*ssp .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(ssp) .ge. TWO*abs(ssm)) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(ssm) .ge. TWO*abs(ssp)) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       ! flatten the parabola AFTER doing the monotonization --
       ! this is the method that Miller & Colella do
       if (ppm_flatten_before_integrals == 2) then
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+1
             sm(i) = flatn(i,j,k3d)*sm(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i) = flatn(i,j,k3d)*sp(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       endif

       ! compute y-component of Ip and Im

       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(s6,sigma)
       do i=ilo1-1,ihi1+1
          s6 = SIX*s(i,j,k3d) - THREE*(sm(i)+sp(i))

          ! v-c wave
          sigma = (u(i,j,k3d,2) - cspd(i,j,k3d))*dtdy
          if (sigma < ZERO) then
             Ip(i,j,kc,2,1) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,2,1) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,2,1) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,2,1) = sm(i) 
          end if

          ! v wave          
          sigma =u(i,j,k3d,2)*dtdy
          if (sigma < ZERO) then
             Ip(i,j,kc,2,2) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,2,2) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,2,2) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,2,2) = sm(i)              
          end if

          ! v+c wave
          sigma = (u(i,j,k3d,2) + cspd(i,j,k3d))*dtdy
          if (sigma < ZERO) then
             Ip(i,j,kc,2,3) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,2,3) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,2,3) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,2,3) = sm(i) 
          end if
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j=ilo2-1,ihi2+1

       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(dsc,dsl,dsr,dsvlm,dsvlp,dsvlz)
       do i=ilo1-1,ihi1+1
          ! compute on slab below
          dsc = FOURTH * (s(i,j,k3d  ) - s(i,j,k3d-2))
          dsl =          (s(i,j,k3d-1) - s(i,j,k3d-2))
          dsr =          (s(i,j,k3d  ) - s(i,j,k3d-1))
          dsvlm = (ONE+sign(ONE,dsl*dsr))*sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on current slab
          dsc = FOURTH * (s(i,j,k3d+1) - s(i,j,k3d-1))
          dsl =          (s(i,j,k3d  ) - s(i,j,k3d-1))
          dsr =          (s(i,j,k3d+1) - s(i,j,k3d  ))
          dsvlz = (ONE+sign(ONE,dsl*dsr))*sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on slab above
          dsc = FOURTH * (s(i,j,k3d+2) - s(i,j,k3d  ))
          dsl =          (s(i,j,k3d+1) - s(i,j,k3d  ))
          dsr =          (s(i,j,k3d+2) - s(i,j,k3d+1))
          dsvlp = (ONE+sign(ONE,dsl*dsr))*sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! interpolate to lo face
          sm(i) = HALF*(s(i,j,k3d)+s(i,j,k3d-1)) - SIXTH*(dsvlz-dsvlm)
          ! make sure sedge lies in between adjacent cell-centered values
          sm(i) = max(sm(i),min(s(i,j,k3d),s(i,j,k3d-1)))
          sm(i) = min(sm(i),max(s(i,j,k3d),s(i,j,k3d-1)))

          ! interpolate to hi face
          sp(i) = HALF*(s(i,j,k3d)+s(i,j,k3d+1)) - SIXTH*(dsvlp-dsvlz)
          ! make sure sedge lies in between adjacent cell-centered values
          sp(i) = max(sp(i),min(s(i,j,k3d),s(i,j,k3d+1)))
          sp(i) = min(sp(i),max(s(i,j,k3d),s(i,j,k3d+1)))
       end do

       ! flatten the parabola BEFORE doing the other                     
       ! monotonization -- this is the method that Flash does       
       if (ppm_flatten_before_integrals == 1) then
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+1
             sm(i) = flatn(i,j,k3d)*sm(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i) = flatn(i,j,k3d)*sp(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       endif

       ! modify using quadratic limiters
       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(ssm,ssp)
       do i=ilo1-1,ihi1+1
          ssm = s(i,j,k3d) - sm(i)
          ssp = sp(i) - s(i,j,k3d)

          if (ssm*ssp .le. ZERO) then
             sp(i) = s(i,j,k3d)
             sm(i) = s(i,j,k3d)
          else if (abs(ssp) .ge. TWO*abs(ssm)) then
             sp(i) = THREE*s(i,j,k3d) - TWO*sm(i)
          else if (abs(ssm) .ge. TWO*abs(ssp)) then
             sm(i) = THREE*s(i,j,k3d) - TWO*sp(i)
          end if
       end do

       ! flatten the parabola AFTER doing the monotonization --
       ! this is the method that Miller & Colella do
       if (ppm_flatten_before_integrals == 2) then
          !DIR$ SIMD vectorlength(BL_SIMD_LEN)
          do i=ilo1-1,ihi1+1
             sm(i) = flatn(i,j,k3d)*sm(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i) = flatn(i,j,k3d)*sp(i) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          end do
       endif

       ! compute z-component of Ip and Im
       !DIR$ SIMD vectorlength(BL_SIMD_LEN) &
       !DIR$ private(s6,sigma)
       do i=ilo1-1,ihi1+1
          s6 = SIX*s(i,j,k3d) - THREE*(sm(i)+sp(i))

          ! w-c wave
          sigma = (u(i,j,k3d,3) - cspd(i,j,k3d))*dtdz
          if (sigma < ZERO) then
             Ip(i,j,kc,3,1) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,3,1) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,3,1) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,3,1) = sm(i) 
          end if

          ! w wave          
          sigma = u(i,j,k3d,3)*dtdz
          if (sigma < ZERO) then
             Ip(i,j,kc,3,2) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,3,2) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,3,2) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,3,2) = sm(i) 
          end if

          ! w+c wave
          sigma = (u(i,j,k3d,3) + cspd(i,j,k3d))*dtdz
          if (sigma < ZERO) then
             Ip(i,j,kc,3,3) = sp(i) + HALF*sigma*(sp(i)-sm(i)-(ONE+TWO3RD*sigma)*s6)
          else
             Ip(i,j,kc,3,3) = sp(i)
          end if
          if (sigma > ZERO) then
             Im(i,j,kc,3,3) = sm(i) + HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          else
             Im(i,j,kc,3,3) = sm(i) 
          end if
       end do
    end do

    deallocate(dsvly,sedgey)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                       u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       flatn,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                       Ip,Im,ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)

    use meth_params_module, only : ppm_type, ppm_flatten_before_integrals
    use bl_constants_module

    implicit none

    integer           s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer           f_l1, f_l2, f_l3, f_h1, f_h2, f_h3
    integer          ilo1,ilo2,ihi1,ihi2
    double precision s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    double precision u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,1:3)
    double precision cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)
    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision dx,dy,dz,dt
    integer          k3d,kc

    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    double precision D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision sgn, sigma, s6
    double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    double precision dachkm, dachkp
    double precision amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:,:)
    double precision, allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:,:)
    double precision, allocatable :: dsvlm(:,:)
    double precision, allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 2) &
         call bl_error("Should have ppm_type = 2 in ppm_type2")

    if (s_l1 .gt. ilo1-3 .or. s_l2 .gt. ilo2-3) then
         print *,'Low bounds of array: ',s_l1, s_l2
         print *,'Low bounds of  loop: ',ilo1 , ilo2
         call bl_error("Need more ghost cells on array in ppm_type2")
    end if

    if (s_h1 .lt. ihi1+3 .or. s_h2 .lt. ihi2+3) then
         print *,'Hi  bounds of array: ',s_h1, s_h2
         print *,'Hi  bounds of  loop: ',ihi1 , ihi2
         call bl_error("Need more ghost cells on array in ppm_type2")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+3
          sedge(i,j) = SEVEN12TH*(s(i-1,j,k3d)+s(i  ,j,k3d)) &
               - TWELFTH*(s(i-2,j,k3d)+s(i+1,j,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i-1,j,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. ZERO) then
             D2  = THREE*(s(i-1,j,k3d)-TWO*sedge(i,j)+s(i,j,k3d))
             D2L = s(i-2,j,k3d)-TWO*s(i-1,j,k3d)+s(i,j,k3d)
             D2R = s(i-1,j,k3d)-TWO*s(i,j,k3d)+s(i+1,j,k3d)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i,j) = HALF*(s(i-1,j,k3d)+s(i,j,k3d)) - SIXTH*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i+1,j)-s(i,j,k3d)
          alpham   = sedge(i  ,j)-s(i,j,k3d)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i-1,j)
             dafacep   = sedge(i+2,j) - sedge(i+1,j)
             dabarm    = s(i,j,k3d) - s(i-1,j,k3d)
             dabarp    = s(i+1,j,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i-2,j,k3d)-TWO*s(i-1,j,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-TWO*s(i+1,j,k3d)+s(i+2,j,k3d)
             D2C    = s(i-1,j,k3d)-TWO*s(i,j,k3d)+s(i+1,j,k3d)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization 
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute x-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))

          ! u-c wave
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt/dx

          if (u(i,j,k3d,1)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,1) = sp(i,j)
          else
             Ip(i,j,kc,1,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,1) = sm(i,j)
          else
             Im(i,j,kc,1,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u wave
          sigma = abs(u(i,j,k3d,1))*dt/dx

          if (u(i,j,k3d,1) <= ZERO) then
             Ip(i,j,kc,1,2) = sp(i,j)
          else
             Ip(i,j,kc,1,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1) >= ZERO) then
             Im(i,j,kc,1,2) = sm(i,j)
          else
             Im(i,j,kc,1,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u+c wave
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt/dx

          if (u(i,j,k3d,1)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,1,3) = sp(i,j) 
          else
             Ip(i,j,kc,1,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,1)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,1,3) = sm(i,j) 
          else
             Im(i,j,kc,1,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

       end do
    end do

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))

    ! compute s at y-edges

    ! interpolate s to y-edges
    do j=ilo2-2,ihi2+3
       do i=ilo1-1,ihi1+1
          sedge(i,j) = SEVEN12TH*(s(i,j-1,k3d)+s(i,j,k3d)) &
               - TWELFTH*(s(i,j-2,k3d)+s(i,j+1,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i,j-1,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. ZERO) then
             D2  = THREE*(s(i,j-1,k3d)-TWO*sedge(i,j)+s(i,j,k3d))
             D2L = s(i,j-2,k3d)-TWO*s(i,j-1,k3d)+s(i,j,k3d)
             D2R = s(i,j-1,k3d)-TWO*s(i,j,k3d)+s(i,j+1,k3d)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i,j) = HALF*(s(i,j-1,k3d)+s(i,j,k3d)) - SIXTH*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i,j+1)-s(i,j,k3d)
          alpham   = sedge(i,j  )-s(i,j,k3d)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i,j-1)
             dafacep   = sedge(i,j+2) - sedge(i,j+1)
             dabarm    = s(i,j,k3d) - s(i,j-1,k3d)
             dabarp    = s(i,j+1,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i,j-2,k3d)-TWO*s(i,j-1,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-TWO*s(i,j+1,k3d)+s(i,j+2,k3d)
             D2C    = s(i,j-1,k3d)-TWO*s(i,j,k3d)+s(i,j+1,k3d)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j-1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j+1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization 
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif


          !
          ! Compute y-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))

          ! v-c wave
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt/dy

          if (u(i,j,k3d,2)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,1) = sp(i,j) 
          else
             Ip(i,j,kc,2,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,1) = sm(i,j) 
          else
             Im(i,j,kc,2,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! v wave
          sigma = abs(u(i,j,k3d,2))*dt/dy

          if (u(i,j,k3d,2) <= ZERO) then
             Ip(i,j,kc,2,2) = sp(i,j) 
          else
             Ip(i,j,kc,2,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2) >= ZERO) then
             Im(i,j,kc,2,2) = sm(i,j) 
          else
             Im(i,j,kc,2,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! v+c wave
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt/dy

          if (u(i,j,k3d,2)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,2,3) = sp(i,j) 
          else
             Ip(i,j,kc,2,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,2)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,2,3) = sm(i,j) 
          else
             Im(i,j,kc,2,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

       end do
    end do

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! interpolate s to z-edges
    do k=k3d-1,k3d+2
       do j=ilo2-1,ihi2+1
          do i=ilo1-1,ihi1+1
             sedgez(i,j,k) = SEVEN12TH*(s(i,j,k-1)+s(i,j,k)) &
                  - TWELFTH*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i,j,k-1)-TWO*sedgez(i,j,k)+s(i,j,k))
                D2L = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
                D2R = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
                sgn = sign(ONE,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
                sedgez(i,j,k) = HALF*(s(i,j,k-1)+s(i,j,k)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    k = k3d
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedgez(i,j,k+1)-s(i,j,k)
          alpham   = sedgez(i,j,k  )-s(i,j,k)
          bigp     = abs(alphap).gt.TWO*abs(alpham)
          bigm     = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
             dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
             dabarm    = s(i,j,k) - s(i,j,k-1)
             dabarp    = s(i,j,k+1) - s(i,j,k)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2     = SIX*(alpham + alphap)
             D2L    = s(i,j,k-2)-TWO*s(i,j,k-1)+s(i,j,k)
             D2R    = s(i,j,k)-TWO*s(i,j,k+1)+s(i,j,k+2)
             D2C    = s(i,j,k-1)-TWO*s(i,j,k)+s(i,j,k+1)
             sgn    = sign(ONE,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(ONE,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j,k-1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(ONE,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j,k+1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -TWO*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k) + alpham
          sp(i,j) = s(i,j,k) + alphap

          if (ppm_flatten_before_integrals > 0) then
             ! flatten the parabola AFTER doing the monotonization (note k = k3d here)
             sm(i,j) = flatn(i,j,k3d)*sm(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
             sp(i,j) = flatn(i,j,k3d)*sp(i,j) + (ONE-flatn(i,j,k3d))*s(i,j,k3d)
          endif

          !
          ! Compute z-component of Ip and Im.
          !
          s6    = SIX*s(i,j,k3d) - THREE*(sm(i,j)+sp(i,j))

          
          ! w-c wave
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt/dz
          
          if (u(i,j,k3d,3)-cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,1) = sp(i,j) 
          else
             Ip(i,j,kc,3,1) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3)-cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,1) = sm(i,j) 
          else
             Im(i,j,kc,3,1) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! w wave
          sigma = abs(u(i,j,k3d,3))*dt/dz

          if (u(i,j,k3d,3) <= ZERO) then
             Ip(i,j,kc,3,2) = sp(i,j)
          else
             Ip(i,j,kc,3,2) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3) >= ZERO) then
             Im(i,j,kc,3,2) = sm(i,j)
          else
             Im(i,j,kc,3,2) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! w+c wave
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt/dz

          if (u(i,j,k3d,3)+cspd(i,j,k3d) <= ZERO) then
             Ip(i,j,kc,3,3) = sp(i,j) 
          else
             Ip(i,j,kc,3,3) = sp(i,j) - &
                  HALF*sigma*(sp(i,j)-sm(i,j)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i,j,k3d,3)+cspd(i,j,k3d) >= ZERO) then
             Im(i,j,kc,3,3) = sm(i,j) 
          else
             Im(i,j,kc,3,3) = sm(i,j) + &
                  HALF*sigma*(sp(i,j)-sm(i,j)+(ONE-TWO3RD*sigma)*s6)
          endif

       end do
    end do

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type2

end module ppm_module

