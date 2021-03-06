! These routines do the characteristic tracing under the parabolic
! profiles in each zone to the edge / half-time.

module trace_ppm_rad_module

  use bl_fort_module, only : rt => c_real
  implicit none

  private

  public trace_ppm_rad

contains

  subroutine trace_ppm_rad(q, qaux, flatn, qd_l1, qd_l2, qd_h1, qd_h2, &
                           dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                           qxm, qxp, qym, qyp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                           srcQ, src_l1, src_l2, src_h1, src_h2, &
                           ilo1, ilo2, ihi1, ihi2, dx, dy, dt)

    use network, only : nspec
    use bl_constants_module
    use meth_params_module, only : NQ, NQAUX, QVAR, QRHO, QU, QV, QREINT, QPRES, &
         QGAME, QC, QCG, QGAMC, QGAMCG, QLAMS, &
         QRADVAR, qrad, qradhi, qptot, qreitot, &
         small_dens, small_pres, &
         ppm_type, ppm_trace_sources, ppm_temp_fix, &
         ppm_reference_eigenvectors, &
         ppm_predict_gammae, &
         npassive, qpass_map
    use rad_params_module, only : ngroups
    use ppm_module, only : ppm

    use bl_fort_module, only : rt => c_real
    implicit none

    integer ilo1,ilo2,ihi1,ihi2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer gc_l1,gc_l2,gc_h1,gc_h2

    real(rt)             q(qd_l1:qd_h1,qd_l2:qd_h2,NQ)
    real(rt)         :: qaux(qd_l1:qd_h1,qd_l2:qd_h2, NQAUX)
    real(rt)         flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    real(rt)         srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    real(rt)         dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)

    real(rt)         qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)
    real(rt)         qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)
    real(rt)         qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)
    real(rt)         qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,NQ)

    real(rt)         dx, dy, dt

    ! Local variables
    integer :: i, j, g
    integer :: n, ipassive

    real(rt)         :: hdt, dtdx, dtdy

    ! To allow for easy integration of radiation, we adopt the
    ! following conventions:
    !
    ! rho : mass density
    ! u, v, w : velocities
    ! p : gas (hydro) pressure
    ! ptot : total pressure (note for pure hydro, this is 
    !        just the gas pressure)
    ! rhoe_g : gas specific internal energy
    ! cgas : sound speed for just the gas contribution
    ! cc : total sound speed (including radiation)
    ! h_g : gas specific enthalpy / cc**2
    ! gam_g : the gas Gamma_1
    ! game : gas gamma_e
    !
    ! for pure hydro, we will only consider:
    !   rho, u, v, w, ptot, rhoe_g, cc, h_g

    real(rt)         :: cc, csq, cgassq, Clag
    real(rt)         :: rho, u, v, p, rhoe_g, h_g, tau
    real(rt)         :: ptot, gam_g, game

    real(rt)         :: drho, dptot, drhoe_g
    real(rt)         :: dge, dtau
    real(rt)         :: dup, dvp, dptotp
    real(rt)         :: dum, dvm, dptotm

    real(rt)         :: rho_ref, u_ref, v_ref, p_ref, rhoe_g_ref, h_g_ref
    real(rt)         :: ptot_ref
    real(rt)         :: tau_ref

    real(rt)         :: cc_ref, csq_ref, Clag_ref, gam_g_ref, game_ref, gfactor

    real(rt)         :: alpham, alphap, alpha0r, alpha0e_g
    real(rt)         :: sourcr,sourcp,source,courn,eta,dlogatmp

    real(rt)         :: tau_s

    real(rt)        , dimension(0:ngroups-1) :: er, der, alphar, sourcer, qrtmp,hr
    real(rt)        , dimension(0:ngroups-1) :: lam0, lamp, lamm

    real(rt)        , dimension(0:ngroups-1) :: er_ref


    real(rt)        , allocatable :: Ip(:,:,:,:,:)
    real(rt)        , allocatable :: Im(:,:,:,:,:)

    real(rt)        , allocatable :: Ip_src(:,:,:,:,:)
    real(rt)        , allocatable :: Im_src(:,:,:,:,:)

    real(rt)         :: er_foo

    if (ppm_type == 0) then
       print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
    end if

    if (ppm_temp_fix > 0) then
       call bl_error("ERROR: ppm_temp_fix > 0 not implemented with radiation")
    endif

    if (ppm_reference_eigenvectors == 1) then
       call bl_error("ERROR: ppm_reference_eigenvectors not implemented with radiation")
    endif

    dtdx = dt/dx
    dtdy = dt/dy

    ! indices: (x, y, dimension, wave, variable)
    allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,NQ))
    allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,NQ))

    if (ppm_trace_sources == 1) then
       allocate(Ip_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
       allocate(Im_src(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QVAR))
    endif

    !allocate(Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))
    !allocate(Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,1))

    hdt = HALF * dt


    !=========================================================================
    ! PPM CODE
    !=========================================================================

    ! This does the characteristic tracing to build the interface
    ! states using the normal predictor only (no transverse terms).
    !
    ! We first fill the Im and Ip arrays -- these are the averages of
    ! the various primitive state variables under the parabolic
    ! interpolant over the region swept out by one of the 3 different
    ! characteristic waves.
    !
    ! Im is integrating to the left interface of the current zone
    ! (which will be used to build the right state at that interface)
    ! and Ip is integrating to the right interface of the current zone
    ! (which will be used to build the left state at that interface).
    !
    ! The indices are: Ip(i, j, dim, wave, var)
    !
    ! The choice of reference state is designed to minimize the
    ! effects of the characteristic projection.  We subtract the I's
    ! off of the reference state, project the quantity such that it is
    ! in terms of the characteristic varaibles, and then add all the
    ! jumps that are moving toward the interface to the reference
    ! state to get the full state on that interface.


    ! Compute Ip and Im -- this does the parabolic reconstruction,
    ! limiting, and returns the integral of each profile under each
    ! wave to each interface
    do n = 1, NQ
       call ppm(q(:,:,n), qd_l1, qd_l2, qd_h1, qd_h2, &
                q(:,:,QU:QV), qaux(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2,&
                flatn, &
                Ip(:,:,:,:,n), Im(:,:,:,:,n), &
                ilo1, ilo2, ihi1, ihi2, dx, dy, dt)
    end do

    ! trace the gas gamma to the edge
    !if (ppm_temp_fix /= 1) then
    !      call ppm(gamc(:,:),gc_l1,gc_l2,gc_h1,gc_h2, &
    !               q(:,:,QU:QV),c,qd_l1,qd_l2,qd_h1,qd_h2, &
    !               flatn, &
    !               Ip_gc(:,:,:,:,1),Im_gc(:,:,:,:,1), &
    !               ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    !   endif

    if (ppm_trace_sources == 1) then
       do n = 1, QVAR
          call ppm(srcQ(:,:,n), src_l1, src_l2, src_h1, src_h2, &
                   q(:,:,QU:QV), qaux(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2, &
                   flatn, &
                   Ip_src(:,:,:,:,n), Im_src(:,:,:,:,n), &
                   ilo1, ilo2, ihi1, ihi2, dx, dy, dt)
       enddo
    endif

    !-------------------------------------------------------------------------
    ! x-direction
    !-------------------------------------------------------------------------

    ! Trace to left and right edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = qaux(i,j,QLAMS+g)
             lamp(g) = qaux(i,j,QLAMS+g)
             lamm(g) = qaux(i,j,QLAMS+g)
          end do

          ! cgassq is the gas soundspeed **2
          ! cc is the total soundspeed **2 (gas + radiation)
          cgassq = qaux(i,j,QCG)**2
          cc = qaux(i,j,QC)
          csq = cc**2

          rho = q(i,j,QRHO)
          tau = ONE/rho  ! should not be needed once reference ev is implemented
          u = q(i,j,QU)
          v = q(i,j,QV)

          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = ( (p+rhoe_g) / rho)/csq

          Clag = rho*cc

          gam_g = qaux(i,j,QGAMCG)
          game = q(i,j,QGAME)

          ptot = q(i,j,qptot)

          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

          !-------------------------------------------------------------------
          ! plus state on face i
          !-------------------------------------------------------------------

          ! set the reference state
          ! this will be the fastest moving state to the left --
          ! this is the method that Miller & Colella and Colella &
          ! Woodward use
          rho_ref  = Im(i,j,1,1,QRHO)
          u_ref    = Im(i,j,1,1,QU)

          p_ref    = Im(i,j,1,1,QPRES)
          rhoe_g_ref = Im(i,j,1,1,QREINT)

          tau_ref = ONE/Im(i,j,1,1,QRHO)

          !gam_g_ref  = Im_gc(i,j,1,1,1)
          game_ref  = Im(i,j,1,1,QGAME)

          ptot_ref = Im(i,j,1,1,qptot)

          er_ref(:) = Im(i,j,1,1,qrad:qradhi)

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! *m are the jumps carried by u-c
          ! *p are the jumps carried by u+c

          dum    = u_ref    - Im(i,j,1,1,QU)
          dptotm = ptot_ref - Im(i,j,1,1,qptot)

          drho    = rho_ref    - Im(i,j,1,2,QRHO)
          dptot   = ptot_ref   - Im(i,j,1,2,qptot)
          drhoe_g = rhoe_g_ref - Im(i,j,1,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,1,2,QRHO)
          der(:)  = er_ref(:)  - Im(i,j,1,2,qrad:qradhi)

          dup    = u_ref    - Im(i,j,1,3,QU)
          dptotp = ptot_ref - Im(i,j,1,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - hdt*Im_src(i,j,1,1,QU)
             dup = dup - hdt*Im_src(i,j,1,3,QU)
          endif


          ! optionally use the reference state in evaluating the
          ! eigenvectors


          if (ppm_predict_gammae == 0) then

             ! (rho, u, p, (rho e)) eigensystem

             ! these are analogous to the beta's from the original
             ! PPM paper (except we work with rho instead of tau).
             ! This is simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e_g = drhoe_g - dptot*h_g

          else

             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             alpham = HALF*( dum - dptotm/Clag)/Clag
             alphap = HALF*(-dup - dptotp/Clag)/Clag
             alpha0r = dtau + dptot/Clag**2

             dge = game_ref - Im(i,j,1,2,QGAME)
             gfactor = (game - ONE)*(game - gam_g)
             alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

          endif ! which tracing method

          alphar(:) = der(:) - dptot/csq*hr

          if (u-cc > ZERO) then
             alpham = ZERO
          else if (u-cc < ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif

          if (u+cc > ZERO) then
             alphap = ZERO
          else if (u+cc < ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif

          if (u > ZERO) then
             alpha0r = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else if (u < ZERO) then
             alpha0r = -alpha0r
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else
             alpha0r = -HALF*alpha0r
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i >= ilo1) then

             if (ppm_predict_gammae == 0) then
                qxp(i,j,QRHO)   = rho_ref + alphap + alpham + alpha0r
                qxp(i,j,QU)     = u_ref + (alphap - alpham)*cc/rho
                qxp(i,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                qxp(i,j,QPRES)  = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

                qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                qxp(i,j,qrad:qradhi) = qrtmp

                qxp(i,j,qptot) = ptot_ref + (alphap + alpham)*csq
                qxp(i,j,qreitot) = qxp(i,j,QREINT) + sum(qrtmp)

             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qxp(i,j,QRHO)   = ONE/tau_s

                qxp(i,j,QU)     = u_ref + (alpham - alphap)*Clag
                qxp(i,j,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamp(:)*alphar(:))

                qxp(i,j,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                qxp(i,j,QREINT) = qxp(i,j,QPRES )/(qxp(i,j,QGAME) - ONE)

                qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                qxp(i,j,qrad:qradhi) = qrtmp

                qxp(i,j,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                qxp(i,j,qreitot) = qxp(i,j,QREINT) + sum(qrtmp)

             endif

             ! enforce small_*
             qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
             qxp(i,j,QPRES) = max(qxp(i,j,QPRES), small_pres)

             do g=0, ngroups-1
                if (qxp(i,j,qrad+g) < ZERO) then
                   er_foo = - qxp(i,j,qrad+g)
                   qxp(i,j,qrad+g) = ZERO
                   qxp(i,j,qptot) = qxp(i,j,qptot) + lamp(g) * er_foo
                   qxp(i,j,qreitot) = qxp(i,j,qreitot) + er_foo
                end if
             end do

             if (qxp(i,j,QPRES) < ZERO) then
                qxp(i,j,QPRES) = p
             end if


             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave

             ! Recall that I already takes the limit of the parabola
             ! in the event that the wave is not moving toward the
             ! interface
             qxp(i,j,QV) = Im(i,j,1,2,QV)

             if (ppm_trace_sources == 1) then
                qxp(i,j,QV) = qxp(i,j,QV) + hdt*Im_src(i,j,1,2,QV)
             endif

          end if


          !-------------------------------------------------------------------
          ! minus state on face i+1
          !-------------------------------------------------------------------

          ! set the reference state
          ! this will be the fastest moving state to the right
          rho_ref  = Ip(i,j,1,3,QRHO)
          u_ref    = Ip(i,j,1,3,QU)

          p_ref      = Ip(i,j,1,3,QPRES)
          rhoe_g_ref = Ip(i,j,1,3,QREINT)

          tau_ref = ONE/Ip(i,j,1,3,QRHO)

          !gam_g_ref    = Ip_gc(i,j,1,3,1)
          game_ref    = Ip(i,j,1,3,QGAME)

          ptot_ref = Ip(i,j,1,3,qptot)

          er_ref(:) = Ip(i,j,1,3,qrad:qradhi)


          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          !  *m are the jumps carried by u-c
          !  *p are the jumps carried by u+c

          dum    = u_ref    - Ip(i,j,1,1,QU)
          dptotm = ptot_ref - Ip(i,j,1,1,qptot)

          drho    = rho_ref    - Ip(i,j,1,2,QRHO)
          dptot   = ptot_ref   - Ip(i,j,1,2,qptot)
          drhoe_g = rhoe_g_ref - Ip(i,j,1,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,1,2,QRHO)
          der(:)  = er_ref(:)  - Ip(i,j,1,2,qrad:qradhi)

          dup    = u_ref    - Ip(i,j,1,3,QU)
          dptotp = ptot_ref - Ip(i,j,1,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dum = dum - hdt*Ip_src(i,j,1,1,QU)
             dup = dup - hdt*Ip_src(i,j,1,3,QU)
          endif

          ! optionally use the reference state in evaluating the
          ! eigenvectors

          if (ppm_predict_gammae == 0) then

             ! (rho, u, p, (rho e)) eigensystem

             ! these are analogous to the beta's from the original
             ! PPM paper (except we work with rho instead of tau).
             ! This is simply (l . dq), where dq = qref - I(q)

             alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e_g = drhoe_g - dptot*h_g

          else

             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             alpham = HALF*( dum - dptotm/Clag)/Clag
             alphap = HALF*(-dup - dptotp/Clag)/Clag
             alpha0r = dtau + dptot/Clag**2

             dge = game_ref - Ip(i,j,1,2,QGAME)
             gfactor = (game - ONE)*(game - gam_g)
             alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

          endif

          alphar(:) = der(:) - dptot/csq*hr

          if (u-cc > ZERO) then
             alpham = -alpham
          else if (u-cc < ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif

          if (u+cc > ZERO) then
             alphap = -alphap
          else if (u+cc < ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif

          if (u > ZERO) then
             alpha0r = -alpha0r
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else if (u < ZERO) then
             alpha0r = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else
             alpha0r = -HALF*alpha0r
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (i <= ihi1) then

             if (ppm_predict_gammae == 0) then

                qxm(i+1,j,QRHO) = rho_ref + alphap + alpham + alpha0r
                qxm(i+1,j,QU) = u_ref + (alphap - alpham)*cc/rho
                qxm(i+1,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                qxm(i+1,j,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))

                qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                qxm(i+1,j,qrad:qradhi) = qrtmp

                qxm(i+1,j,qptot) = ptot_ref + (alphap + alpham)*csq
                qxm(i+1,j,qreitot) = qxm(i+1,j,QREINT) + sum(qrtmp)

             else

                tau_s = tau_ref + alphap + alpham + alpha0r
                qxm(i+1,j,QRHO)   = ONE/tau_s

                qxm(i+1,j,QU)     = u_ref + (alpham - alphap)*Clag
                qxm(i+1,j,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:))

                qxm(i+1,j,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                qxm(i+1,j,QREINT) = qxm(i+1,j,QPRES )/(qxm(i+1,j,QGAME) - ONE)

                qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                qxm(i+1,j,qrad:qradhi) = qrtmp

                qxm(i+1,j,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                qxm(i+1,j,qreitot) = qxm(i+1,j,QREINT) + sum(qrtmp)

             endif

             ! enforce small_*
             qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
             qxm(i+1,j,QPRES) = max(qxm(i+1,j,QPRES), small_pres)

             do g=0, ngroups-1
                if (qxm(i+1,j,qrad+g) < ZERO) then
                   er_foo = - qxm(i+1,j,qrad+g)
                   qxm(i+1,j,qrad+g) = ZERO
                   qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + lamm(g) * er_foo
                   qxm(i+1,j,qreitot) = qxm(i+1,j,qreitot) + er_foo
                end if
             end do

             if (qxm(i+1,j,QPRES) < ZERO) then
                qxm(i+1,j,QPRES) = p
             end if

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave
             qxm(i+1,j,QV) = Ip(i,j,1,2,QV)

             if (ppm_trace_sources == 1) then
                qxm(i+1,j,QV) = qxm(i+1,j,QV) + hdt*Ip_src(i,j,1,2,QV)
             endif

          end if

          !-------------------------------------------------------------------
          ! geometry source terms
          !-------------------------------------------------------------------

          if (dloga(i,j) /= 0) then
             courn = dtdx*(cc+abs(u))
             eta = (ONE-courn)/(cc*dt*abs(dloga(i,j)))
             dlogatmp = min(eta,ONE)*dloga(i,j)
             sourcr = -HALF*dt*rho*dlogatmp*u
             sourcp = sourcr*cgassq
             source = sourcr*h_g*csq
             sourcer(:) = -HALF*dt*dlogatmp*u*(lam0(:)+ONE)*er(:)
             if (i <= ihi1) then
                qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
                qxm(i+1,j,qrad:qradhi) = qxm(i+1,j,qrad:qradhi) + sourcer(:)
                !              qxm(i+1,j,qptot ) = sum(lamm(:)*qxm(i+1,j,qrad:qradhi)) + qxm(i+1,j,QPRES)
                qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                qxm(i+1,j,qreitot) = sum(qxm(i+1,j,qrad:qradhi))  + qxm(i+1,j,QREINT)
             end if
             if (i >= ilo1) then
                qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
                qxp(i,j,qrad:qradhi) = qxp(i,j,qrad:qradhi) + sourcer(:)
                !              qxp(i,j,qptot ) = sum(lamp(:)*qxp(i,j,qrad:qradhi)) + qxp(i,j,QPRES)
                qxp(i,j,qptot) = qxp(i,j,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                qxp(i,j,qreitot) = sum(qxp(i,j,qrad:qradhi))  + qxp(i,j,QREINT)
             end if
          endif

       end do
    end do


    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------

    ! We do all passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,QU)

             ! We have
             !
             ! q_l = q_ref - Proj{(q_ref - I)}
             !
             ! and Proj{} represents the characteristic projection.
             ! But for these, there is only 1-wave that matters, the u
             ! wave, so no projection is needed.  Since we are not
             ! projecting, the reference state doesn't matter.

             if (u > ZERO) then
                qxp(i,j,n) = q(i,j,n)    ! we might want to change this to
                                         ! the limit of the parabola
             else if (u < ZERO) then
                qxp(i,j,n) = Im(i,j,1,2,n)
             else
                qxp(i,j,n) = q(i,j,n) + HALF*(Im(i,j,1,2,n) - q(i,j,n))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,QU)

             if (u > ZERO) then
                qxm(i+1,j,n) = Ip(i,j,1,2,n)
             else if (u < ZERO) then
                qxm(i+1,j,n) = q(i,j,n)
             else
                qxm(i+1,j,n) = q(i,j,n) + HALF*(Ip(i,j,1,2,n) - q(i,j,n))
             endif
          enddo

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! y-direction
    !-------------------------------------------------------------------------

    ! Trace to bottom and top edges using upwind PPM
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          do g=0, ngroups-1
             lam0(g) = qaux(i,j,QLAMS+g)
             lamp(g) = qaux(i,j,QLAMS+g)
             lamm(g) = qaux(i,j,QLAMS+g)
          end do

          ! cgassq is the gas soundspeed **2
          ! cc is the total soundspeed **2 (gas + radiation)
          cgassq = qaux(i,j,QCG)**2
          cc = qaux(i,j,QC)
          csq = cc**2

          rho = q(i,j,QRHO)
          tau = ONE/rho
          u = q(i,j,QU)
          v = q(i,j,QV)

          p = q(i,j,QPRES)
          rhoe_g = q(i,j,QREINT)
          h_g = ( (p+rhoe_g) / rho)/csq

          Clag = rho*cc

          gam_g = qaux(i,j,QGAMCG)
          game = q(i,j,QGAME)

          ptot = q(i,j,qptot)

          er(:) = q(i,j,qrad:qradhi)
          hr(:) = (lam0+ONE)*er/rho

          !-------------------------------------------------------------------
          ! plus state on face j
          !-------------------------------------------------------------------

          ! set the reference state
          ! this will be the fastest moving state to the left
          rho_ref  = Im(i,j,2,1,QRHO)
          v_ref    = Im(i,j,2,1,QV)

          p_ref = Im(i,j,2,1,QPRES)
          rhoe_g_ref = Im(i,j,2,1,QREINT)

          tau_ref = ONE/Im(i,j,2,1,QRHO)

          game_ref = Im(i,j,2,1,QGAME)

          ptot_ref = Im(i,j,2,1,qptot)

          er_ref(:) = Im(i,j,2,1,qrad:qradhi)

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c
          
          dvm    = v_ref    - Im(i,j,2,1,QV)
          dptotm = ptot_ref - Im(i,j,2,1,qptot)

          drho    = rho_ref    - Im(i,j,2,2,QRHO)
          dptot   = ptot_ref   - Im(i,j,2,2,qptot)
          drhoe_g = rhoe_g_ref - Im(i,j,2,2,QREINT)
          dtau  = tau_ref  - ONE/Im(i,j,2,2,QRHO)
          der(:)  = er_ref(:)  - Im(i,j,2,2,qrad:qradhi)

          dvp    = v_ref    - Im(i,j,2,3,QV)
          dptotp = ptot_ref - Im(i,j,2,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - hdt*Im_src(i,j,2,1,QV)
             dvp = dvp - hdt*Im_src(i,j,2,3,QV)
          endif

          ! optionally use the reference state in evaluating the
          ! eigenvectors

          if (ppm_predict_gammae == 0) then

             ! (rho, u, p, (rho e)) eigensystem

             ! these are analogous to the beta's from the original PPM
             ! paper (except we work with rho instead of tau).  This
             ! is simply (l . dq), where dq = qref - I(q)
             alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e_g = drhoe_g - dptot*h_g

          else
             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             alpham = HALF*( dvm - dptotm/Clag)/Clag
             alphap = HALF*(-dvp - dptotp/Clag)/Clag
             alpha0r = dtau + dptot/Clag**2

             dge = game_ref - Im(i,j,2,2,QGAME)
             gfactor = (game - ONE)*(game - gam_g)
             alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

          endif

          alphar(:) = der(:) - dptot/csq*hr

          if (v-cc > ZERO) then
             alpham = ZERO
          else if (v-cc < ZERO) then
             alpham = -alpham
          else
             alpham = -HALF*alpham
          endif

          if (v+cc > ZERO) then
             alphap = ZERO
          else if (v+cc < ZERO) then
             alphap = -alphap
          else
             alphap = -HALF*alphap
          endif

          if (v > ZERO) then
             alpha0r = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else if (v < ZERO) then
             alpha0r = -alpha0r
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else
             alpha0r = -HALF*alpha0r
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (j >= ilo2) then
             if (ppm_predict_gammae == 0) then
                qyp(i,j,QRHO) = rho_ref + alphap + alpham + alpha0r
                qyp(i,j,QV) = v_ref + (alphap - alpham)*cc/rho
                qyp(i,j,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                qyp(i,j,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

                qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                qyp(i,j,qrad:qradhi) = qrtmp

                qyp(i,j,qptot) = ptot_ref + (alphap + alpham)*csq
                qyp(i,j,qreitot) = qyp(i,j,QREINT) + sum(qrtmp)

             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qyp(i,j,QRHO)   = ONE/tau_s

                qyp(i,j,QV)     = v_ref + (alpham - alphap)*Clag
                qyp(i,j,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamp(:)*alphar(:))

                qyp(i,j,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                qyp(i,j,QREINT) = qyp(i,j,QPRES )/(qyp(i,j,QGAME) - ONE)

                qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                qyp(i,j,qrad:qradhi) = qrtmp

                qyp(i,j,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                qyp(i,j,qreitot) = qyp(i,j,QREINT) + sum(qrtmp)

             endif

             ! enforce small_*
             qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
             qyp(i,j,QPRES) = max(qyp(i,j,QPRES), small_pres)

             do g=0, ngroups-1
                if (qyp(i,j,qrad+g) < ZERO) then
                   er_foo = - qyp(i,j,qrad+g)
                   qyp(i,j,qrad+g) = ZERO
                   qyp(i,j,qptot) = qyp(i,j,qptot) + lamp(g) * er_foo
                   qyp(i,j,qreitot) = qyp(i,j,qreitot) + er_foo
                end if
             end do

             if (qyp(i,j,QPRES) < ZERO) then
                qyp(i,j,QPRES) = p
             end if

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave
             qyp(i,j,QU)     = Im(i,j,2,2,QU)

             if (ppm_trace_sources == 1) then
                qyp(i,j,QU) = qyp(i,j,QU) + hdt*Im_src(i,j,2,2,QU)
             endif

          end if

          !-------------------------------------------------------------------
          ! minus state on face j+1
          !-------------------------------------------------------------------

          ! set the reference state
          ! this will be the fastest moving state to the right
          rho_ref = Ip(i,j,2,3,QRHO)
          v_ref   = Ip(i,j,2,3,QV)

          p_ref      = Ip(i,j,2,3,QPRES)             
          rhoe_g_ref = Ip(i,j,2,3,QREINT)

          tau_ref = ONE/Ip(i,j,2,3,QRHO)

          !gam_g_ref    = Ip_gc(i,j,2,3,1)
          game_ref    = Ip(i,j,2,3,QGAME)

          ptot_ref = Ip(i,j,2,3,qptot)

          er_ref(:) = Ip(i,j,2,3,qrad:qradhi)

          rho_ref = max(rho_ref,small_dens)
          p_ref = max(p_ref,small_pres)

          ! *m are the jumps carried by v-c
          ! *p are the jumps carried by v+c

          dvm    = v_ref    - Ip(i,j,2,1,QV)
          dptotm = ptot_ref - Ip(i,j,2,1,qptot)

          drho    = rho_ref    - Ip(i,j,2,2,QRHO)
          dptot   = ptot_ref   - Ip(i,j,2,2,qptot)
          drhoe_g = rhoe_g_ref - Ip(i,j,2,2,QREINT)
          dtau  = tau_ref  - ONE/Ip(i,j,2,2,QRHO)
          der(:)  = er_ref(:)  - Ip(i,j,2,2,qrad:qradhi)

          dvp    = v_ref    - Ip(i,j,2,3,QV)
          dptotp = ptot_ref - Ip(i,j,2,3,qptot)

          ! if we are doing source term tracing, then we add the force to
          ! the velocity here, otherwise we will deal with this in the
          ! trans_X routines
          if (ppm_trace_sources == 1) then
             dvm = dvm - hdt*Ip_src(i,j,2,1,QV)
             dvp = dvp - hdt*Ip_src(i,j,2,3,QV)
          endif

          ! optionally use the reference state in evaluating the
          ! eigenvectors

          if (ppm_predict_gammae == 0) then

             ! (rho, u, p, (rho e) eigensystem

             ! these are analogous to the beta's from the original PPM
             ! paper.  This is simply (l . dq), where dq = qref - I(q)
             alpham = HALF*(dptotm/(rho*cc) - dvm)*rho/cc
             alphap = HALF*(dptotp/(rho*cc) + dvp)*rho/cc
             alpha0r = drho - dptot/csq
             alpha0e_g = drhoe_g - dptot*h_g

          else

             ! (tau, u, p, game) eigensystem

             ! this is the way things were done in the original PPM
             ! paper -- here we work with tau in the characteristic
             ! system.

             alpham = HALF*( dvm - dptotm/Clag)/Clag
             alphap = HALF*(-dvp - dptotp/Clag)/Clag
             alpha0r = dtau + dptot/Clag**2

             dge = game_ref - Ip(i,j,2,2,QGAME)
             gfactor = (game - ONE)*(game - gam_g)
             alpha0e_g = gfactor*dptot/(tau*Clag**2) + dge

          endif

          alphar(:) = der(:)- dptot/csq*hr

          if (v-cc > ZERO) then
             alpham = -alpham
          else if (v-cc < ZERO) then
             alpham = ZERO
          else
             alpham = -HALF*alpham
          endif

          if (v+cc > ZERO) then
             alphap = -alphap
          else if (v+cc < ZERO) then
             alphap = ZERO
          else
             alphap = -HALF*alphap
          endif

          if (v > ZERO) then
             alpha0r = -alpha0r
             alpha0e_g = -alpha0e_g
             alphar(:) = -alphar(:)
          else if (v < ZERO) then
             alpha0r = ZERO
             alpha0e_g = ZERO
             alphar(:) = ZERO
          else
             alpha0r = -HALF*alpha0r
             alpha0e_g = -HALF*alpha0e_g
             alphar(:) = -HALF*alphar(:)
          endif

          ! the final interface states are just
          ! q_s = q_ref - sum (l . dq) r
          if (j <= ihi2) then
             if (ppm_predict_gammae == 0) then
                qym(i,j+1,QRHO) = rho_ref + alphap + alpham + alpha0r
                qym(i,j+1,QV) = v_ref + (alphap - alpham)*cc/rho
                qym(i,j+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                qym(i,j+1,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
                
                qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                qym(i,j+1,qrad:qradhi) = qrtmp

                qym(i,j+1,qptot) = ptot_ref + (alphap + alpham)*csq
                qym(i,j+1,qreitot) = qym(i,j+1,QREINT) + sum(qrtmp)

             else
                tau_s = tau_ref + alphap + alpham + alpha0r
                qym(i,j+1,QRHO)   = ONE/tau_s

                qym(i,j+1,QV)     = v_ref + (alpham - alphap)*Clag
                qym(i,j+1,QPRES)  = p_ref - (alphap + alpham)*(cgassq/tau**2) - sum(lamm(:)*alphar(:))

                qym(i,j+1,QGAME) = game_ref + gfactor*(alpham + alphap)/tau + alpha0e_g
                qym(i,j+1,QREINT) = qym(i,j+1,QPRES )/(qym(i,j+1,QGAME) - ONE)

                qrtmp = er_ref(:) - (alphap + alpham)*hr/tau**2 + alphar(:)
                qym(i,j+1,qrad:qradhi) = qrtmp

                qym(i,j+1,qptot) = ptot_ref - (alphap + alpham)*Clag**2
                qym(i,j+1,qreitot) = qym(i,j+1,QREINT) + sum(qrtmp)

             endif

             ! enforce small_*
             qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
             qym(i,j+1,QPRES) = max(qym(i,j+1,QPRES), small_pres)

             do g=0, ngroups-1
                if (qym(i,j+1,qrad+g) < ZERO) then
                   er_foo = - qym(i,j+1,qrad+g)
                   qym(i,j+1,qrad+g) = ZERO
                   qym(i,j+1,qptot) = qym(i,j+1,qptot) + lamm(g) * er_foo
                   qym(i,j+1,qreitot) = qym(i,j+1,qreitot) + er_foo
                end if
             end do

             if (qym(i,j+1,QPRES) < ZERO) then
                qym(i,j+1,QPRES) = p
             end if

             ! transverse velocity -- there is no projection here, so
             ! we don't need a reference state.  We only care about
             ! the state traced under the middle wave
             qym(i,j+1,QU)     = Ip(i,j,2,2,QU)
             
             if (ppm_trace_sources == 1) then
                qym(i,j+1,QU) = qym(i,j+1,QU) + hdt*Ip_src(i,j,2,2,QU)
             endif

          end if

       end do
    end do


    !-------------------------------------------------------------------------
    ! Now do the passively advected quantities
    !-------------------------------------------------------------------------

    ! do all of the passively advected quantities in one loop
    do ipassive = 1, npassive
       n = qpass_map(ipassive)
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,QV)

             if (v > ZERO) then
                qyp(i,j,n) = q(i,j,n)
             else if (v < ZERO) then
                qyp(i,j,n) = Im(i,j,2,2,n)
             else
                qyp(i,j,n) = q(i,j,n) + HALF*(Im(i,j,2,2,n) - q(i,j,n))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,QV)

             if (v > ZERO) then
                qym(i,j+1,n) = Ip(i,j,2,2,n)
             else if (v < ZERO) then
                qym(i,j+1,n) = q(i,j,n)
             else
                qym(i,j+1,n) = q(i,j,n) + HALF*(Ip(i,j,2,2,n) - q(i,j,n))
             endif
          enddo

       enddo
    enddo

    deallocate(Ip,Im)
    if (ppm_trace_sources == 1) then
       deallocate(Ip_src,Im_src)
    endif

  end subroutine trace_ppm_rad

end module trace_ppm_rad_module
