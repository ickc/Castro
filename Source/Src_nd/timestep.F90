module timestep_module

  use bl_fort_module, only : rt => c_real
  implicit none

  public

contains

  ! Courant-condition limited timestep

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="ca_estdt")

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX
    use prob_params_module, only: dim
    use bl_constants_module
#ifdef ROTATION
    use meth_params_module, only: do_rotation, state_in_rotating_frame
    use rotation_module, only: inertial_to_rotational_velocity
    use amrinfo_module, only: amr_time
#endif

    use bl_fort_module, only : rt => c_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt)         :: dx(3), dt

    real(rt)         :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3
    integer          :: i, j, k

    type (eos_t) :: eos_state

#ifdef ROTATION
    real(rt)         :: vel(3)
#endif

    ! Call EOS for the purpose of computing sound speed

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho = u(i,j,k,URHO )
             eos_state % T   = u(i,j,k,UTEMP)
             eos_state % e   = u(i,j,k,UEINT) * rhoInv
             eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             ! Compute velocity and then calculate CFL timestep.

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

#ifdef ROTATION
             if (do_rotation .eq. 1 .and. state_in_rotating_frame .ne. 1) then
                vel = [ux, uy, uz]
                call inertial_to_rotational_velocity([i, j, k], amr_time, vel)
                ux = vel(1)
                uy = vel(2)
                uz = vel(3)
             endif
#endif
             
             c = eos_state % cs

             dt1 = dx(1)/(c + abs(ux))
             if (dim .ge. 2) then
                dt2 = dx(2)/(c + abs(uy))
             else
                dt2 = dt1
             endif
             if (dim .eq. 3) then
                dt3 = dx(3)/(c + abs(uz))
             else
                dt3 = dt1
             endif

             dt  = min(dt,dt1,dt2,dt3)

          enddo
       enddo
    enddo

  end subroutine ca_estdt

  ! Reactions-limited timestep

#ifdef REACTIONS
  subroutine ca_estdt_burning(sold, so_lo, so_hi, &
                              snew, sn_lo, sn_hi, &
                              rold, ro_lo, ro_hi, &
                              rnew, rn_lo, rn_hi, &
                              lo, hi, dx, dt_old, dt) &
                              bind(C, name="ca_estdt_burning")

    use bl_constants_module, only: ONE
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, dtnuc_e, dtnuc_X, dtnuc_mode
    use prob_params_module, only : dim
#if naux > 0
    use meth_params_module, only : UFX
#endif
    use actual_rhs_module, only: actual_rhs
    use eos_module
    use burner_module, only: ok_to_burn
    use burn_type_module
    use eos_type_module
    use extern_probin_module, only: small_x

    use bl_fort_module, only : rt => c_real
    implicit none

    integer          :: so_lo(3), so_hi(3)
    integer          :: sn_lo(3), sn_hi(3)
    integer          :: ro_lo(3), ro_hi(3)
    integer          :: rn_lo(3), rn_hi(3)
    integer          :: lo(3), hi(3)
    real(rt)         :: sold(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),NVAR)
    real(rt)         :: snew(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
    real(rt)         :: rold(ro_lo(1):ro_hi(1),ro_lo(2):ro_hi(2),ro_lo(3):ro_hi(3),nspec+2)
    real(rt)         :: rnew(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3),nspec+2)
    real(rt)         :: dx(3), dt, dt_old

    real(rt)         :: e, X(nspec), dedt, dXdt(nspec)
    integer          :: i, j, k

    type (burn_t)    :: state_old, state_new
    type (eos_t)     :: eos_state
    real(rt)         :: rhooinv, rhoninv

    ! We want to limit the timestep so that it is not larger than
    ! dtnuc_e * (e / (de/dt)).  If the timestep factor dtnuc is
    ! equal to 1, this says that we don't want the
    ! internal energy to change by any more than its current
    ! magnitude in the next timestep. 
    !
    ! If dtnuc is less than one, it controls the fraction we will
    ! allow the internal energy to change in this timestep due to
    ! nuclear burning, provided that the last timestep's burning is a
    ! good estimate for the current timestep's burning.
    !
    ! We do the same thing for the species, using a timestep
    ! limiter dtnuc_X * (X_k / (dX_k/dt)).
    !
    ! To estimate de/dt and dX/dt, we are going to call the RHS of the
    ! burner given the current state data. We need to do an EOS
    ! call before we do the RHS call so that we have accurate
    ! values for the thermodynamic data like abar, zbar, etc.
    ! But we will call in (rho, T) mode, which is inexpensive.

    if (dtnuc_e > 1.e199_rt .and. dtnuc_X > 1.e199_rt) return

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhooinv = ONE / sold(i,j,k,URHO)
             rhoninv = ONE / snew(i,j,k,URHO)

             state_old % rho = sold(i,j,k,URHO)
             state_old % T   = sold(i,j,k,UTEMP)
             state_old % e   = sold(i,j,k,UEINT) * rhooinv
             state_old % xn  = sold(i,j,k,UFS:UFS+nspec-1) * rhooinv
#if naux > 0
             state_old % aux = sold(i,j,k,UFX:UFX+naux-1) * rhooinv
#endif

             state_new % rho = snew(i,j,k,URHO)
             state_new % T   = snew(i,j,k,UTEMP)
             state_new % e   = snew(i,j,k,UEINT) * rhoninv
             state_new % xn  = snew(i,j,k,UFS:UFS+nspec-1) * rhoninv
#if naux > 0
             state_new % aux = snew(i,j,k,UFX:UFX+naux-1) * rhoninv
#endif

             if (.not. ok_to_burn(state_new)) cycle

             e    = state_new % e
             X    = max(state_new % xn, small_x)

             if (dtnuc_mode == 1) then

                call burn_to_eos(state_new, eos_state)
                call eos(eos_input_rt, eos_state)
                call eos_to_burn(eos_state, state_new)

                state_new % dx = minval(dx(1:dim))

                call actual_rhs(state_new)

                dedt = state_new % ydot(net_ienuc)
                dXdt = state_new % ydot(1:nspec) * aion

             else if (dtnuc_mode == 2) then

                dedt = rnew(i,j,k,nspec+1)
                dXdt = rnew(i,j,k,1:nspec)

             else if (dtnuc_mode == 3) then

                dedt = HALF * (rold(i,j,k,nspec+1) + rnew(i,j,k,nspec+1))
                dXdt = HALF * (rold(i,j,k,1:nspec) + rnew(i,j,k,1:nspec))

             else if (dtnuc_mode == 4) then

                dedt = (state_new % e - state_old % e) / dt_old
                dXdt = (state_new % xn - state_old % xn) / dt_old

             else

                call bl_error("Error: unrecognized burning timestep limiter mode in timestep.F90.")

             endif

             dedt = max(abs(dedt), 1.e-50_rt)
             dXdt = max(abs(dXdt), 1.e-50_rt)

             dt = min(dt, dtnuc_e * e / dedt)
             dt = min(dt, dtnuc_X * minval(X / dXdt))

          enddo
       enddo
    enddo

  end subroutine ca_estdt_burning
#endif

  ! Diffusion-limited timestep

#ifdef DIFFUSION
  subroutine ca_estdt_temp_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) &
       bind(C, name="ca_estdt_temp_diffusion")

    use network, only: nspec, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
         diffuse_cutoff_density
    use prob_params_module, only: dim
    use bl_constants_module
    use conductivity_module

    use bl_fort_module, only : rt => c_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    real(rt)         :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)         :: dx(3), dt

    real(rt)         :: dt1, dt2, dt3, rho_inv
    integer          :: i, j, k
    real(rt)         :: cond, D

    type (eos_t) :: eos_state

    ! dt < 0.5 dx**2 / D
    ! where D = k/(rho c_v), and k is the conductivity

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) > diffuse_cutoff_density) then

                rho_inv = ONE/state(i,j,k,URHO)

                ! we need cv
                eos_state % rho = state(i,j,k,URHO )
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   = state(i,j,k,UEINT) * rho_inv

                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rho_inv
                eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rho_inv

                call eos(eos_input_re, eos_state)

                ! we also need the conductivity
                call thermal_conductivity(eos_state, cond)

                ! maybe we should check (and take action) on negative cv here?
                D = cond*rho_inv/eos_state%cv

                dt1 = HALF*dx(1)**2/D

                if (dim >= 2) then
                   dt2 = HALF*dx(2)**2/D
                else
                   dt2 = dt1
                endif

                if (dim == 3) then
                   dt3 = HALF*dx(3)**2/D
                else
                   dt3 = dt1
                endif

                dt  = min(dt,dt1,dt2,dt3)

             endif

          enddo
       enddo
    enddo

  end subroutine ca_estdt_temp_diffusion

  subroutine ca_estdt_enth_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) &
       bind(C, name="ca_estdt_enth_diffusion")

    use network, only: nspec, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
         diffuse_cutoff_density
    use prob_params_module, only: dim
    use bl_constants_module
    use conductivity_module

    use bl_fort_module, only : rt => c_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    real(rt)         :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)         :: dx(3), dt

    real(rt)         :: dt1, dt2, dt3, rho_inv
    integer          :: i, j, k
    real(rt)         :: cond, D

    type (eos_t) :: eos_state

    ! dt < 0.5 dx**2 / D
    ! where D = k/(rho c_p), and k is the conductivity

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) > diffuse_cutoff_density) then

                rho_inv = ONE/state(i,j,k,URHO)

                ! we need cv
                eos_state % rho = state(i,j,k,URHO )
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   = state(i,j,k,UEINT) * rho_inv

                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rho_inv
                eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rho_inv

                call eos(eos_input_re, eos_state)

                ! we also need the conductivity
                call thermal_conductivity(eos_state, cond)

                D = cond*rho_inv/eos_state%cp

                dt1 = HALF*dx(1)**2/D

                if (dim >= 2) then
                   dt2 = HALF*dx(2)**2/D
                else
                   dt2 = dt1
                endif

                if (dim == 3) then
                   dt3 = HALF*dx(3)**2/D
                else
                   dt3 = dt1
                endif

                dt  = min(dt,dt1,dt2,dt3)

             endif

          enddo
       enddo
    enddo

  end subroutine ca_estdt_enth_diffusion
#endif


  ! Check whether the last timestep violated any of our stability criteria.
  ! If so, suggest a new timestep which would not.

  subroutine ca_check_timestep(s_old, so_lo, so_hi, &
                               s_new, sn_lo, sn_hi, &
#ifdef REACTIONS
                               r_old, ro_lo, ro_hi, &
                               r_new, rn_lo, rn_hi, &
#endif
                               lo, hi, &
                               dx, dt_old, dt_new) &
                               bind(C, name="ca_check_timestep")

    use bl_constants_module, only: HALF, ONE
    use meth_params_module, only: NVAR, URHO, UTEMP, UEINT, UFS, UFX, UMX, UMZ, &
                                  cfl, do_hydro
#ifdef REACTIONS
    use meth_params_module, only: dtnuc_e, dtnuc_X, do_react
#endif
    use prob_params_module, only: dim
    use network, only: nspec, naux
    use extern_probin_module, only: small_x
    use eos_module

    use bl_fort_module, only : rt => c_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: so_lo(3), so_hi(3)
    integer          :: sn_lo(3), sn_hi(3)
#ifdef REACTIONS
    integer          :: ro_lo(3), ro_hi(3)
    integer          :: rn_lo(3), rn_hi(3)
#endif
    real(rt)         :: s_old(so_lo(1):so_hi(1),so_lo(2):so_hi(2),so_lo(3):so_hi(3),NVAR)
    real(rt)         :: s_new(sn_lo(1):sn_hi(1),sn_lo(2):sn_hi(2),sn_lo(3):sn_hi(3),NVAR)
#ifdef REACTIONS
    real(rt)         :: r_old(ro_lo(1):ro_hi(1),ro_lo(2):ro_hi(2),ro_lo(3):ro_hi(3),nspec+2)
    real(rt)         :: r_new(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3),nspec+2)
#endif
    real(rt)         :: dx(3), dt_old, dt_new

    integer          :: i, j, k
    real(rt)         :: rhooinv, rhoninv
#ifdef REACTIONS
    real(rt)         :: X_old(nspec), X_new(nspec), X_avg(nspec), X_dot(nspec)
    real(rt)         :: e_old, e_new, e_avg, e_dot
    real(rt)         :: tau_X, tau_e
#endif
    real(rt)         :: tau_CFL


    real(rt)         :: v(3), c
    type (eos_t)     :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhooinv = ONE / s_old(i,j,k,URHO)
             rhoninv = ONE / s_new(i,j,k,URHO)

             ! CFL hydrodynamic stability criterion

             ! If the timestep violated (v+c) * dt / dx > 1,
             ! suggest a new timestep such that (v+c) * dt / dx <= CFL,
             ! where CFL is the user's chosen timestep constraint.
             ! We don't use the CFL choice in the test for violation
             ! because if we did, then even a small increase in velocity
             ! over the timestep would be enough to trigger a retry
             ! (and empirically we have found that this can happen
             ! basically every timestep, which can greatly increase
             ! the cost of a simulation for not much benefit);
             ! but these types of issues are exactly why a user picks a
             ! safe buffer value like CFL = 0.8 or CFL = 0.5. We only
             ! want to trigger a retry if the timestep strongly violated
             ! the stability criterion.

             if (do_hydro .eq. 1) then

                eos_state % rho = s_new(i,j,k,URHO )
                eos_state % T   = s_new(i,j,k,UTEMP)
                eos_state % e   = s_new(i,j,k,UEINT) * rhoninv
                eos_state % xn  = s_new(i,j,k,UFS:UFS+nspec-1) * rhoninv
                eos_state % aux = s_new(i,j,k,UFX:UFX+naux-1) * rhoninv

                call eos(eos_input_re, eos_state)

                v = HALF * (s_old(i,j,k,UMX:UMZ) * rhooinv + s_new(i,j,k,UMX:UMZ) * rhoninv)

                c = eos_state % cs

                tau_CFL = minval(dx(1:dim) / (c + abs(v(1:dim))))

                if (dt_old > tau_CFL) then
                   dt_new = min(dt_new, cfl * tau_CFL)
                endif

             endif

#ifdef REACTIONS
             ! Burning stability criterion

             if (do_react .eq. 1) then

                X_old = s_old(i,j,k,UFS:UFS+nspec-1) * rhooinv
                X_new = s_new(i,j,k,UFS:UFS+nspec-1) * rhoninv
                X_avg = max(small_x, HALF * (X_old + X_new))
                X_dot = HALF * (r_old(i,j,k,1:nspec) + r_new(i,j,k,1:nspec))

                X_dot = max(abs(X_dot), 1.e-50_rt)
                tau_X = minval( X_avg / X_dot )

                e_old = s_old(i,j,k,UEINT) * rhooinv
                e_new = s_new(i,j,k,UEINT) * rhoninv
                e_avg = HALF * (e_old + e_new)
                e_dot = HALF * (r_old(i,j,k,nspec+1) + r_new(i,j,k,nspec+1))

                e_dot = max(abs(e_dot), 1.e-50_rt)
                tau_e = e_avg / e_dot

                if (dt_old > dtnuc_e * tau_e) then

                   dt_new = min(dt_new, dtnuc_e * tau_e)

                endif

                if (dt_old > dtnuc_X * tau_X) then

                   dt_new = min(dt_new, dtnuc_X * tau_X)

                endif

             endif
#endif

          enddo
       enddo
    enddo

  end subroutine ca_check_timestep

end module timestep_module
