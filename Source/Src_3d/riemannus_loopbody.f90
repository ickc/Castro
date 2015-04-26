
          rl = max(ql(i,j,kc,QRHO),small_dens)

          ! pick left velocities based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,kc,QU)
             v1l = ql(i,j,kc,QV)
             v2l = ql(i,j,kc,QW)
          elseif(idir.eq.2) then
             ul  = ql(i,j,kc,QV)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QW)
          else
             ul  = ql(i,j,kc,QW)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QV)
          endif

          pl  = max(ql(i,j,kc,QPRES ),small_pres)
          rel =     ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO),small_dens)

          ! pick right velocities based on direction
          if(idir.eq.1) then
             ur  = qr(i,j,kc,QU)
             v1r = qr(i,j,kc,QV)
             v2r = qr(i,j,kc,QW)
          elseif(idir.eq.2) then
             ur  = qr(i,j,kc,QV)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QW)
          else
             ur  = qr(i,j,kc,QW)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QV)
          endif

          pr  = max(qr(i,j,kc,QPRES),small_pres)
          rer =     qr(i,j,kc,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
          pstar = max(pstar,small_pres)

          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl(i,j)
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
          endif
          ro = max(small_dens,ro)

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          entho = (reo/ro + po/ro)/co**2
          rstar = ro + (pstar - po)/co**2
          rstar = max(small_dens,rstar)
          estar = reo + (pstar - po)*entho
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)
          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar .gt. ZERO) then
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             v1gdnv = v1r
             v2gdnv = v2r
          else
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
          endif
          rgdnv = frac*rstar + (ONE - frac)*ro

          ugdnv(i,j,kc) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j,kc) = frac*pstar + (ONE - frac)*po

          regdnv = frac*estar + (ONE - frac)*reo
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j,kc) = uo
             pgdnv(i,j,kc) = po
             regdnv = reo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j,kc) = ustar
             pgdnv(i,j,kc) = pstar
             regdnv = estar
          endif

          gegdnv(i,j,kc) = pgdnv(i,j,kc)/regdnv + 1.0d0

          pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                  (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                  (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                  (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                  (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 3) then
             if (k3d.eq.domlo(3) .and. &
                  (physbc_lo(3) .eq. Symmetry .or.  physbc_lo(3) .eq. SlipWall .or. &
                  physbc_lo(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (k3d.eq.domhi(3)+1 .and. &
                  (physbc_hi(3) .eq. Symmetry .or.  physbc_hi(3) .eq. SlipWall .or. &
                  physbc_hi(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

          if(idir.eq.1) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          elseif(idir.eq.2) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          else
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          endif

          rhoetot = regdnv + HALF*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))
          uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv

#ifdef SGS
          ! Treat K as a passively advected quantity but allow it to affect fluxes of (rho E) and momenta.
          if (UESGS .gt. -1) then
        xxx     n  = UESGS
             nq = QESGS
             if (ustar .gt. ZERO) then
                qavg = ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                qavg = qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
             endif

             uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg

             rho_K_contrib =  TWO3RD * rgdnv * qavg

             if(idir.eq.1) then
                uflx(i,j,kflux,UMX) = uflx(i,j,kflux,UMX) + rho_K_contrib
             elseif(idir.eq.2) then
                uflx(i,j,kflux,UMY) = uflx(i,j,kflux,UMY) + rho_K_contrib
             elseif(idir.eq.3) then
                uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,UMZ) + rho_K_contrib
             endif

             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
          end if
#endif

          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
