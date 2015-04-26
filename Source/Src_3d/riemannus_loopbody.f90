
          rl = max(ql(i,j,kc,QRHO),small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ),small_pres)
          rel =     ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO),small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

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
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
             v1gdnv = v1r
             v2gdnv = v2r
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
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
             if ( (i.eq.domlo(1)   .and. zerov_lo) .or. &
                  (i.eq.domhi(1)+1 .and. zerov_hi) ) then
                zerov_fac = ZERO
             else
                zerov_fac = ONE
             end if
          end if

          ugdnv(i,j,kc) = ugdnv(i,j,kc) * zerov_fac

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*v1gdnv
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*v2gdnv

          rhoetot = regdnv + HALF*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))
          uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv

          rgd1d(i) = rgdnv
          us1d(i) = ustar

