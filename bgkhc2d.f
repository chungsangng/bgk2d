c----------------------------------------------------------------------2
      program bgkhc2d
c
c by C. S. Ng 
c
c Version: 2015/10/16
c
c  This program solve for a 2D BGK solution for the helical solution
c  of the Possioin-Laplace-Ampere equations.
c
c  connect with tridia solution for large radial distance
c
c  based on bgk2d19.f --- works for + and - potentials
c  
c
      integer nmax
      parameter (nmax = 1000000)
      real p(0:nmax),q(0:nmax)
      real az(0:nmax), bp(0:nmax), ap(0:nmax), bz(0:nmax)
c
c p is the electric potential (psi in the notes BGK-2D...SI-notes)
c q is dp/dr (r is the cylindrical radial coordinate rho)
c bp is the phi component of the magnetic field B
c az is the z component of the magnetic potential: bp = -d(az)/dr
c bz is the z component of the magnetic field B
c ap is the phi component of the magnetic potential: bz = d(ap)/dr+ap/r
c
      real a(nmax),b(nmax),c(nmax),f(nmax),t(nmax),t1(nmax)
      real c1,c2,c3,g1,p1,p2,ep1,f1,e1,g0
      real cc1,cc2,cc3
      real h0,k,b0,xi,az0,be2,rmax,r,dr,eps,fj,eps1
      real xi12,h1,kz2,r1,dr1
      real p0U,p0L,p0n,g,dqdr0,dbp0,p00,q00
      real pt1,pt2,qt1,qt2,rm
      real azt1,azt2,bpt1,bpt2
      real apt1,apt2,bzt1,bzt2
      real tsum,tdiff,ttol,tmax
      real rhs1,rhs2,rhs3
      real ne,vp,vz,te,tez,tmp,tmp1,e0
      real rpl,vmin,vmax,dv,fe,cfe,psi,azi,api,pi
      integer irpl,nv,nrho
      integer tcount,timax
      real perr,psum
      common/gpara/h0,b0,k,xi,az0,be2,xi12
      integer i,j,n,iout,ib,n1
      namelist/in/n,h0,k,b0,xi,az0,be2,rmax,eps,eps1,iout,fj,ttol,
     .    timax,tmax,rpl,vmin,vmax,nv
      read(*,in)
c
      nrho = n/iout/nv
c
      write(*,in)
      write(*,*)n
      write(*,*)iout
      write(*,*)h0
      write(*,*)k
      write(*,*)b0
      write(*,*)be2
      write(*,*)xi
      write(*,*)az0
      write(*,*)rpl
      write(*,*)vmin
      write(*,*)vmax
      write(*,*)nv
c
c h0 is a parameter that is < 1 (positive for electron holes, and 
c negative for negative potential solution --- magnetic holes)
c k and xi are positive parameters
c b0 is the z component of the magnetic field at r=0 
c az0 is the az at r=0
c be2 is the square of the ratio of electron thermal speed to c
c
c n is the number of r steps in the integration from r=0 to rmax
c dr is the fixed integration step in r
c
      if (h0.ge.1.) then
        write(*,*) "h0 >= 1"
        go to 1000
      endif
      if (h0.eq.0.) then
        write(*,*) "h0 = 0"
        go to 1000
      endif
      if (n.gt.nmax) then
        write(*,*) "n > nmax"
        go to 1000
      endif
      if (k.lt.0.) then
        write(*,*) "k < 0"
        go to 1000
      endif
      if (xi.lt.0.) then
        write(*,*) "xi < 0"
        go to 1000
      endif
c
      dr = rmax/n
      write(*,*)dr
      q(0) = 0.
      ap(0) = 0.
      bz(0) = b0
      bz(1) = bz(0)
      ap(1) = bz(0)*dr/2.
      az(0) = az0
      bp(0) = 0.
      xi12 = 1.+2.*xi
      h1 = 1.-h0*exp(-xi*az0**2/xi12)/sqrt(xi12)
      if (h0.gt.0.) then
        p0U = 1.01*log(1./h1)
        p00 = abs(p0U)
        p0L = 0.
      else
        p0L = 1.01*log(1./h1)
        p00 = abs(p0L)
        p0U = 0.
      endif
c
c      p0U = 1.01*log(1./h1)
c      p00 = p0U
c**
c**      write(*,*)"p0U= ", p0U
c**
c      p0L = 0.
c
c--begin loop 1
c
1     p0n = p0L + 0.5*(p0U-p0L) 
c**
c**        write(*,*)p0n
c**
      if ((abs(p0U-p0L)/abs(p0U+p0L).lt.eps).and.(ib.eq.1)) go to 100
c
c p0 converges but p still blows up before rho gets to the end of box
c
      do i = 1,n
        p(i) = 0.
        q(i) = 0.
      enddo
      dqdr0 = (exp(p0n)*h1 - 1.)/2.      
      p(0) = p0n
      p00 = abs(p(0))
      q(1) = dqdr0*dr
      q00 = abs(dqdr0)
      p(1) = p(0) + 0.5*q(1)*dr
      kz2 = 2.*be2*h0*xi*exp(p0n-xi*az0**2/xi12)/xi12**(1.5)
      dbp0 = kz2*az0/2.
      bp(1) = dbp0*dr
      az(1) = az0 - 0.5*bp(1)*dr
c
c--begin loop i -- integrate along rho
c
      do i = 2,n
        r = i*dr
        dr1 = 1. + dr/r
        pt1 = p(i-1) + q(i-1)*dr
        azt1 = az(i-1) - bp(i-1)*dr
        apt1 = (ap(i-1)+bz(i-1)*dr)/dr1
        qt1 = (q(i-1) + rhs1(r,pt1,azt1,apt1)*dr)/dr1
        bpt1 = (bp(i-1) + rhs2(r,pt1,azt1,apt1)*dr)/dr1
        bzt1 = bz(i-1) + rhs3(r,pt1,azt1,apt1)*dr
c
c--begin loop 10 -- iterate to make it semi-implicit
c    This is to make sure blow-up is not numerical
c
        tcount = 0
10      tcount = tcount + 1
        tsum = pt1**2+azt1**2+apt1**2+qt1**2+bpt1**2+bzt1**2
        pt2 = p(i-1) + qt1*dr
        azt2 = az(i-1) - bpt1*dr
        apt2 = (ap(i-1)+bzt1*dr)/dr1
        qt2 = (q(i-1) + rhs1(r,pt2,azt2,apt2)*dr)/dr1
        bpt2 = (bp(i-1) + rhs2(r,pt2,azt2,apt2)*dr)/dr1
        bzt2 = bz(i-1) + rhs3(r,pt2,azt2,apt2)*dr
        tdiff = (pt2-pt1)**2+(azt2-azt1)**2+(apt2-apt1)**2
     .    +(qt2-qt1)**2+(bpt2-bpt1)**2+(bzt2-bzt1)**2
        if (tcount.gt.timax) then 
c          write(*,*) "iteration > timax",i,tsum,tdiff
          go to 20
c
c cannot converge => use current values but printout a warning
c
        endif
c
        if ((abs(pt2-pt1).lt.ttol*abs(pt1)).or.((abs(pt1).lt.ttol).and.
     .     (abs(pt2).lt.ttol))) then
        if ((abs(azt2-azt1).lt.ttol*abs(azt1)).or.((abs(azt1).lt.ttol)
     .     .and.(abs(azt2).lt.ttol))) then
        if ((abs(apt2-apt1).lt.ttol*abs(apt1)).or.((abs(apt1).lt.ttol)
     .     .and.(abs(apt2).lt.ttol))) then
        if ((abs(qt2-qt1).lt.ttol*abs(qt1)).or.((abs(qt1).lt.ttol).and.
     .     (abs(qt2).lt.ttol))) then
        if ((abs(bpt2-bpt1).lt.ttol*abs(bpt1)).or.((abs(bpt1).lt.ttol)
     .     .and.(abs(bpt2).lt.ttol))) then
        if ((abs(bzt2-bzt1).lt.ttol*abs(bzt1)).or.((abs(bzt1).lt.ttol)
     .     .and.(abs(bzt2).lt.ttol))) go to 20
c
        endif
        endif
        endif
        endif
        endif
c
        pt1 = pt2
        azt1 = azt2
        apt1 = apt2
        qt1 = qt2
        bpt1 = bpt2
        bzt1 = bzt2
        go to 10
c
c--end loop 10
c
20      p(i) = pt2
        az(i) = azt2
        ap(i) = apt2
        q(i) = qt2
c        if (q(i).lt.q00) then
c          q00 = q(i)
c          write(*,*) "q00= ", q00
c        endif
        bp(i) = bpt2
        bz(i) = bzt2
          
c
        if (h0.gt.0.) then
          if (q(i).gt.tmax*q00) then
            p0U = p0n
            j = i
            ib = 1
            go to 1
c
c--end loop 1 --- p blows up positive => reduce p0
c
          endif
          if (p(i).lt.-tmax*p00) then
            p0L = p0n
            j = i
            ib = 2
            go to 1
c
c--end loop 1 --- p blows up negative => increase p0
c
          endif
        else
          if (q(i).lt.-tmax*q00) then
            p0L = p0n
            j = i
            ib = 1
            go to 1
c
c--end loop 1 --- p blows up positive => reduce p0
c
          endif
          if (p(i).gt.tmax*p00) then
            p0U = p0n
            j = i
            ib = 2
            go to 1
c
c--end loop 1 --- p blows up negative => increase p0
c
          endif
        endif
      enddo
c
c--end loop i --- get to here if p doesn't blow up either way
c
      j = n
c
c--begin loop 100 -- connect to large rho solution
c
100   j = fj*j
c
c backup some distance in rho before connection
c
      n1 = n - j
      cc3 = ap(j)/(j*dr)
c
c A_phi -- oo --> cc3 * rho
c
      cc1 = (az(j+1)-az(j))/(log((j+1)*dr)-log(j*dr))
      cc2 = az(j) - cc1*log(j*dr)
c
c A_z -- oo --> cc1 * log(rho) + cc2 
c
      do i = 1,n1
        ap(j+i) = cc3*(j+i)*dr
        az(j+i) = cc1*log((j+i)*dr) + cc2
      enddo
c
c Set the zero_th order A_phi and A_z
c
      c1 = 1./dr/dr
      c2 = 2.*c1
      c3 = 1./2./dr
      p1 = p(j)
      p2 = 1.-p1
      ep1 = exp(p1)
      r = (j+1)*dr
      g1 = ep1*g(r,az(j+1),ap(j+1))
      f(1) = g1*p2 - 1. - (c1 -c3/r)*p(j)
      b(1) = -(c2 + g1)
      c(1) = c1 + c3/r
      do i = 2,n1-1
        r = (j+i)*dr
        g0 = g(r,az(j+i),ap(j+i))
        p1 = log(1./g0)
        p2 = 1. - p1
        ep1 = exp(p1)
        g1 = ep1*g0
        a(i) = c1 - c3/r
        b(i) = -(c2 + g1)
        c(i) = c1 + c3/r
        f(i) = g1*p2 - 1.
      enddo
      r = n*dr
      g0 = g(r,az(n),ap(n))
      p1 = log(1./g0)
      p2 = 1. - p1
      ep1 = exp(p1)
      g1 = ep1*g0
      f(n1) = g1*p2 - 1.
      a(n1) = c1 - c3/r
      if (g0.ne.1.) then
        f1 = log(1./g((n+1)*dr,cc1*log((n+1)*dr) + cc2,
     .       cc3*(n+1)*dr))/log(1./g0)
        if (f1.gt.exp(-dr)) then
          e1 = f1
        else
          e1 = exp(-dr)
        endif
      else
        e1 = exp(-dr)
      endif
      b(n1) = (c1 + c3/r)*e1 - (c2 + g1)
      call tridia(a,b,c,f,n1,t)
c
c--begin loop 200
c
c 
      tcount = 0
200   r = j*dr
      tcount = tcount + 1
      az(j+1) = az(j) - bp(j)*dr
c
c bp = -d(az)/dr
c
      ap(j+1) = ap(j)*(1.-dr/r) + bz(j)*dr
c
c bz = d(ap)/dr+ap/r
c
      bp(j+1) = bp(j)*(1.-dr/r) + rhs2(r,p(j),az(j),ap(j))*dr
      bz(j+1) = bz(j) + rhs3(r,p(j),az(j),ap(j))*dr
      do i = 2,n1
        r = (j+i-1)*dr
        az(j+i) = az(j+i-1) - bp(j+i-1)*dr
        ap(j+i) = ap(j+i-1)*(1.-dr/r) + bz(j+i-1)*dr
        bp(j+i) = bp(j+i-1)*(1.-dr/r) + rhs2(r,t(i-1),
     .            az(j+i-1),ap(j+i-1))*dr
        bz(j+i) = bz(j+i-1) + rhs3(r,t(i-1),
     .            az(j+i-1),ap(j+i-1))*dr
      enddo
      cc3 = ap(n)/(n*dr)
      cc1 = (az(n)-az(n-1))/(log(n*dr)-log((n-1)*dr))
      cc2 = az(n) - cc1*log(n*dr)
c       
      do i = 1,n1
        t1(i) = t(i)
      enddo
      p1 = t1(1)
      p2 = 1.-p1
      ep1 = exp(p1)
      r = (j+1)*dr
      g1 = ep1*g(r,az(j+1),ap(j+1))
      f(1) = g1*p2 - 1. - (c1 -c3/r)*p(j)
      b(1) = -(c2 + g1)
      c(1) = c1 + c3/r
      do i = 2,n1-1
        r = (j+i)*dr
        g0 = g(r,az(j+i),ap(j+i))
        p1 = t1(i)
        p2 = 1. - p1
        ep1 = exp(p1)
        g1 = ep1*g0
        a(i) = c1 - c3/r
        b(i) = -(c2 + g1)
        c(i) = c1 + c3/r
        f(i) = g1*p2 - 1.
      enddo
      r = n*dr
      g0 = g(r,az(n),ap(n))
      p1 = t1(n1)
      p2 = 1. - p1
      ep1 = exp(p1)
      g1 = ep1*g0
      f(n1) = g1*p2 - 1.
      a(n1) = c1 - c3/r
      if (g0.ne.1.) then
        f1 = log(1./g((n+1)*dr,cc1*log((n+1)*dr) + cc2,
     .       cc3*(n+1)*dr))/log(1./g0)
        if (f1.gt.exp(-dr)) then
          e1 = f1
        else
          e1 = exp(-dr)
        endif
      else
        e1 = exp(-dr)
      endif
      b(n1) = (c1 + c3/r)*e1 - (c2 + g1)
      call tridia(a,b,c,f,n1,t)
      perr = 0.
      psum = 0.
      do i = 1,n1
        perr = perr + abs(t(i)-t1(i))
        psum = psum + abs(t(i))
      enddo
c**      write(*,*)"perr= ",perr," psum= ",psum,tcount
      if ((perr/psum.gt.eps1).and.(tcount.le.timax)) go to 200
      write(*,*)"perr= ",perr," psum= ",psum
c
c--end loop 200
c      
      do i = 1,n1
        p(j+i) = t(i)
      enddo
      do i = 1,n1-1
        q(j+i) = c3*(p(j+i+1)-p(j+i-1))
      enddo
      q(n) = (p(n)-p(n-1))/dr
c
      write(*,*)n/iout
      write(*,*)"p0U = ",p0U
      write(*,*)"p0L = ",p0L
      write(*,*)"err = ",abs(p0U-p0L)/abs(p0U+p0L)
      write(*,500)"   r               ",
     .            "   p               ",
     .            "   q               ",
     .            "   g               ",
     .            "   rho             ",
     .            "   az              ",
     .            "   bp              ",
     .            "   ap              ",
     .            "   bz              ",
     .            "   ne              ",
     .            "   vp              ",
     .            "   vz              ",
     .            "   te              ",
     .            "   tez             "
c      write(*,500)"   r                     ",
c     .            "   p                     ",
c     .            "   q                     ",
c     .            "   g                     ",
c     .            "   rho                   ",
c     .            "   az                    ",
c     .            "   bp                    ",
c     .            "   ap                    ",
c     .            "   bz                    "
      r = 0.
      g0 = g(r,az(0),ap(0))
      ne = exp(p(0))*g0
      tmp1 = 4.*k*r*r
      tmp = 1. + 2.*tmp1
      e0 = exp(-tmp1*ap(0)*ap(0)/tmp-xi*az(0)*az(0)/xi12)
      vp = 0.
      vz = -2.*h0*xi*az(0)*e0/xi12/(sqrt(tmp*xi12)-h0*e0)
      te = (exp(p(0))*(3.-h0*(1.+(tmp+(2.*tmp1*ap(0))**2)/tmp**2
     .  + (xi12+(2.*xi*az(0))**2)/xi12**2)*e0/sqrt(tmp*xi12))/ne 
     .  - vp*vp - vz*vz)/3.
      tez = exp(p(0))*(1.-h0*((xi12+(2.*xi*az(0))**2)/xi12**2)
     .  *e0/sqrt(tmp*xi12))/ne - vz*vz
      write(*,600) 0.,p(0),0.,g0,1.-ne,
     .   az(0),bp(0),ap(0),bz(0),ne,vp,vz,te,tez
      do i = iout,n,iout
        r = i*dr
        g0 = g(r,az(i),ap(i))
        ne = exp(p(i))*g0
        tmp1 = 4.*k*r*r
        tmp = 1. + 2.*tmp1
        e0 = exp(-tmp1*ap(i)*ap(i)/tmp-xi*az(i)*az(i)/xi12)
        vp = -2.*h0*tmp1*ap(i)*e0/tmp/(sqrt(tmp*xi12)-h0*e0)
        vz = -2.*h0*xi*az(i)*e0/xi12/(sqrt(tmp*xi12)-h0*e0)
        te = (exp(p(i))*(3.-h0*(1.+(tmp+(2.*tmp1*ap(i))**2)/tmp**2
     .    + (xi12+(2.*xi*az(i))**2)/xi12**2)*e0/sqrt(tmp*xi12))/ne 
     .    - vp*vp - vz*vz)/3.
        tez = exp(p(i))*(1.-h0*((xi12+(2.*xi*az(i))**2)/xi12**2)
     .    *e0/sqrt(tmp*xi12))/ne - vz*vz
        write(*,600)i*dr,p(i),q(i),g0,1.-ne,
     .    az(i),bp(i),ap(i),bz(i),ne,vp,vz,te,tez
      enddo
      write(*,*)
      dv = (vmax-vmin)/nv
      psi = p(0)
      azi = az(0)
      pi = 4.*atan(1.)
      cfe = (2.*pi)**(-1.5)
      do j=0,nv
        vz = vmin + j*dv
        do i=0,nv
          vp = vmin + i*dv
          fe = cfe*exp(-(vp*vp+vz*vz)/2.+psi)
     .         *(1.-h0*exp(-xi*(vz-azi)**2))
          write(*,700)log(fe)
        enddo
      enddo
      write(*,*)
      irpl = rpl/dr
      rpl = irpl*dr
      psi = p(irpl)
      azi = az(irpl)
      api = ap(irpl)
      do j=0,nv
        vz = vmin + j*dv
        do i=0,nv
          vp = vmin + i*dv
          fe = cfe*exp(-(vp*vp+vz*vz)/2.+psi)
     .         *(1.-h0*exp(-k*(2.*rpl*(vp-api))**2-xi*(vz-azi)**2))
          write(*,700)log(fe)
        enddo
      enddo
      write(*,*)
      do j=0,nv
        vz = vmin + j*dv
        vp = 0.
        do i=0,nv
          r = i*iout*nrho*dr
          psi = p(i*iout*nrho)
          azi = az(i*iout*nrho)
          api = ap(i*iout*nrho)
          fe = cfe*exp(-(vp*vp+vz*vz)/2.+psi)
     .         *(1.-h0*exp(-k*(2.*r*(vp-api))**2-xi*(vz-azi)**2))
          write(*,700)log(fe)
        enddo
      enddo
      write(*,*)
      do j=0,nv
        vz = 0.
        vp = vmin + j*dv
        do i=0,nv
          r = i*iout*nrho*dr
          psi = p(i*iout*nrho)
          azi = az(i*iout*nrho)
          api = ap(i*iout*nrho)
          fe = cfe*exp(-(vp*vp+vz*vz)/2.+psi)
     .         *(1.-h0*exp(-k*(2.*r*(vp-api))**2-xi*(vz-azi)**2))
          write(*,700)log(fe)
        enddo
      enddo
c          
500   format(14(a19))
600   format(14(e18.8,x))
700   format(e18.8)
c500   format(9(a25))
c600   format(9(e24.14,x))
1000  stop
      end
c----------------------------------------------------------------------2
      function g(r,az,ap)
      real g
      real r,az,ap
      real tmp,xi12,tmp1
      real h0,b0,k,xi,az0,be2
      common/gpara/h0,b0,k,xi,az0,be2,xi12
c
      tmp1 = 4.*k*r*r
      tmp = 1. + 2.*tmp1
      xi12 = 1. + 2.*xi
      g = 1. - h0*exp(-tmp1*ap*ap/tmp-xi*az*az/xi12)/sqrt(tmp*xi12)
      return
      end
c----------------------------------------------------------------------2
      function rhs1(r,p,az,ap)
      real rhs1
      real r,p,az,ap
      real h0,b0,k,xi,az0,be2,xi12
      real g
      common/gpara/h0,b0,k,xi,az0,be2,xi12
c
      rhs1 = exp(p)*g(r,az,ap) - 1.
      return
      end
c----------------------------------------------------------------------2
      function rhs2(r,p,az,ap)
      real rhs2
      real r,p,az,ap
      real tmp,xi12,tmp1
      real h0,b0,k,xi,az0,be2
      common/gpara/h0,b0,k,xi,az0,be2,xi12
c
      tmp1 = 4.*k*r*r
      tmp = 1. + 2.*tmp1
      xi12 = 1. + 2.*xi
      rhs2 = 2.*be2*h0*xi*az*exp(p-tmp1*ap*ap/tmp-xi*az*az/xi12)/
     .      sqrt(tmp*xi12**3)
      return
      end
c----------------------------------------------------------------------2
      function rhs3(r,p,az,ap)
      real rhs3
      real r,p,az,ap
      real tmp,xi12,tmp1
      real h0,b0,k,xi,az0,be2
      common/gpara/h0,b0,k,xi,az0,be2,xi12
c
      tmp1 = 4.*k*r*r
      tmp = 1. + 2.*tmp1
      xi12 = 1. + 2.*xi
      rhs3 = -2.*be2*h0*tmp1*ap*exp(p-tmp1*ap*ap/tmp-xi*az*az/xi12)/
     .      sqrt(xi12*tmp**3)
      return
      end
c----------------------------------------------------------------------2
c tridiagonal matrix algorithm
c
c copied from Y. Jaluria K. Torrance ,computational heat transfer pp336
c
c a,b and c are the three elements in each row. With b at the diagonal
c f is the constant on the right-hand side of each equation. n is the
c number of equations and t is the variable to be computed.
c
      subroutine tridia(a,b,c,f,n,t)
c
      integer nmax
      parameter (nmax = 1000000)
      integer m,n,i,j
      parameter (m=100)
c
c      real a(2:m-1),b(m-1),c(m-2),f(m-1),t(m-1),d
      real a(nmax),b(nmax),c(nmax),f(nmax),t(nmax),d
c
c Elimiante the a's by Gaussian elimination and determine
c the new coefficients.
c
      do i = 2,n
        d = a(i)/b(i-1)
        b(i) = b(i) - c(i-1)*d
        f(i) = f(i) - f(i-1)*d
      enddo
c
c back subsititution
c
      t(n) = f(n)/b(n)
      do i = 1,n-1
        j = n-i
        t(j) = (f(j) - c(j)*t(j+1))/b(j)
      enddo
      return
      end
c----------------------------------------------------------------------2
