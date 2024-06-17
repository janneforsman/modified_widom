      implicit double precision (a-h,o-z)      
      include 't.inc'    
      pi = acos(-1.d0)
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      elch = 1.602D-19
      avno = 6.02214D23
      bk = 1.38066D-23
      faraday = 8.85418782D-12
      dielc = 78.3d0
      temp = 298.d0
      mmm= 47
      lll= 49
      kkk = 45
      isl = 30
      open (mmm,file='input',form='formatted')
      open (kkk,file='coordm',form='formatted')
      open (lll,file='coords',form='formatted')      
      open (isl,file='slump',form='formatted')
      rewind isl
      read(isl,*) islu
      rewind mmm
      read(mmm,*) 
      read(mmm,*) Nsub
      read(mmm,*) 
      read(mmm,*) Nm
      Ns = Nm
      write(*,*) 'NOTE: Ns = Nm',Nm,Ns
      read(mmm,*) 
      read(mmm,*) ed2
      read(mmm,*) 
      read(mmm,*) dielc
      read(mmm,*) 
      read(mmm,*) dhspp,dhsnn
      read(mmm,*) 
      read(mmm,*) bltrams
      read(mmm,*) 
      read(mmm,*) kread,knsamp
      read(mmm,*) 
      read(mmm,*) dr
      read(mmm,*) 
      read(mmm,*) dpsp,dpsn
      read(mmm,*) 
      read(mmm,*) fracclm,Rclmax,cldps
      write(*,*) 'OBS.! Varying cluster radius, Rclmax = ',Rclmax          
      write(*,*) 'cldps = ',cldps
      write(*,*) 'fracclm = ',fracclm                  
      
      bjerrum = elch*elch/(fourpi*bk*temp*faraday*dielc*1.e-10)
      write(*,*) 'bjerrum/AA = ',bjerrum
      ed = 0.5d0*ed2
      edsq = ed*ed
      red = 1.d0/ed
      red2 = 1.d0/ed2
c      dhs2 = dhs**2
      dhspn = 0.5d0*(dhspp+dhsnn)      
      dhspp2 = dhspp**2
      dhspn2 = dhspn**2
      dhsnn2 = dhsnn**2
c      shuss = 2.d0/ed
      shuss = 7.d0/(4.d0*ed)
      y5 = 21.d0/(4.d0*ed**5)
      y6 = 7.d0/ed**6
      y7 = 2.5d0/ed**7
      redsq = 1.d0/edsq
      write(*,*) 'ed,shuss = ',ed,shuss
      write(*,*) 'centrally placed charges!!! '   
      write(*,*) 'dhspp,dhsnn = ',dhspp,dhsnn
      write(*,*) 'dhspn = ',dhspn      
      write(*,*) 'Nm = ',Nm
      write(*,*) 'Nsub = ',Nsub
      write(*,*) 'lateral length of box (ed2): ',ed2
      write(*,*) 'displacement parameter, dpsp = ',dpsp
      write(*,*) 'displacement parameter, dpsn = ',dpsn      
      write(*,*) 'knsamp (rdf distr.) = ',knsamp
      write(*,*) 'dr (rdf distr.) = ',dr
      rNm = dfloat(Nm)
      rNs = rNm
      etam = pi*dhspp**3/6.d0
      etas = pi*dhsnn**3/6.d0      
      eta = (etam+etas)*rNm/ed2**3
      write(*,*) 'volume fraction (linear) = ',eta      
      if (kread.eq.1) then
      rewind kkk
      rewind lll
      do i = 1,Nm
      read(kkk,*) xhp(i),yhp(i),zhp(i)
      read(lll,*) xhn(i),yhn(i),zhn(i)      
      enddo
      elseif (kread.eq.2) then
      do i = 1,Nm
      read(kkk,*) xhp(i),yhp(i),zhp(i)
      read(kkk,*) t1,t2,t3
      read(lll,*) xhn(i),yhn(i),zhn(i)
      read(lll,*) t1,t2,t3
      enddo
      else
      CALL cinit
      endif
      Nmac = 10
      rNm = dfloat(Nm)
      rrNm = 1./rNm
      rdr = 1.d0/dr
      CALL ucalc
      utot = bjerrum*utot
      write(*,*) 'utot(init) = ',utot
      irdf = 70
      open(irdf,file='rdf',form='formatted')
      rewind irdf
      rnss = 0.
      rnmm = 0.
      rnms = 0.
      asampnz = 0.
      do j = 1,maxnz
      grhphp(i) = 0.d0
      grhphn(i) = 0.d0
      grhnhn(i) = 0.d0
      enddo
c      ia =  90
c      open (ia,file='a.pqr',form='formatted')
c      rewind ia
c      do i = 1,Nm
c      write(ia,'(A,i7,A,i6,1f12.2,2f9.2,2f6.2)')
c     *'ATOM',i,' unk unk ',i,xhp(i),yhp(i),zhp(i),0.d0,2.d0
c     *'ATOM',i,' uns uns ',i,xhn(i),yhn(i),zhn(i),0.d0,2.d0
c      enddo
c      stop
      rhphn = 0.
      attdp = 0.d0
      accdp = 0.d0
      rejdp = 0.d0
      attdn = 0.d0
      accdn = 0.d0
      rejdn = 0.d0
      attclp = 0.
      attcln = 0.
      accclp = 0.
      acccln = 0.
      rejclp = 0.
      rejcln = 0.      
      do kmac = 1,Nmac
      do ksub = 1,Nsub         
      ks = int(1.d0+rNm*ran2(islu))
c      ks = 40
      if (ks.gt.Nm) ks = Nm
      if (ks.lt.1) ks = 1
      if (ran2(islu).lt.0.5d0) then
c     attempted displacement of cation

      if (ran2(islu).lt.fracclm) then
c     attempted cluster displacement of cation
      Rclust = dhspp+ran2(islu)*(Rclmax-dhspp)
      Rclust2 = Rclust*Rclust
      attclp = attclp+1.
c      write(*,*) 'attclp = ',attclp
      call clp
      if (idonk.eq.1) goto 927
      dut = bjerrum*(udn-udo) 
      if (dut.gt.100.d0) goto 927      
      if (dut.lt.0.d0) goto 9162
      if (dexp(-dut).lt.ran2(islu)) goto 927
 9162 accclp = accclp+1.      
      utot = utot+dut
      xhp(ks) = txs
      yhp(ks) = tys
      zhp(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xhn(km) = xcn(km)
      yhn(km) = ycn(km)
      zhn(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xhp(km) = xcp(km)
      yhp(km) = ycp(km)
      zhp(km) = zcp(km)
      enddo
      goto 6
 927  rejclp = rejclp+1.
      goto 6
      else         
c     attempted ordinary displacement of cation                  
      attdp = attdp+1.
      call dispp
      if (idonk.eq.1) goto 203
      delta = bjerrum*(un-uo) 
      if (delta.lt.0.) goto 527
      if (dexp(-delta).lt.ran2(islu)) go to 203
 527  accdp = accdp+1.
      xhp(ks) = tx
      yhp(ks) = ty
      zhp(ks) = tz
      utot = utot+delta
      goto 6
 203  rejdp = rejdp+1.
      goto 6
      endif
      else
c     attempted displacement of anion         

      if (ran2(islu).lt.fracclm) then
c     attempted cluster displacement of anion
      Rclust = dhspp+ran2(islu)*(Rclmax-dhspp)
      Rclust2 = Rclust*Rclust
      attcln = attcln+1.         
      call cln

      if (idonk.eq.1) goto 727
      dut = bjerrum*(udn-udo) 
      if (dut.gt.100.d0) goto 727      
      if (dut.lt.0.d0) goto 9762
      if (dexp(-dut).lt.ran2(islu)) goto 727
 9762 acccln = acccln+1.      
      utot = utot+dut
      xhn(ks) = txs
      yhn(ks) = tys
      zhn(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xhn(km) = xcn(km)
      yhn(km) = ycn(km)
      zhn(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xhp(km) = xcp(km)
      yhp(km) = ycp(km)
      zhp(km) = zcp(km)
      enddo
      goto 6
 727  rejcln = rejcln+1.
      goto 6
      else
c     attempted ordinary displacement of anion                           
      attdn = attdn+1.
      call dispn
      if (idonk.eq.1) goto 403
      delta = bjerrum*(un-uo) 
      if (delta.lt.0.) goto 327
      if (dexp(-delta).lt.ran2(islu)) go to 403
 327  accdn = accdn+1.
      xhn(ks) = tx
      yhn(ks) = ty
      zhn(ks) = tz
      utot = utot+delta
      goto 6
 403  rejdn = rejdn+1.
      goto 6
      endif
      endif
 6    continue

      if (mod(ksub,knsamp).eq.0) then 
      samp = samp+1.
      do i = 1,Nm
      tx = xhp(i)
      ty = yhp(i)
      tz = zhp(i)
      do k = i+1,Nm
      rr = dist(tx,ty,tz,xhp(k),yhp(k),zhp(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      grhphp(idr) = grhphp(idr)+1.
      endif
      enddo
      do k = 1,Nm
      rr = dist(tx,ty,tz,xhn(k),yhn(k),zhn(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      grhphn(idr) = grhphn(idr)+1.
      rhphn = rhphn+1.
      endif
      enddo
      enddo
      do i = 1,Nm
      tx = xhn(i)
      ty = yhn(i)
      tz = zhn(i)
      do k = i+1,Nm
      rr = dist(tx,ty,tz,xhn(k),yhn(k),zhn(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      grhnhn(idr) = grhnhn(idr)+1.
      endif
      enddo
      enddo
      endif
      
      enddo
      write(*,*) 
      write(*,*) kmac
      sutot = utot
      CALL ucalc
      utot = bjerrum*utot
      duu = utot-sutot
      write(*,*) utot,sutot
      write(*,*) 'utot,diff in utot (tot,rel.) = ',utot,duu,duu/utot
      write(*,*) 'faccdp,faccdn = ',accdp/attdp,accdn/attdn
      write(*,*) 'faccclp = ',accclp/attclp
      write(*,*) 'facccln = ',acccln/attcln
      write(*,*) attclp,accclp,rejclp
      write(*,*) attcln,acccln,rejcln      
      enddo
      write(*,*) 'samp = ',samp
      rewind kkk
      rewind lll
      do i = 1,Nm
      write(kkk,*) xhp(i),yhp(i),zhp(i)
      write(lll,*) xhn(i),yhn(i),zhn(i)
      enddo
      
      dV = ed2*ed2*ed2
      write(*,*) 'samp = ',samp
      nr = int(ed/dr)      
      avNm = dfloat(Nm)
      avNs = dfloat(Ns)
      refmm = 0.5d0*avNm*(avNm-1.)*samp/ed2**3
      refms = avNm*avNs*samp/ed2**3
      refss = 0.5d0*avNs*(avNs-1.)*samp/ed2**3
      write(*,*) grhphp(1),grhphp(2),grhphp(3)
      r = -dr*0.5d0
      rewind irdf
      do j = 1,nr
      r = r+dr
      rup = r+0.5d0*dr
      rlow = r-0.5d0*dr
      dV = 4.d0*pi*(rup**3-rlow**3)/3.d0
      refnomm = refmm*dV
      refnoms = refms*dV
      refnoss = refss*dV      
      rdfmm = grhphp(j)/refnomm
      rdfms = grhphn(j)/refnoms
      rdfss = grhnhn(j)/refnoss
      write(irdf,'(5e16.7)') r,rdfmm,rdfms,rdfss
      enddo

      rewind isl
      write(isl,*) islu
      STOP                                                              
      END

      subroutine cln
      implicit double precision (a-h,o-z)      
      include 't.inc'                              
      idonk = 0
      dpx = cldps*(ran2(islu)-0.5d0)
      dpy = cldps*(ran2(islu)-0.5d0)
      dpz = cldps*(ran2(islu)-0.5d0)
      txs = xhn(ks)+dpx
      tys = yhn(ks)+dpy
      tzs = zhn(ks)+dpz
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif
      if (tzs.gt.ed) then
      tzs = tzs-ed2
      elseif (tzs.lt.-ed) then
      tzs = tzs+ed2
      endif
      kclp = 0
      ttxs = xhn(ks)
      ttys = yhn(ks)
      ttzs = zhn(ks)
      udo = 0.d0
      do k = 1,Nm
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz
      else
      udo = -ucr(rmm2)+udo               
      endif
      enddo
      
      kcln = 0
      do k = 1,ks-1
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz

      else
      udo = ucr(rmm2)+udo         
      endif
      enddo

      do k = ks+1,Ns
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz

      else
      udo = ucr(rmm2)+udo         
      endif
      enddo

      udn = 0.d0
      kkclp = 0
      do k = 1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      udn = -ucr(rmm2)+udn         
      endif   
      enddo

      kkcln = 0
      do k = 1,ks-1
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      udn = ucr(rmm2)+udn         
      endif   
      enddo
      do k = ks+1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      udn = ucr(rmm2)+udn         
      endif
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcn(ks) = txs
      ycn(ks) = tys
      zcn(ks) = tzs      
      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xhp(km)
      tty = yhp(km)
      ttz = zhp(km)
      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo            
      enddo

      do k = 1,Ns      
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      udn = -ucr(rmm2)+udn
      udo = -uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo       
      enddo
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xhn(km)
      tty = yhn(km)
      ttz = zhn(km)      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo            
      enddo

      do k = 1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      udn = -ucr(rmm2)+udn
      udo = -uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo       
      enddo
      enddo
      goto 929
 927  idonk = 1
 929  continue
      return
      end                  

      subroutine clp
      implicit double precision (a-h,o-z)      
      include 't.inc'                              
      idonk = 0
      dpx = cldps*(ran2(islu)-0.5d0)
      dpy = cldps*(ran2(islu)-0.5d0)
      dpz = cldps*(ran2(islu)-0.5d0)
      txs = xhp(ks)+dpx
      tys = yhp(ks)+dpy
      tzs = zhp(ks)+dpz
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif
      if (tzs.gt.ed) then
      tzs = tzs-ed2
      elseif (tzs.lt.-ed) then
      tzs = tzs+ed2
      endif
      kcln = 0
      ttxs = xhp(ks)
      ttys = yhp(ks)
      ttzs = zhp(ks)
      udo = 0.d0
      do k = 1,Ns
      xcn(k) = xhn(k)
      ycn(k) = yhn(k)
      zcn(k) = zhn(k)
      dx = dabs(ttxs-xhn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xhn(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = yhn(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zhn(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz
      else
      udo = -ucr(rmm2)+udo               
      endif
      enddo
      
      kclp = 0
      do k = 1,ks-1
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz

      else
      udo = ucr(rmm2)+udo         
      endif
      enddo

      do k = ks+1,Nm
      xcp(k) = xhp(k)
      ycp(k) = yhp(k)
      zcp(k) = zhp(k)
      dx = dabs(ttxs-xhp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-yhp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zhp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xhp(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = yhp(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zhp(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz

      else
      udo = ucr(rmm2)+udo         
      endif
      enddo            

      udn = 0.d0
      kkcln = 0
      do k = 1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
      udn = -ucr(rmm2)+udn         
      endif   
      enddo

      kkclp = 0
      do k = 1,ks-1
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      udn = ucr(rmm2)+udn         
      endif   
      enddo
      do k = ks+1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
      udn = ucr(rmm2)+udn
      endif   
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcp(ks) = txs
      ycp(ks) = tys
      zcp(ks) = tzs      
      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xhn(km)
      tty = yhn(km)
      ttz = zhn(km)
      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo            
      enddo

      do k = 1,Nm      
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      udn = -ucr(rmm2)+udn
      udo = -uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo       
      enddo
      enddo

      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xhp(km)
      tty = yhp(km)
      ttz = zhp(km)      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) goto 927
      udn = ucr(rmm2)+udn
      udo = uc(ttx,tty,ttz,xhp(k),yhp(k),zhp(k))+udo            
      enddo

      do k = 1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) goto 927
      udn = -ucr(rmm2)+udn
      udo = -uc(ttx,tty,ttz,xhn(k),yhn(k),zhn(k))+udo       
      enddo
      enddo
      goto 929
 927  idonk = 1
 929  continue
      return
      end      
      

      subroutine cinit
      implicit double precision (a-h,o-z)      
      include 't.inc'                        
      do i = 1,Nm
      xhp(i) = (ran2(islu)-0.5d0)*ed2
      yhp(i) = (ran2(islu)-0.5d0)*ed2
      zhp(i) = (ran2(islu)-0.5d0)*ed2
      xhn(i) = (ran2(islu)-0.5d0)*ed2
      yhn(i) = (ran2(islu)-0.5d0)*ed2
      zhn(i) = (ran2(islu)-0.5d0)*ed2
      enddo
      write(*,*) 'Ions successfully inserted'                            
      return
      end

      subroutine dispp
      implicit double precision (a-h,o-z)
      include 't.inc' 
      idonk = 0
      if (ran2(islu).lt.0.1d0) then
      tz = ed2*(ran2(islu)-0.5d0)
      tx = ed2*(ran2(islu)-0.5d0)
      ty = ed2*(ran2(islu)-0.5d0)
      else         
      dtz = dpsp*(ran2(islu)-0.5d0)
      dtx = dpsp*(ran2(islu)-0.5d0)
      dty = dpsp*(ran2(islu)-0.5d0)
      tx = xhp(ks)+dtx
      if (tx.gt.ed) then
      tx = tx-ed2
      elseif (tx.lt.-ed) then
      tx = tx+ed2
      endif
      ty = yhp(ks)+dty
      if (ty.gt.ed) then
      ty = ty-ed2
      elseif (ty.lt.-ed) then
      ty = ty+ed2
      endif
      tz = zhp(ks)+dtz
      if (tz.gt.ed) then
      tz = tz-ed2
      elseif (tz.lt.-ed) then
      tz = tz+ed2
      endif
      endif
      
      do k = 1,Nm
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhspn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = 1,ks-1
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspp2)
      if (idonk.eq.1) goto 98
      enddo
      do k = ks+1,Nm
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspp2)
      if (idonk.eq.1) goto 98
      enddo

      un = 0.d0
      uo = 0.d0
      qttx = xhp(ks)
      qtty = yhp(ks)
      qttz = zhp(ks)
      do k = 1,Nm
      un = -uc(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = -uc(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo
      do k = 1,ks-1
      un = uc(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = uc(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
      do k = ks+1,Nm
      un = uc(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = uc(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
 98   continue
      return
      end

      subroutine dispn
      implicit double precision (a-h,o-z)
      include 't.inc' 
      idonk = 0
      if (ran2(islu).lt.0.1d0) then
      tz = ed2*(ran2(islu)-0.5d0)
      tx = ed2*(ran2(islu)-0.5d0)
      ty = ed2*(ran2(islu)-0.5d0)
      else               
      dtz = dpsn*(ran2(islu)-0.5d0)
      dtx = dpsn*(ran2(islu)-0.5d0)
      dty = dpsn*(ran2(islu)-0.5d0)
      tx = xhn(ks)+dtx
      if (tx.gt.ed) then
      tx = tx-ed2
      elseif (tx.lt.-ed) then
      tx = tx+ed2
      endif
      ty = yhn(ks)+dty
      if (ty.gt.ed) then
      ty = ty-ed2
      elseif (ty.lt.-ed) then
      ty = ty+ed2
      endif
      tz = zhn(ks)+dtz
      if (tz.gt.ed) then
      tz = tz-ed2
      elseif (tz.lt.-ed) then
      tz = tz+ed2
      endif
      endif

      do k = 1,Nm
      idonk = kdonk(tx,ty,tz,xhp(k),yhp(k),zhp(k),dhspn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = 1,ks-1
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhsnn2)
      if (idonk.eq.1) goto 98
      enddo
      do k = ks+1,Nm
      idonk = kdonk(tx,ty,tz,xhn(k),yhn(k),zhn(k),dhsnn2)
      if (idonk.eq.1) goto 98
      enddo

      un = 0.d0
      uo = 0.d0
      qttx = xhn(ks)
      qtty = yhn(ks)
      qttz = zhn(ks)
      do k = 1,Nm
      un = -uc(tx,ty,tz,xhp(k),yhp(k),zhp(k))+un
      uo = -uc(qttx,qtty,qttz,xhp(k),yhp(k),zhp(k))+uo      
      enddo
      do k = 1,ks-1
      un = uc(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = uc(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo
      do k = ks+1,Nm
      un = uc(tx,ty,tz,xhn(k),yhn(k),zhn(k))+un
      uo = uc(qttx,qtty,qttz,xhn(k),yhn(k),zhn(k))+uo      
      enddo
 98   continue
      return
      end

      subroutine ucalc
      implicit double precision (a-h,o-z)
      include 't.inc'         
      utot = 0.d0 
      do i = 1,Nm
      tx = xhp(i)
      ty = yhp(i)
      tz = zhp(i)
      if (abs(tx).gt.ed) write(*,*) '|XM|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YM|.GT.ED!',i,ty
      if (abs(tz).gt.ed) write(*,*) '|ZM|.GT.ED!',i,tz
      do k = i+1,Nm
      dx = abs(tx-xhp(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhp(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhp(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspp2) write(*,*) 'RMM2.LT.DHSPP2!',i,k,rmm2
      enddo
      do k = 1,Nm
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhspn2) write(*,*) 'RMM2.LT.DHSPN2!',i,k,rmm2
      enddo
      enddo

      do i = 1,Nm
      tx = xhp(i)
      ty = yhp(i)
      tz = zhp(i)
      if (abs(tx).gt.ed) write(*,*) '|XM|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YM|.GT.ED!',i,ty
      if (abs(tz).gt.ed) write(*,*) '|ZM|.GT.ED!',i,tz
      do k = i+1,Nm
      dx = abs(tx-xhp(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhp(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhp(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.edsq) then
      r = dsqrt(rmm2)
c      utot = 1.d0/r+r*redsq-shuss+utot
      sij = 1.d0-r*shuss+y5*r**5-y6*r**6+y7*r**7
      utot = sij/r+utot
      endif
      enddo
      do k = 1,Nm
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.edsq) then
      r = dsqrt(rmm2)
c      utot = -(1.d0/r+r*redsq-shuss)+utot
      sij = 1.d0-r*shuss+y5*r**5-y6*r**6+y7*r**7
      utot = -sij/r+utot
      endif
      enddo
      enddo

      do i = 1,Nm
      tx = xhn(i)
      ty = yhn(i)
      tz = zhn(i)
      if (abs(tx).gt.ed) write(*,*) '|XM|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YM|.GT.ED!',i,ty
      if (abs(tz).gt.ed) write(*,*) '|ZM|.GT.ED!',i,tz
      do k = i+1,Nm
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dhsnn2) write(*,*) 'RMM2.LT.DHSNN2!',im,k,rmm2
      enddo
      enddo
      do i = 1,Nm
      tx = xhn(i)
      ty = yhn(i)
      tz = zhn(i)
      if (abs(tx).gt.ed) write(*,*) '|XM|.GT.ED!',i,tx
      if (abs(ty).gt.ed) write(*,*) '|YM|.GT.ED!',i,ty
      if (abs(tz).gt.ed) write(*,*) '|ZM|.GT.ED!',i,tz
      do k = i+1,Nm
      dx = abs(tx-xhn(k))
      dx = dx-aint(dx/ed)*ed2
      dy = abs(ty-yhn(k))
      dy = dy-aint(dy/ed)*ed2
      dz = abs(tz-zhn(k))
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.edsq) then
      r = dsqrt(rmm2)
c      utot = 1.d0/r+r*redsq-shuss+utot
      sij = 1.d0-r*shuss+y5*r**5-y6*r**6+y7*r**7
      utot = sij/r+utot
      endif
      enddo
      enddo
      return
      end

      function kdonk(xx,yy,zz,xxx,yyy,zzz,dref2)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx/ed)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy/ed)*ed2
      dz = abs(zz-zzz)
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dref2) then
      kdonk = 1
      else
      kdonk = 0
      endif
      return
      end

      function ucr(r2)
      implicit double precision (a-h,o-z)
      include 't.inc'
      uu = 0.d0
      if (r2.lt.edsq) then      
      r = dsqrt(r2)
      r4 = r2**2
      uu = 1.d0/r-shuss+y5*r4-y6*r4*r+y7*r4*r2
      endif
      ucr = uu
      return
      end
      
      function uc(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx/ed)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy/ed)*ed2
      dz = abs(zz-zzz)
      dz = dz-aint(dz/ed)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      uu = 0.d0
      if (rmm2.lt.edsq) then
      r = dsqrt(rmm2)
c      uu = 1.d0/r+r*redsq-shuss
      r4 = rmm2**2
      uu = 1.d0/r-shuss+y5*r4-y6*r4*r+y7*r4*rmm2
      endif
      uc = uu
      return
      end


      function dist(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 't.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx/ed)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy/ed)*ed2
      dz = abs(zz-zzz)
      dz = dz-aint(dz/ed)*ed2
      dist = dsqrt(dx*dx+dy*dy+dz*dz)
      return
      end








