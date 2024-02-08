c This program written to compute the impulse response of a Maxwell
c Earth .... February, 1991
c modified (JXM) April 1991
c modified (JXM) May 1991
c modified (JXM) August 1991
c modified (JXM after AMF) April 1992
c modified (JXM ) October 1998
c .................................................................
c Program is based on Wu and Peltier (1992) ... most relevant
c equations can be found there. Inversion of Laplace transform is 
c performed using a full collocation technique. 
c Non-dimensionalization is performed as in Wu's thesis (MSc)
c Output file formats follow those used by WRP's codes so that
c I don't have to change all my peripheral programs !!!!
c last checked (7/92) for model 12012 with output from WRP code
c
c COLLOCATION VERSION - 06/05
c
c ....................................................................
       program maxwell
       implicit double precision (a-h, o-z)
       parameter (nla=850)
       character*80 char
       dimension rad(nla),rhod(nla),vel_p(nla),vel_s(nla),gacc(nla),
     .           visc(nla),zmu(nla),zlam(nla),
     .           zmu_s(nla),zlam_s(nla)
       dimension elast(3),viscous(3),prout(3)
       dimension elast_tide(3),viscous_tide(3),prout_tide(3)
       dimension spole(1000),amp(3,1000),amp_tide(3,1000)
       dimension slog_begin_array(500),slog_end_array(500),
     .            slog_delta_array(500),ztol_array(500)
       integer    mode_max_array(500),num_half_array(500)
c collocation stuff 
       dimension valh(200),vall(200),valk(200)
       dimension valh_tide(200),vall_tide(200),valk_tide(200)
       dimension zalh(200),zall(200),zalk(200)
       dimension zalh_tide(200),zall_tide(200),zalk_tide(200)
       dimension coll_mat(200,200),coll_mati(200,200)
       dimension coll_mat2(200,200)
       integer indx(200),modess
c common blocks ...
       common/earth/ rad,rhod,vel_p,vel_s,
     .                gacc,visc,zmu,zlam,ntot,ncore,naes
       common/searth/ zmu_s,zlam_s
       common/propog8/ prout,prout_tide,sdet,i
       common/bound/ bound1,bound2
       common/asympv/iasymp
c ....................................................................
c set up input/output files
       call askfil(3, ' input Earth model ')
       call askfil(7, ' output Love numbers etc .... ')
c ....................................................................
c set up some constants
       s_large=-1.0d12
       s_small=-1.0d-12
c .....................................................................
c read in Earth model
       ict=1
       ict2=1
 100   continue
       read(3,*,end=101) rad(ict),rhod(ict),vel_p(ict),vel_s(ict),
     .                   gacc(ict),visc(ict)
       if(rhod(ict).gt.7000.) ict2=ict2+1
       ict=ict+1
       go to 100
 101   continue
c compute the number of nodes in total and the number in the core 
       ntot=ict-1
       ncore=ict2-1
c      write(7,*) ntot,ncore
c ......................................................................
c first step is to non-dimensionalize the equations
c this follows Pat Wu's thesis !!!
       call dimen(rad,rhod,vel_p,vel_s,gacc,visc,ntot,ncore,naes)
c compute the elastic lame parameters
       call plam(vel_p,vel_s,rhod,zmu,zlam,ntot)
c .....................................................................
c read in header
       read(5,'(a)') char
       write(7,900) char
 900   format(a80)
c Next we read in the required degree range - this program permits
c different s-searches for each degree range and therefore as many of
c of these ranges as one specifies are possible ...
  990  continue
c read in required degree range
       read(5,*,end=991) l_begin,l_end,l_delta,nsval
      do i=1,nsval
        read(5,*) slog_begin_array(i),slog_end_array(i),
     .            slog_delta_array(i)
        enddo
c start the loop over spherical harmonic degree
      do i=l_begin,l_end,l_delta 
      write(6,*) i
c ...................................................................
c find elastic and viscous asymptotes
      iasymp=1
c for elastic asymptote set transform variable very high !!!
      call slame(s_large)
      call propo(s_large,3)
      do ijk=1,3
        elast(ijk)=prout(ijk)
        enddo
      do ijk=1,3
        elast_tide(ijk)=prout_tide(ijk)
        enddo
c for viscous asymptote set transform variable very small
      call slame(s_small)
      call propo(s_small,3)
      do ijk=1,3
        viscous(ijk)=prout(ijk)
        enddo
      do ijk=1,3
        viscous_tide(ijk)=prout_tide(ijk)
        enddo
c ....................................................................
      iasymp=0
      modes_found=0
c next loop over the Laplace transform domain 
       do kkk=1,nsval
c read in required transform variable range 
       slog_begin=slog_begin_array(kkk)
       slog_end=slog_end_array(kkk)
       slog_delta=slog_delta_array(kkk)
c maximum number of modes and number of times to half interval
c when a mode has been bounded
       numlog=aint(abs(slog_end-slog_begin)/slog_delta+0.0001)
c
        do j=1,numlog
          slval1=slog_begin+float(j-1)*slog_delta
c         sval1=-10.0**slval1
          sval1=10.0**slval1
          call slame(sval1)
          call propo(sval1,3)
          modes_found=modes_found+1
          spole(modes_found)=sval1
c loop over the different types of love numbers
              do ijk=1,3
                amp(ijk,modes_found)=prout(ijk)
                enddo
              do ijk=1,3
                amp_tide(ijk,modes_found)=prout_tide(ijk)
                enddo
          enddo
c next end kkk (different s-ranges) loop
        enddo

c now add the collocation scheme to compute the modal amplitudes
        do iv=1,modes_found
          do jv=1,modes_found
            coll_mat(iv,jv)=1./(spole(iv)+spole(jv))
c           write(6,*) iv,jv,coll_mat(iv,jv)
            enddo
          enddo

c invert the matrix
c coll_mat ... coll_mati

c       do iv=1,modes_found
c         do jv=1,modes_found
c           coll_mat2(iv,jv)=coll_mat(iv,jv)
c           enddo
c         enddo

        modess=modes_found
        call matinv(coll_mat,coll_mati,indx,200,modess)

c       write(6,*) i,i,viscous(1)
        do iv=1,modes_found
          valh(iv)=amp(1,iv)-elast(1)
c         write(6,*) spole(iv),hsval(iv)+elast(1)
          vall(iv)=amp(2,iv)-elast(2)
          valk(iv)=amp(3,iv)-elast(3)
          valh_tide(iv)=amp_tide(1,iv)-elast_tide(1)
          vall_tide(iv)=amp_tide(2,iv)-elast_tide(2)
          valk_tide(iv)=amp_tide(3,iv)-elast_tide(3)
          enddo
c       write(6,*) i,i,elast(1)

       do iv=1,modes_found
        amp(1,iv)=0.0
        amp(2,iv)=0.0
        amp(3,iv)=0.0
        amp_tide(1,iv)=0.0
        amp_tide(2,iv)=0.0
        amp_tide(3,iv)=0.0
        do jv=1,modes_found
         amp(1,iv)=amp(1,iv)+coll_mati(iv,jv)*valh(jv)
         amp(2,iv)=amp(2,iv)+coll_mati(iv,jv)*vall(jv)
         amp(3,iv)=amp(3,iv)+coll_mati(iv,jv)*valk(jv)
         amp_tide(1,iv)=amp_tide(1,iv)+coll_mati(iv,jv)*valh_tide(jv)
         amp_tide(2,iv)=amp_tide(2,iv)+coll_mati(iv,jv)*vall_tide(jv)
         amp_tide(3,iv)=amp_tide(3,iv)+coll_mati(iv,jv)*valk_tide(jv)
         coll_mat(iv,jv)=0.0
         enddo
        enddo

c check collocation
c      do iv=1,modes_found
c       zalh(iv)=0.0
c       zall(iv)=0.0
c       zalk(iv)=0.0
c       zalh_tide(iv)=0.0
c       zall_tide(iv)=0.0
c       zalk_tide(iv)=0.0
c       do jv=1,modes_found
c        zalh(iv)=zalh(iv)+coll_mat2(iv,jv)*amp(1,jv)
c        zall(iv)=zall(iv)+coll_mat2(iv,jv)*amp(2,jv)
c        zalk(iv)=zalk(iv)+coll_mat2(iv,jv)*amp(3,jv)
c        zalh_tide(iv)=zalh_tide(iv)+coll_mat2(iv,jv)*amp_tide(1,jv)
c        zall_tide(iv)=zall_tide(iv)+coll_mat2(iv,jv)*amp_tide(2,jv)
c        zalk_tide(iv)=zalk_tide(iv)+coll_mat2(iv,jv)*amp_tide(3,jv)
c        enddo
c       enddo
c
c      write(6,*) 'h Love number '
c      do iv=1,modes_found
c        write(6,*) spole(iv),zalh(iv),valh(iv)
c        enddo
c      write(6,*) 'l Love number '
c      do iv=1,modes_found
c        write(6,*) spole(iv),zall(iv),vall(iv)
c        enddo
c      write(6,*) 'k Love number '
c      do iv=1,modes_found
c        write(6,*) spole(iv),zalk(iv),valk(iv)
c        enddo

c
c write out everything for this degree
c this format is intended to mimic WRP's code so that I can
c use seaiter, 3-d, etc. I have not ye implemented tidal Love
c numbers but this is extremely easy !!!
      zdeg=float(i)
      zer=0.0
      write(7,1003) i, modes_found
      write(7,1020) (spole(ii),ii=1,modes_found) 
      write(7,1020) elast(1),elast(2),zer,zer,elast(3)
      write(7,1020) viscous(1),viscous(2),zer,zer,viscous(3)
      write(7,1020) (amp(1,iii),iii=1,modes_found)
      write(7,1020) (amp(2,iii),iii=1,modes_found)
      write(7,1020) (amp(3,iii),iii=1,modes_found)
      if(zdeg.gt.1.0) then
      write(7,1020) elast_tide(1),elast_tide(2),zer,zer,elast_tide(3)
      write(7,1020) viscous_tide(1),viscous_tide(2),zer,zer,
     .              viscous_tide(3)
      write(7,1020) (amp_tide(1,iii),iii=1,modes_found)
      write(7,1020) (amp_tide(2,iii),iii=1,modes_found)
      write(7,1020) (amp_tide(3,iii),iii=1,modes_found)
      endif
 1003 format(i10,1x,i5)
 1020 format(5e16.8)
c move to next degree
        enddo
c next degree range
      go to 990
  991 continue
      stop
      end

       subroutine dimen(rad,rhod,vel_p,vel_s,gacc,visc,ntot,ncore,naes)
c for details see Pat Wu's master's thesis (p.35) (also appendix)
       implicit double precision (a-h,o-z)
       parameter (nla=850)
       dimension rad(nla),rhod(nla),vel_p(nla),vel_s(nla),
     .    gacc(nla),visc(nla)
c set up some constants ...
c average Earth density
       rho_average=5520.0
       big_g=6.67d-11
       rad_earth=6371000.
       zfac=sqrt(3.14159*big_g*rho_average*rad_earth*rad_earth)
       zfac2=3.14159*big_g*rho_average*rad_earth
       zfac3=3.14159*big_g*rho_average*rho_average*rad_earth*
     .       rad_earth*3.153D10
       ict_core=0
       ict_aes=0
       do i=1,ntot
        if(visc(i).lt. 10.0) ict_core=ict_core+1
        if(visc(i).lt. 1.d30) ict_aes=ict_aes+1
        rad(i)=rad(i)/rad_earth
        rhod(i)=rhod(i)/rho_average
        vel_p(i)=vel_p(i)/zfac
        vel_s(i)=vel_s(i)/zfac
        gacc(i)=gacc(i)/zfac2
        visc(i)=visc(i)/zfac3
        enddo
c      write(7,*) ict_core,ict_aes
        ncore=ict_core
        naes=ict_aes
       return
       end
c subroutine plam simply computes the lame parameters from the 
c input p and s wave velocities
       subroutine plam(vel_p,vel_s,rhod,zmu,zlam,ntot)
       implicit double precision (a-h,o-z)
       parameter (nla=850)
       dimension vel_p(nla),vel_s(nla),rhod(nla)
       dimension zmu(nla), zlam(nla)
       do i=1,ntot
        zmu(i)=rhod(i)*vel_s(i)*vel_s(i)
        zlam(i)=rhod(i)*vel_p(i)*vel_p(i)-2.0*zmu(i)
        enddo
       return
       end
c subroutine slam computes the transform dependent (visco-elastic)
c lame parameters
       subroutine slame(svalue)
       implicit double precision (a-h, o-z)
       parameter (nla=850)
       dimension rad(nla),rhod(nla),vel_p(nla),vel_s(nla),gacc(nla),
     .           visc(nla),zmu(nla),zlam(nla),
     .           zmu_s(nla),zlam_s(nla)
       common /earth/rad,rhod,vel_p,vel_s,
     .                gacc,visc,zmu,zlam,ntot,ncore,naes
       common /searth/ zmu_s,zlam_s
        do i=1,ntot
         zk=zlam(i)+2./3.*zmu(i)
         ratio=zmu(i)/visc(i)
         zmu_s(i)=zmu(i)*svalue/(svalue+ratio)
         zlam_s(i)=(zlam(i)*svalue+zk*ratio)/(svalue+ratio)
         enddo
c      endif
        return
        end
c subroutine propo calls prop_soln which does the progogation etc..
        subroutine propo(sval,ii)
        implicit double precision (a-h, o-z)
       parameter (nla=850)
       dimension rad(nla),rhod(nla),vel_p(nla),vel_s(nla),gacc(nla),
     .           visc(nla),zmu(nla),zlam(nla),
     .           zmu_s(nla),zlam_s(nla)
       dimension elast(3),viscous(3),prout(3)
       dimension elast_tide(3),viscous_tide(3),prout_tide(3)
       dimension spole(1000),amp(3,1000),array(6),array2(6)
c common blocks
       common /earth/ rad,rhod,vel_p,vel_s,
     .                gacc,visc,zmu,zlam,ntot,ncore,naes
       common /searth/ zmu_s,zlam_s
       common /propog8/ prout,prout_tide,sdet,ideg
       common /bound/ bound1,bound2
c ...................................................................
c need to set up a check for the propogation - if the transform
c variable is very small then we assume the mantle below the 
c lithosphere is inviscid.
         zab=2.0*float(ideg)+1.0
         if(abs(sval).lt.1.d-10) then
           icheck=0
         else 
           icheck=1
         endif
c set up boundary conditions ...
         bound1=zab*gacc(ntot)
         bound2=-zab*gacc(ntot)*gacc(ntot)/4.0
c prop_soln is the program which does the propogation and boundary
c condition match. It outputs the determinant of hte boundary 
c condition matrix, as well as the final propogated vector
        if (ii.eq.3) then
          call prop_soln(icheck,ideg,sdet,array,array2)
          prout(1)=array(1)
          prout(2)=array(2)
          prout(3)=float(ideg)*(array(5)/gacc(ntot)-1.0)
          prout_tide(1)=array2(1)
          prout_tide(2)=array2(2)
          prout_tide(3)=float(ideg)*(array2(5)/gacc(ntot)-1.0)
          endif
        if (ii.eq.1) then
          call prop_soln(icheck,ideg,sdet,array,array2)
          endif
        if (ii.eq.2)  then
          call prop_soln(icheck,ideg,sdet,array,array2)
          prout(1)=array(1)*sdet
          prout(2)=array(2)*sdet
          prout(3)=float(ideg)*(array(5)/gacc(ntot))*sdet
          prout_tide(1)=array2(1)*sdet
          prout_tide(2)=array2(2)*sdet
          prout_tide(3)=float(ideg)*(array2(5)/gacc(ntot))*sdet
          endif
        return
        end
       subroutine prop_soln(icheck,ideg,det,soln,soln2)
       implicit double precision (a-h, o-z)
       parameter (nla=850)
c this subroutine gets starting solutions, propogates them to the surface,
c combines them to satisfy the boundary conditions, and computes the 
c determinant of the boundary condition matrix
      dimension crst(2,1),zmanst(6,3),bc(3),weight(3),bcm(3,3),
     .  soln(6),yy(6,3,nla),soln2(6),bcm2(3,3)
      common/earth/ rad(nla),rhod(nla),vel_p(nla),
     .             vel_s(nla),gacc(nla),vsc(nla),ztemp(nla),ztemp2(nla),
     .             ntot,ncore,naes
      common/searth/ zmu_s(nla),zlam_s(nla) 
      common/bound/ bound1,bound2
      common/asympv/iasymp
c
c ..................................................................
c compute the fundamental boundary for the propogation ...
      if(icheck.eq.0) then
       nbot=naes
       ntop=naes+1
      else
       nbot=ncore
       ntop=ncore+1
      endif
c start the propogations !!!
c a simple calculation tell us if we start in the core or the mantle
      xval=10.**(-6./float(ideg))
c In the following I force the propogation within the core.
c These are the following steps:
c 1. Calculate starting solution within the core (everything below is
c    assumed homogenious)
c 2. Propogate through the core 
c 3. Match solution to mantle
c 4. Propogate through mantle
c I have not written code to start in 
c the mantle. Therefore I force the start within the core. 
c The code is written with dummy programs in place so all you have to 
c do is write them and delete one line below (bounded by ***'s)
c ...................................................................
c zero out the eigenfunction vectors
c       write(6,*) ntot,ntot,ntot
        do  i=1,ntot
         do  j=1,3
          do  k=1,6
           yy(k,j,i)=0.0
           enddo
          enddo
         enddo
c ************
c     if(xval .gt. rad(ncore)) xval=rad(ncore)
c ************
      if(xval.le.rad(nbot)) then
c we start the propogation in the core
        xcore=xval
        if(xcore.lt. 0.15) xcore=0.15
c first thing to do is to find the starting node
        call rnode(xcore,inbeg,rbeg)
        rval=rad(inbeg)
c next get the starting vector in the core
        call core_st(inbeg,crst)
        yy(1,1,inbeg)=crst(1,1)/gacc(inbeg)
        yy(5,1,inbeg)=crst(1,1)
        yy(6,1,inbeg)=crst(2,1)
c       write(6,*) inbeg,ntop,nbot
c next propogate through the rest of the core
        inbeg=inbeg +1
        if(inbeg.gt.nbot) go to 20 
        do i=inbeg,nbot
         call cprop8(i,crst)
         yy(1,1,i)=crst(1,1)/gacc(i)
         yy(5,1,i)=crst(1,1)
         yy(6,1,i)=crst(2,1)
         enddo
 20     continue
c need to match vectors at cmb (or base of lithosphere)
       call matcm(nbot,crst,zmanst)
        do i=1,3
         do j=1,6
           yy(j,i,ntop)=zmanst(j,i)
           enddo
         enddo
c next propogate through mantle
       rval=rad(ntop)
       inbeg=ntop+1
       do i=inbeg,ntot
        call zmprop8(i,zmanst)
        do j=1,3        
         do k=1,6        
          yy(k,j,i)=zmanst(k,j)
          enddo
         enddo
        enddo
c       do i=1,ntot
c       write(6,*) (yy(k,j,i),k=1,6)
c922    format(6f12.6)
c        enddo
       else
c in this case you start the propogation in the mantle !!
c mantle_st NOT YET IMPLEMENTED
       xmantle=xval
       call rnode(xmantle,inbeg,rbeg)
       rval=rad(inbeg)
       call mantle_st(inbeg,zmanst)
       do i=1,3
        do j=1,6
          yy(j,i,inbeg)=zmanst(j,i)
          enddo
        enddo
      inbeg=inbeg+1
      do i=inbeg,ntot
        call zmprop8(i,zmanst)
        do j=1,3
          do k=1,6
           yy(k,j,i)=zmanst(k,j)
          enddo
         enddo
       enddo
      endif
C
c at this point we check to see if we have degree 1 or not
c
      if(ideg.eq.1) then
c
c degree = 1 case
c
c for l=1 we only need two vectors: either 1,3 or 2,3
c we will choose 2,3; furthermore, we only need two
c boundary conditions so we choose the first and second
c
c set up the boundary element matrix 
      do i=1,2
       bcm(1,i)=yy(3,i+1,ntot)
       bcm(2,i)=yy(4,i+1,ntot)
      enddo
c set up the boundary condition vector
      bc(1)=bound2
      bc(2)=0.
      do i=1,3
        weight(i)=0.0
        enddo
c the determinant and inverse can be done analytically for a 2X2 !
c     do i=1,2
c      write(7,*) bcm(i,1),bcm(i,2)
c     enddo
      det=bcm(1,1)*bcm(2,2)-bcm(1,2)*bcm(2,1)
      weight(1)=1./det*(bcm(2,2)*bc(1)-bcm(1,2)*bc(2))
      weight(2)=1./det*(-bcm(2,1)*bc(1)+bcm(1,1)*bc(2))
c combine the vectors for the final solution ...
      do i=1,6
        soln(i)=0.0
       do j=2,3
        soln(i)=soln(i)+yy(i,j,ntot)*weight(j-1)
        enddo
       enddo
c next compute the various shifts
      if(iasymp.eq.1)zshift=soln(5)/gacc(ntot)-1.0
      if(iasymp.eq.0)zshift=soln(5)/gacc(ntot)
c This is the correct version. The old version subtracted zshift:
c I subtracted because it agreed with Dahlen's numbers, but I think his
c numbers (for l=1) are wrong
      soln(1)=soln(1)-zshift
      soln(2)=soln(2)-zshift
      if(iasymp.eq.1) soln(5)=gacc(ntot)
      if(iasymp.eq.0) soln(5)=0.0
c     soln(5)=gacc(ntot)
c
      do i=1,6
        soln2(i)=0.0
       enddo
c
      else
c
c degree not equal to 1
c set up the boundary element matrix 
      do i=1,3
       bcm(1,i)=yy(3,i,ntot)
       bcm(2,i)=yy(4,i,ntot)
       bcm(3,i)=yy(6,i,ntot)
       bcm2(1,i)=yy(3,i,ntot)
       bcm2(2,i)=yy(4,i,ntot)
       bcm2(3,i)=yy(6,i,ntot)
      enddo
c set up the boundary condition vector
      bc(1)=bound2
      bc(2)=0.
      bc(3)=bound1
      do i=1,3
        weight(i)=0.0
        enddo
c Matrix inverts bcm and multiplies by bc to get the vector weights. It
c also computes the determinant of bcm
      call matrix(bcm,bc,weight,det)
c combine the vectors for the final solution ...
      do i=1,6
        soln(i)=0.0
       do j=1,3
        soln(i)=soln(i)+yy(i,j,ntot)*weight(j)
        enddo
       enddo
c set up the boundary condition vector
      bc(1)=0.0      
      bc(2)=0.
      bc(3)=bound1
      do i=1,3
        weight(i)=0.0
        enddo
c Matrix inverts bcm and multiplies by bc to get the vector weights. It
c also computes the determinant of bcm
      call matrix(bcm2,bc,weight,det2)
c combine the vectors for the final solution ...
      do i=1,6
        soln2(i)=0.0
       do j=1,3
        soln2(i)=soln2(i)+yy(i,j,ntot)*weight(j)
        enddo
       enddo
      endif
c
      return  
      end
c rnode derives the node counter for a radius just below the input rval
       subroutine rnode(rval,inode,rnval)
       implicit double precision (a-h, o-z)
       parameter (nla=850)
       common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .                vel_s(nla),gacc(nla),visc(nla),zmu(nla),zlam(nla),
     .                ntot,ncore,naes
       inode=0
       do i=1,ntot
        if(rad(i).le.rval) then
         inode=inode+1
         endif
        if(rad(i).eq.rval) go to 222
        enddo
 222    continue
        rnval=rad(inode)
        return
        end
        subroutine matcm(nval,cmat,dmat)
        implicit double precision (a-h, o-z)
       parameter (nla=850)
        common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .                vel_s(nla),gacc(nla),visc(nla),zmu(nla),zlam(nla),
     .                ntot,ncore,naes
        dimension cmat(2,1),dmat(6,3)
c this routine matches propogated solutions at the core to
c those in the mantle at the cmb: see Pat Wu's thesis
        do i=1,3
         do j=1,6
          dmat(j,i)=0.0
          enddo
          enddo 
c vector 1
         dmat(1,1)=cmat(1,1)/gacc(nval)
         dmat(5,1)=cmat(1,1)
         dmat(6,1)=cmat(2,1)
c vector 2
         dmat(2,2)=1.0
c vector 3
         dmat(1,3)=1.
         dmat(3,3)=rhod(nval)*gacc(nval)
         dmat(6,3)=-4.*rhod(nval)
       return
       end
c subroutine to compute the vector weights and the determinant of
c the boundary element matrix
      subroutine matrix (a,bc,c,d)
      implicit double precision (a-h, o-z)
      dimension a(3,3), bc(3), c(3) 
      dimension b(3,3)
      real*8 d

      integer  indx(3), i, j
     

      do i = 1,3
	do j =1, 3
          b(i,j) = 0.d0
        end do
	b(i,i) = 1.d0
      end do

      call ludcmp (a, 3, 3, indx, d)
      do j = 1,3
      call lubksb(a,3, 3, indx, b(1,j))
      end do

        do10 i=1,3 
         sum=0.0
         do20 j=1,3 
          sum=sum+b(i,j)*bc(j)
 20       continue
         c(i)=sum
 10      continue
        
      do j= 1,3
	d = d*a(j,j)
      end do
      return
      end

c these are numerical recipes programs .....

      subroutine ludcmp(a,n,np,indx,d)
      integer nmax
      real*8 tiny

      parameter(nmax = 200, tiny = 1.d-30)

      real*8 a(np,np),  vv(nmax), aamax, d, sum, dum

      integer indx(np), n, np, i, j, k, imax

      d = 1.d0
      do i = 1, n
	aamax = 0.d0
	do j = 1, n
	 if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
        end do
c        if (aamax .eq. 0.d0) pause 'singular matrix'
	if (aamax .eq. 0.d0) write(6,*) 'singular matrix'
	vv(i) = 1.d0/aamax
      end do

      do j = 1, n
	do i = 1, j-1
	  sum = a(i,j)
	  do k = 1, i-1
	    sum = sum - a(i,k)*a(k,j)
          end do
	  a(i,j) = sum
        end do

	aamax = 0.d0
	do i = j, n
	  sum = a(i,j)
	  do k = 1, j-1
	    sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
	  dum = vv(i)*abs(sum)
	  if (dum .ge. aamax) then
	     imax = i
	     aamax = dum
          end if
        end do

	if (j .ne. imax) then
	  do k = 1, n
	    dum = a(imax, k)
	    a(imax,k) = a(j,k)
	    a(j,k) = dum
          end do
	  d = -d
	  vv(imax) = vv(j)
        end if

	indx(j) = imax
	if (a(j,j) .eq. 0.d0)a(j,j) =tiny
	   if (j .ne. n)then
	     dum = 1.d0/a(j,j)
	     do i = j+1, n
	       a(i,j) = a(i,j)*dum
             end do
           end if
      end do
      return
      end

      subroutine lubksb(a, n, np, indx, b)
      
      real*8 a(np, np), b(np), sum
      integer indx(np), i, ll, j, ii, n, np

      ii = 0
      do i = 1, n
	ll = indx(i)
	sum = b(ll)
	b(ll) = b(i)
	if (ii .ne. 0) then
	   do j = ii, i-1
	     sum = sum - a(i,j)*b(j)
           end do
        else if (sum .ne. 0.d0) then
	    ii = i
        end if
	b(i) = sum
      end do

      do i = n, 1, -1
	 sum = b(i)
	 do j = i+1, n
	    sum = sum - a(i,j)*b(j)
         end do
	 b(i) = sum/a(i,i)
      end do
      return
      end

c this is a program to propogate a vector through the core

      subroutine cprop8(itop,cvec)
      implicit double precision (a-h,o-z)
      parameter (nla=850)
      dimension rkw(6),rka(3)
      dimension amat(2,2),cvec(2,1),temp(2,1),rkmat(3,2,1)
      dimension temp2(2,1)
      common /earth/rad(nla),rhod(nla),vel_p(nla),
     .         vel_s(nla),gacc(nla),vsc(nla),ztemp(nla),ztemp2(nla),
     .         ntot,itemp,itemp2
      common /propog8/ prout(3),prout_tide(3),sdet,ideg
      data rka/0.0,0.5,1.0/
      data rkw/0.5,-1.0,2.0,0.166667,0.666667,0.166667/
c ............................................
      ibot=itop-1
      rad_bot=rad(ibot)
      rad_top=rad(itop)
      del_rad=rad_top-rad_bot
      if(rad_bot .ge. rad_top) return
c zero out the 3-dim vector
      do i=1,3
        do j=1,2
          do k=1,1
            rkmat(i,j,k)=0.0
            enddo
          enddo
        enddo
c get a copy of cvec
       do j=1,1
         do k=1,2
           temp2(k,j)=cvec(k,j)
           enddo
         enddo
c start the prpogation (simple 3-ord)
      do i=1,3
c get the radial value
        radv=rad_bot + rka(i)*del_rad
c get the propogator matrix at this value
        call prmat_core(radv,itop,amat)
c multiply the propogator matrix times the soln vectors
        do j=1,1
          do k=1,2
            temp(k,j)=0.0
            do jk=1,2
              temp(k,j)=temp(k,j)+amat(k,jk)*cvec(jk,j)
              enddo
            enddo
          enddo
c copy the results
       do j=1,1
         do k=1,2
           rkmat(i,k,j)=temp(k,j)
           enddo
         enddo
c reset cvec
       do j=1,1
         do k=1,2
           cvec(k,j)=temp2(k,j)
           enddo
         enddo
c 
       do j=1,1
         do k=1,2
           do jk=1,i
             ict=i*(i-1)/2 + jk
             cvec(k,j)=cvec(k,j)+rkw(ict)*del_rad*rkmat(jk,k,j)
             enddo
            enddo
          enddo
        enddo
      return
      end

      subroutine core_st(inbeg,crst)
      implicit double precision  (a-h,o-z)
      parameter (nla=850)
c this program gives the starting solution for a propogation
c begun in the core. It is very simple - see Wu & Peltier, p448
c
       dimension crst(2,1)
       common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .                vel_s(nla),gacc(nla),visc(nla),zmu(nla),zlam(nla),
     .                ntot,ncore,naes
       common/propog8/ prout(3),prout_tide(3),sdet,ideg
       zi=float(ideg) 
       val=rad(inbeg)**ideg
       val2=rad(inbeg)**(ideg-1)
       do i=1,2
        do j=1,1
         crst(i,j)=0.0
         enddo
        enddo
      crst(1,1)=val
      crst(2,1)=2.0*(float(ideg)-1.0)*val2
      return
      end

c this program propogates themantle vectors through that region

      subroutine zmprop8(itop,zmvec)
      implicit double precision (a-h,o-z)
       parameter (nla=850)
      dimension rkw(6), rka(3)
      dimension amat(6,6),zmvec(6,3),temp(6,3),rkmat(3,6,3)
      dimension temp2(6,3)
      common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .         vel_s(nla),gacc(nla),vsc(nla),ztemp(nla),ztemp2(nla),
     .         ntot,itemp,itemp2
      common /searth/ zmu_s(nla),zlam_s(nla) 
      common /propog8/  prout(3),prout_tide(3),sdet,ideg
c     common/countx/iclo
      data rka/0.0,0.5,1.0/
      data rkw/0.5,-1.0,2.0,0.166667,0.666667,0.166667/
c ............................................
      iclo=iclo+1
      ibot=itop-1
      rad_bot=rad(ibot)
      rad_top=rad(itop)
      del_rad=rad_top-rad_bot
      if(rad_bot .ge. rad_top) return
c zero out the 3-dim vector
      do i=1,3
        do j=1,6
          do k=1,3
            rkmat(i,j,k)=0.0
            enddo
          enddo
        enddo
c get a copy of zmvec
       do j=1,3
         do k=1,6
           temp2(k,j)=zmvec(k,j)
           enddo
         enddo
c start the prpogation (simple 3-ord)
      do i=1,3
c get the radial value
        radv=rad_bot + rka(i)*del_rad
c get the propogator matrix at this value
        call prmat(radv,itop,amat)
c multiply the propogator matrix times the soln vectors
        do j=1,3
          do k=1,6
            temp(k,j)=0.0
            do jk=1,6
              temp(k,j)=temp(k,j)+amat(k,jk)*zmvec(jk,j)
              enddo
            enddo
          enddo
c copy the results
       do j=1,3
         do k=1,6
           rkmat(i,k,j)=temp(k,j)
           enddo
         enddo
c reset zmvec
       do j=1,3
         do k=1,6
           zmvec(k,j)=temp2(k,j)
           enddo
         enddo
c 
       do j=1,3
         do k=1,6
           do jk=1,i
             ict=i*(i-1)/2 + jk
             zmvec(k,j)=zmvec(k,j)+rkw(ict)*del_rad*rkmat(jk,k,j)
             enddo
            enddo
          enddo
        enddo
c      do i=1,3
c        do j=1,6
c          write(6,*) i,j,zmvec(j,i)
c          enddo
c        enddo
c      write(6,*) zzq,iclo,zmvec(2,2)
       return
       end

c this is simply the propogator matrix: see Wu MSc appendix
      subroutine prmat(rval,itop,amat)
      implicit double precision (a-h, o-z)
      parameter (nla=850)
      dimension amat(6,6)
      common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .         vel_s(nla),gacc(nla),vsc(nla),ztemp(nla),ztemp2(nla),
     .         ntot,itemp,itemp2
      common /searth/ zmu_s(nla),zlam_s(nla) 
      common /propog8/ prout(3),prout_tide(3),sdet,ideg
c ............................................
      ibot=itop-1
      zdeg=float(ideg)
c     if(ideg.eq.1)zdeg=zdeg+0.001
      del=(rval-rad(ibot))/(rad(itop)-rad(ibot))
      den=rhod(ibot)+del*(rhod(itop)-rhod(ibot))
      ymu=zmu_s(ibot)+del*(zmu_s(itop)-zmu_s(ibot))
      ylam=zlam_s(ibot)+del*(zlam_s(itop)-zlam_s(ibot))
      grav=gacc(ibot)+del*(gacc(itop)-gacc(ibot))
      comb=(ylam+2.0*ymu)
      gam=ymu*(3.0*ylam+2.0*ymu)/(ylam+2.0*ymu)
c end of preliminaries ...
      amat(1,1)=-2.0*ylam/(comb*rval)
      amat(1,2)=ylam*(zdeg+1.0)/(comb*rval)
      amat(1,3)=1.0/comb
      amat(1,4)=0.0
      amat(1,5)=0.0
      amat(1,6)=0.0
      amat(2,1)=-zdeg/rval
      amat(2,2)=1.0/rval
      amat(2,3)=0.0
      amat(2,4)=1.0/ymu
      amat(2,5)=0.0
      amat(2,6)=0.0
      amat(3,1)=4.0/rval*(gam/rval-den*grav)
      amat(3,2)=(den*grav-2.0*gam/rval)*(zdeg+1.0)/rval
      amat(3,3)=-4.0*ymu/(comb*rval)
      amat(3,4)=(zdeg+1.0)/rval
      amat(3,5)=den*(zdeg+1.0)/rval
      amat(3,6)=-den
      amat(4,1)=zdeg/rval*(den*grav-2.0*gam/rval)
      amat(4,2)=-1.0/(rval*rval)*(2.0*ymu-(gam+ymu)*zdeg*(zdeg+1.0))
      amat(4,3)=-zdeg*ylam/(comb*rval)
      amat(4,4)=-3.0/rval
      amat(4,5)=-zdeg*den/rval
      amat(4,6)=0.0
      amat(5,1)=4.0*den
      amat(5,2)=0.0
      amat(5,3)=0.0
      amat(5,4)=0.0
      amat(5,5)=-(zdeg+1.0)/rval
      amat(5,6)=1.0
      amat(6,1)=4.0*den*(zdeg+1.0)/rval
      amat(6,2)=-4.0*den*(zdeg+1.0)/rval
      amat(6,3)=0.0
      amat(6,4)=0.0
      amat(6,5)=0.0
      amat(6,6)=(zdeg-1.0)/rval
      return
      end
c this is simply the propogator matrix: see Wu MSc appendix
      subroutine prmat_core(rval,itop,amat)
      implicit double precision (a-h, o-z)
      parameter (nla=850)
      dimension amat(2,2)
      common /earth/ rad(nla),rhod(nla),vel_p(nla),
     .         vel_s(nla),gacc(nla),vsc(nla),ztemp(nla),ztemp2(nla),
     .         ntot,itemp,itemp2
      common /propog8/ prout(3),prout_tide(3),sdet,ideg
c ............................................
      ibot=itop-1
      zdeg=float(ideg)
      del=(rval-rad(ibot))/(rad(itop)-rad(ibot))
      den=rhod(ibot)+del*(rhod(itop)-rhod(ibot))
      grav=gacc(ibot)+del*(gacc(itop)-gacc(ibot))
c end of preliminaries ...
      amat(1,1)=4.*den/grav-(zdeg+1.)/rval
      amat(1,2)=1.0
      amat(2,1)=8.*den*(zdeg-1.)/(grav*rval)
      amat(2,2)=(zdeg-1.)/rval-4.*den/grav
c ...........................................
      return
      end
C	*******************************
	SUBROUTINE ASKFIL(CHANN,STRING)
C	*******************************
C
	INTEGER*4 CH,CHANN,ERROR
	INTEGER*4 IERR
	CHARACTER*80 NAME
	CHARACTER *(*) STRING
C
C  Write out prompting string.
C
	WRITE(*,'(A,$)') STRING
C
C  Read name of file which is to be opened.
C
	READ(*,'(A)') NAME
C
C  Open file.
C
	IF(IABS(CHANN).EQ.5 .OR. IABS(CHANN).EQ.6) THEN
	    PRINT *,'USING CHANNELS 5 AND 6 IS NOT ADVISED'
	END IF
	IF(CHANN.GT.0) THEN
	  OPEN(CHANN,FILE=NAME,BLANK='ZERO')
	ELSE IF(CHANN.EQ.0) THEN
	  PRINT *, 'CHANNEL 0 IS NOT ALLOWED'
	  STOP
	ELSE
          OPEN(-CHANN,FILE=NAME,STATUS='OLD',BLANK='ZERO')
	  PRINT *, 'READ ONLY CHANNEL OPENED'
	END IF
C
	RETURN
	END

        subroutine matinv(a,y,indx,np,n)
        implicit double precision (a-h,o-z)
        dimension a(np,np),y(np,np)
        integer indx(np)
        do12 i=1,n
         do11 j=1,n
          y(i,j)=0.0
 11       continue
         y(i,i)=1.0
 12      continue
        call ludcmp(a,n,np,indx,d)
        do13 j=1,n
         call lubksb(a,n,np,indx,y(1,j))
 13      continue
        return
        end
c
c MANTLE_ST
c this is a program to compute mantle starting solutions
c this version assumes an incompressible homogeneous starting
c model
      subroutine mantle_st(inbeg,zmanst)
      implicit double precision (a-h,o-z)
      parameter(nla=850)
      dimension zmanst(6,3)

      common/searth/ zmu_s(nla),zlam_s(nla)
      common/propog8/ prout(3),prout_tide(3),sdet,ideg
      common/earth/ rad(nla),rhod(nla),vel_p(nla),vel_s(nla),
     .      gacc(nla),visc(nla),zmu(nla),zlam(nla),ntot,ncore,naes

       do i=1,6
        do j=1,3
         zmanst(i,j)=0.0
         enddo
        enddo
c     write(6,*) 4.*rhod(inbeg),gacc(inbeg)/rad(inbeg)
      zp=rad(inbeg)**(ideg+1)
      z0=rad(inbeg)**ideg
      z1=rad(inbeg)**(ideg-1)
      z2=rad(inbeg)**(ideg-2)
      zsqr=sqrt(float(ideg)*(float(ideg)+1.))/float(ideg)
c vector 1
      zmanst(1,1)=0.0
      zmanst(2,1)=0.0
      zmanst(3,1)=rhod(inbeg)*z0
      zmanst(4,1)=0.0
      zmanst(5,1)=z0
      zmanst(6,1)=(2.*ideg+1)*z1
c vector 2
      zmanst(1,2)=z1 
      zmanst(2,2)=z1*zsqr
      zmanst(3,2)=rhod(inbeg)*gacc(inbeg)*z1+
     .            2.*(ideg-1)*zmu_s(inbeg)*z2
      zmanst(4,2)=2.*(ideg-1)*zmu_s(inbeg)*z2/float(ideg)
      zmanst(5,2)=0.0
      zmanst(6,2)=3.*gacc(inbeg)*z2
c vector 3
      zmanst(1,3)=ideg*zp/(2.*(2.*ideg+3))
      zmanst(2,3)=ideg*(ideg+3)/(2.*(2.*ideg+3)*(ideg+1))*zp*zsqr
      zmanst(3,3)=(rhod(inbeg)*gacc(inbeg)*ideg*zp+
     .            2.*(ideg**2-ideg-3)*zmu_s(inbeg)*z0)/
     .            (2.*(2.*ideg+3))
      zmanst(4,3)=2.*(ideg)*(ideg+2)*zmu_s(inbeg)*z0/
     .            (2.*(2.*ideg+3)*(ideg+1))
      zmanst(5,3)=0.0
      zmanst(6,3)=3.*gacc(inbeg)*z0*ideg/
     .            (2.*(2.*ideg+3))
c     do i=1,6
c      write(6,*) (zmanst(i,j),j=1,3)
c      enddo
      return
      end

