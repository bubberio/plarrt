!  Console1.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

       program debris
       implicit double precision (a-h ,o-z)
       dimension x(6), v(6), y(6)
       dimension xx(10,6), vv(10,6), ff(10,6), gg(10,6)
       dimension xs(8), xm(4), orb(7), adata(25)
       dimension conv(1)
       double precision tsid, pi2, degree, ageo
       double precision mj2000

       pi2=8.d0*datan(1.d0)
       rad= pi2/360.d0
       degree= pi2/360.d0
       ageo=42164.1696d0

        open(1,file='fli_a_M.plt',
     !  status='unknown',  form='formatted')
     
         open(2, file='fli_a_M_v2.plt',
     !  status='unknown',  form='formatted')


        open(40,file='sparam.plt',
     !  status='unknown',  form='formatted')


*************************************************************************
*      INPUT DATA from file input.txt
*************************************************************************

        open (6, file='inputs.txt', form='formatted')


        read(6,*)
        read(6,*)
        read(6,*)
        read(6,*)
        read(6,*)
        read(6,*) totdays
        read(6,*) mj2000
        read(6,*) andays
         read(6,*)
        read(6,*) ! a
        read(6,*) e
        read(6,*) ain
        read(6,*) ! amean
        read(6,*) om
        read(6,*) bom
        read(6,*) are2
        read(6,*) gam
        read(6,*)
        read(6,*) ipr
        read(6,*) isrp
        read(6,*) ieh
        read(6,*) isun
        read(6,*) imoon

       close(6, STATUS='KEEP')

        write(40,*) totdays
        write(40,*) mj2000
        write(40,*) andays
        write(40,*) !a
        write(40,*) e
        write(40,*) ain
        write(40,*) om
        write(40,*) bom
        write(40,*) ! amean
        write(40,*) are2
        write(40,*) gam
        write(40,*) ipr
        write(40,*) isrp
        write(40,*) ieh
        write(40,*) isun
        write(40,*) imoon


***********************************************************************
*      TIME
***********************************************************************
*      N=the integer part of the days (the total time of the integration)
        N=int(totdays)


*     The step of integration. "andays" represents the number of minutes:
*     We transform andays introduced above in our unit of time
       h=pi2*dble(andays)/(24.d0*60.d0)


*      Julian date


       dateju=mj2000+2400000.5d0+51545.d0

*       call djulian(iyear, month, nday, jhour, minu, seco, dateju)
*      Rescale the time in our unit of time. Choose the origin  of time (time=0) at dateju= 2451545 (or Epoch J2000)
       time0=(dateju-2451545.d0)*pi2*366.242196d0/365.242196d0
*********************************************************************
*      Parameters of the problem

      call param(adata, are2, ipr, isrp,ieh, isun, imoon)
***********************************************************************
*************************************************************************
*      Initial position of debris

*     bom=longitude of ascending node, (Big OMega=bom)
*     om=argument of perihelion,  (OMega=om)
*     ain=inclination,
*     amean=mean anomaly
*     a=semimajor axis in km
*     e=eccentricity

*      The kilometers are transformed in our units
*      Degrees are transformed in radians and kilometers in "g" units

        bom=bom*degree
        om= om*degree
        ain=ain*degree
        !  amean=amean*degree     
        !  a=a/ageo

*************************************************************************
       Npoint=20
       do iaa=1,Npoint+1
       a=(42120.d0+(dble(iaa-1)*100.d0)/dble(Npoint))/ageo
!       a=(32120.d0+(dble(iaa-1)*20000.d0)/dble(Npoint))/ageo

       do iee=1,Npoint+1
!	write(*,*)iaa,iee
       amean=(-20.d0+(dble(iee-1)*190.d0)/dble(Npoint))*rad
              
        
*       Solve the Kepler's equation
*       eec is the eccentric anomaly and fan is the true anomaly
        eec=0.d0
        do iii=1,100
        eec=amean+e*dsin(eec)
        end do
        fan=2.d0*datan(dsqrt((1.d0+e)/(1.d0-e))*dtan(eec/2.d0))


        write(40, *) iaa, iee
**************************************************************************************

*Initial position and velocity of the Satellite (or Debris) (in cartesian coordinates)

       call indeb(x, a, e, bom, om, ain, fan)

*        write(*,*) x(1), x(2), x(3)

*********************************************************************

* Initial vector for the variational problem (here, we are not interested in computing FLI,
* so please do not care about this)

        do ii=1,6
       y(ii)=x(ii)
       end do


        call rukugf(adata, y, v,  time0, h)



        v(1)=((y(2)-x(2))*cos(ain)+(y(3)-x(3))*cos(bom)*sin(ain))
       v(2)=((y(3)-x(3))*sin(bom)*sin(ain)-(y(1)-x(1))*cos(ain))
       v(3)=(-(y(1)-x(1))*cos(bom)*sin(ain)
     c       -(y(2)-x(2))*sin(bom)*sin(ain))
        vr=dsqrt(v(1)**2+v(2)**2+v(3)**2)
        v(1)=v(1)/vr
       v(2)=v(2)/vr
       v(3)=v(3)/vr
       v(4)=0.d0
       v(5)=0.d0
       v(6)=0.d0


       pnorm=0.d0
       do ii=1,6
         pnorm=pnorm+v(ii)**2
       enddo
       pnorm = dsqrt(pnorm)
       do ii=1,6
         v(ii)=v(ii)/pnorm
       enddo

       suppnorm = -10.d0

!        write(*,*) v(1), v(2), v(3),v(4),v(5),v(6)
***************************************************************************************


*           call orbital(x, orb)
*           call sun(xs, time0)
*           call moon(xm, time0)
           

*        write(1,5) time0, orb(1), dmod(orb(6), 360.d0)
*  5   format (' ', 3(F20.10), ' ')
*        write(2,10) time0, orb(2), dmod(orb(4), 360.d0)
*  10   format (' ', 3(F20.10), ' ')
*        write(3,15) time0, orb(3), dmod(orb(5), 360.d0)
*  15   format (' ', 3(F20.10), ' ')


*********************************************************************
*First data are computed by using the Kunge-Kutta method.

       call ainit(adata, x, v, xx, vv, ff, gg, h, time0)

***********************************************************************************
*  A-B-M (Adams-Bashforth-Multon) method for the initial value problem
*  and the variational problem.


       eps=0.1d-8
       iter=20
       count=0.d0


       tsid=time0
       
         do i=1,300000000
         


         call sun(xs, tsid)
         call moon(xm, tsid)



       call adams(adata, xx, vv, ff, gg, xm, xs, tsid, h,
     -  iter, eps, conv)

        tdays=(tsid-time0)/pi2      ! sidereal days
        tdays2=tdays*365.24219d0/366.24219d0            ! solar days



       if (conv(1) .lt. 1.d0) then
       write(40, *) " Adams-Moulton not converge at time =", tdays
       end if

       



       !   if   (tdays2 .ge. count) then
       !    count=count+1.d0

        ! do iii=1,6
        ! x(iii)=xx(8,iii)
        ! end do
           
         !  call orbital(x, orb)

       


        !res21=dmod(orb(6)-2.d0*tsid/degree+2*orb(5)+2.d0*orb(4), 360.d0)
        !if (res21 .lt. 0.d0) then
        !res21=res21+360.d0
        !end if

        ! tdays3=tdays2 +mj2000         ! sidereal days

        !write(1,30) tdays3, orb(1)*ageo, dmod(orb(6), 360.d0)
        
 !  30   format (' ', 3(F20.10), ' ')
           
 !          write(2,35) tdays3, orb(2), dmod(orb(4), 360.d0)
 !  35   format (' ', 3(F20.10), ' ')
!*           write(3,*) tyears, dsqrt(amu*(1-beta)*orb(1)/semaxJ), anbe
!           write(3,40) tdays3, orb(3), dmod(orb(5), 360.d0)
!   40   format (' ', 3(F20.10), ' ')

!           write(40,55) tdays3, conv(1)
!   55    format (' ', 2(F10.5), ' ')

!          end if


         tsid=tsid+h


         if (tdays .ge. totdays) then
         go to 342
         end if


         
         
            pnorm=0.d0
	  do ii=1,6
	     pnorm=pnorm+vv(8,ii)**2
	  enddo
	  pnorm=0.5d0*dlog10(pnorm)
	  suppnorm = max(dabs(pnorm),suppnorm)


          if (i .lt. 500) then
          suppnorm=0.d0
          end if

          if (suppnorm .gt. 7.d0) then
          goto 342
          end if


        
         
         if (tdays .ge. totdays) then
         go to 342
         end if


         end do  ! over i=1,N
         
         

 342     continue
 
 
        write(1,*) amean/rad, a*ageo, suppnorm
        
         write(2,*) amean/rad, a*ageo, suppnorm
        
        end do
        write(1,*)

        end do

        end







*************************************************************************
       subroutine ainit(adata, x, v, xx, vv, ff, gg, h, time0)
       implicit double precision (a-h, o-z)
       dimension v(6), x(6)
       dimension xx(10,6), vv(10, 6), ff(10, 6), gg(10,6)
       dimension ga11(6), xs(8), xm(4), adata(25)
       dimension xtemp(6), vtemp(6), ftemp(6), gtemp(6)


        hr=-h/10000.d0
        tsid=time0
        do ii=1,6
         xtemp(ii)=x(ii)
         vtemp(ii)=v(ii)
        end do
        
        
        do it=1,9
        ait=dble(it-1)
        
         call sun(xs, tsid)
         call moon(xm, tsid)

         call  fun(adata, xtemp, xm, xs,  tsid, ftemp)
         call fung(adata, xtemp, vtemp, xm, xs,  tsid, gtemp)

         do ii=1,6
         xx(10-it, ii)=xtemp(ii)
         vv(10-it, ii)=vtemp(ii)
         ff(10-it, ii)=ftemp(ii)
         gg(10-it, ii)=gtemp(ii)
         end do

        do i=1,10000
        call rukugf(adata, xtemp, vtemp, tsid, hr)
        tsid=time0-ait*h+hr*dble(i)
        end do

        end do ! end it

        return
        end



****************************************************************************
*                        The Adams method
****************************************************************************

       subroutine adams(adata, xx, vv, ff, gg, xm, xs, tsid, h,
     -  iter, eps, conv)
       implicit double precision (a-h, o-z)
       dimension xx(10,6), vv(10, 6), ff(10, 6), gg(10,6)
       dimension xs(8), xm(4), adata(25)
       dimension conv(1)

       call adba(adata, xx, vv, ff, gg, xm, xs, tsid, h)
       call admo(adata, xx, vv, ff, gg, xm, xs, tsid, h,
     -  iter, eps, conv)
     

       do jj=1,9
       do ii=1,6
       xx(jj,ii)=xx(jj+1,ii)
       vv(jj,ii)=vv(jj+1,ii)
       ff(jj,ii)=ff(jj+1,ii)
       gg(jj,ii)=gg(jj+1,ii)
       end do
       end do
       

       return
       end





****************************************************************************
*                        The Adams-Bashforth method
****************************************************************************
       subroutine adba(adata, xx, vv, ff, gg, xm, xs, tsid, h)
       implicit double precision (a-h, o-z)
       dimension xx(10,6), vv(10, 6), ff(10, 6), gg(10,6)
       dimension xs(8), xm(4), adata(25)
       dimension ff10(6), gg10(6), xx10(6), vv10(6)
       




        a0=14097247.d0
        a1=-43125206.d0
        a2= 95476786.d0
        a3=-139855262.d0
        a4= 137968480.d0
        a5= -91172642.d0
        a6= 38833486.d0
        a7= -9664106.d0
        a8=1070017.d0





       do ii=1,6
        xx(10,ii)=xx(9,ii)+(h/3628800.d0)*(a0*ff(9,ii)
     -   +a1*ff(8,ii)+a2*ff(7,ii) +a3*ff(6,ii)+a4*ff(5,ii)
     -   +a5*ff(4,ii)+a6*ff(3,ii) +a7*ff(2,ii)+a8*ff(1,ii))
        end do


       do ii=1,6
       vv(10,ii)=vv(9,ii)+(h/3628800.d0)*(a0*gg(9,ii)
     -   +a1*gg(8,ii)+a2*gg(7,ii) +a3*gg(6,ii)+a4*gg(5,ii)
     -   +a5*gg(4,ii)+a6*gg(3,ii) +a7*gg(2,ii)+a8*gg(1,ii))
        end do


        do ii=1,6
        xx10(ii)=xx(10, ii)
        end do
        
        do ii=1,6
        vv10(ii)=vv(10, ii)
        end do

       call fun(adata, xx10, xm, xs, tsid, ff10)
       call fung(adata, xx10, vv10, xm, xs, tsid, gg10)
       
       do ii=1,6
        ff(10, ii)=ff10(ii)
        end do

        do ii=1,6
        gg(10, ii)=gg10(ii)
        end do

       return
       end
       

****************************************************************************
*                        The Adams-Moulton method
****************************************************************************

       subroutine admo(adata, xx, vv, ff, gg, xm, xs, tsid, h,
     -  iter, eps, conv)
       implicit double precision (a-h, o-z)
       dimension xx(10,6), vv(10, 6), ff(10, 6), gg(10,6)
       dimension xs(8), xm(4), adata(25)
       dimension ff10(6), gg10(6), xx10(6), vv10(6)
       dimension xold(6), vold(6), xtemp(6), vtemp(6)
       dimension conv(1)
       


        
        b0=1070017.d0
        b1= 4467094.d0
        b2= -4604594.d0
        b3= 5595358.d0
        b4= -5033120.d0
        b5=3146338.d0
        b6= -1291214.d0
        b7= 312874.d0
        b8= -33953.d0
        
        


       do ii=1,6
        xtemp(ii)=xx(9,ii)+(h/3628800.d0)*(b1*ff(9,ii)
     -   +b2*ff(8,ii) +b3*ff(7,ii)+b4*ff(6,ii)
     -   +b5*ff(5,ii)+b6*ff(4,ii) +b7*ff(3,ii)+b8*ff(2,ii))
        end do


       do ii=1,6
        vtemp(ii)=vv(9,ii)+(h/3628800.d0)*(b1*gg(9,ii)
     -   +b2*gg(8,ii) +b3*gg(7,ii)+b4*gg(6,ii)
     -   +b5*gg(5,ii)+b6*gg(4,ii) +b7*gg(3,ii)+b8*gg(2,ii))
        end do


        do it=1,iter

        do ii=1,6
         xold(ii)=xx(10,ii)
         vold(ii)=vv(10,ii)
         xx(10,ii)=xtemp(ii)+(h/3628800.d0)*b0*ff(10,ii)
         vv(10,ii)=vtemp(ii)+(h/3628800.d0)*b0*gg(10,ii)
        end do


        do ii=1,6
        xx10(ii)=xx(10, ii)
        end do

        do ii=1,6
        vv10(ii)=vv(10, ii)
        end do

       call fun(adata, xx10, xm, xs, tsid, ff10)
       call fung(adata, xx10, vv10, xm, xs, tsid, gg10)

       do ii=1,6
        ff(10, ii)=ff10(ii)
        end do

        do ii=1,6
        gg(10, ii)=gg10(ii)
        end do

        if (hconv(xx10,xold, eps) .gt. 0.1d0)  then
*         conv(1)=dabs(xx10(1)-xold(1))+dabs(xx10(2)-xold(2))
          conv(1)=dble(it)
         go to 255
        else
          conv(1)=0.d0
        end if
        

       end do  ! end iter

 255   continue

       return
       end
*****************************************************************************
*                 hconv

*  This function returns 1.d0 (true) if the relative difference between the current
*  and previous estimate of xx21 is less than eps when | xx21 | and | xold | are greater than 1 or
*     if the absolute difference of xx21 and xold is less than eps when
*     | xx21 | and | xold | is less than 1.  If neither condition is true the
*    function returns 0.d0 (false).
*****************************************************************************
       function hconv(xx10,xold, eps)
       implicit double precision (a-h ,o-z)
       dimension xx10(6), xold(6)
       
        xxnorm=dsqrt(xx10(1)**2+xx10(2)**2+xx10(3)**2)
        xonorm=dsqrt(xold(1)**2+xold(2)**2+xold(3)**2)
        xn=dsqrt((xold(1)-xx10(1))**2+(xold(2)-xx10(2))**2+
     -  (xold(3)-xx10(3))**2)

*       if((xxnorm .gt. 1.d0) .and. (xonorm .gt. 1.d0)) then
*       eps=eps*xxnorm
*       end if

       if ((dabs(xold(1)-xx10(1)) .lt. eps) .and.
     -   (dabs(xold(2)-xx10(2)) .lt. eps) .and.
     -    (dabs(xold(3)-xx10(3)) .lt. eps) .and.
     - (dabs(xold(4)-xx10(4)) .lt. eps) .and.
     -   (dabs(xold(5)-xx10(5)) .lt. eps) .and.
     -    (dabs(xold(6)-xx10(6)) .lt. eps)) then
       hconv=1.d0
       else
       hconv=0.d0
       end if

       return
       end function
       
**************************************************************************
*********************************************************************
       subroutine indeb(x, a, e, bom, om, ain, fan)
       implicit double precision(a-h, o-z)
       dimension x(6)

C       pi2=8.d0*datan(1.d0)
C       degree= pi2/360.d0
C       eps=23.4392911d0*degree

       r=a*(1.d0-e*e)/(1.d0+e*cos(fan))
       rdot=(1.d0/dsqrt(a*(1.d0-e*e)))*e*sin(fan)
       rfdot=(1.d0/dsqrt(a*(1.d0-e*e)))*(1.d0+e*cos(fan))
       x(1)=r*(cos(bom)*cos(om+fan)-sin(bom)*sin(om+fan)*cos(ain))

         yy=r*(sin(bom)*cos(om+fan)+cos(bom)*sin(om+fan)*cos(ain))
         zz=r*sin(om+fan)*sin(ain)
         x(2)=yy
         x(3)=zz
*       x(2)=yy*cos(eps)-zz*sin(eps)
*       x(3)=yy*sin(eps)+zz*cos(eps)

       x(4)=rdot*(cos(bom)*cos(om+fan)-sin(bom)*sin(om
     !   +fan)*cos(ain))-rfdot*(cos(bom)*sin(om+fan)
     !   +sin(bom)*cos(om+fan)*cos(ain))
       yydot=rdot*(sin(bom)*cos(om+fan)+cos(bom)*sin(om
     !   +fan)*cos(ain))-rfdot*(sin(bom)*sin(om+fan)
     !   -cos(bom)*cos(om+fan)*cos(ain))
       zzdot=rdot*sin(om+fan)*sin(ain)
     !   +rfdot*cos(om+fan)*sin(ain)
       x(5)=yydot
       x(6)=zzdot


       return
       end
*********************************************************************
******************************************************************************

        subroutine rukugf(adata, x, v, tsid, h)
        implicit double precision(a-h,o-z)
        dimension x(6), v(6)
        dimension z(6), adata(25)
        dimension bfive(6), bsix(6), bseven(6)
        dimension bone(6), btwo(6), bthree(6), bfour(6), btemp(6)
        dimension b2temp(6), b3temp(6), b4temp(6)
        dimension b5temp(6), b6temp(6), b7temp(6)
        dimension xs(8), xsh3(8), xs2h3(8)
        dimension xsh2(8), xsh(8)
        dimension xm(4), xmh3(4), xm2h3(4)
        dimension xmh2(4), xmh(4)


         do ii=1,6
         z(ii)=x(ii)
         end do

         call sun(xs, tsid)

         call sun(xsh3, tsid+h/3.d0)

         call sun(xs2h3, tsid+2.d0*h/3.d0)

         call sun(xsh2, tsid+h/2.d0)

         call sun(xsh, tsid+h)


         call moon(xm, tsid)

         call moon(xmh3, tsid+h/3.d0)

         call moon(xm2h3, tsid+2.d0*h/3.d0)

         call moon(xmh2, tsid+h/2.d0)

         call moon(xmh, tsid+h)



         call fun(adata, z, xm, xs, tsid, bone)


          do ii=1,6
          b2temp(ii)=z(ii)+(h/3.d0)*bone(ii)
          end do
          call fun(adata, b2temp , xmh3, xsh3, tsid+h/3.d0,  btwo)


          do ii=1,6
          b3temp(ii)=z(ii)+(2.d0*h/3.d0)*btwo(ii)
          end do
        call fun(adata, b3temp, xm2h3, xs2h3, tsid+2.d0*h/3.d0,  bthree)


          do ii=1,6
          b4temp(ii)=z(ii)+h*(bone(ii)+4.d0*btwo(ii)-bthree(ii))/12.d0
          end do
          call fun(adata, b4temp, xmh3, xsh3,  tsid+h/3.d0,  bfour)



          do ii=1,6
          b5temp(ii)=z(ii)+h*(-bone(ii)+18.d0*btwo(ii)-
     -     3.d0*bthree(ii)-6.d0*bfour(ii))/16.d0
          end do
          call fun(adata, b5temp, xmh2, xsh2, tsid+h/2.d0,  bfive)




          do ii=1,6
          b6temp(ii)=z(ii)+h*(9.d0*btwo(ii)-
     -     3.d0*bthree(ii)-6.d0*bfour(ii)+4.d0*bfive(ii))/8.d0
          end do
          call fun(adata, b6temp, xmh2, xsh2, tsid+h/2.d0,  bsix)


          do ii=1,6
          b7temp(ii)=z(ii)+h*(9.d0*bone(ii)-36.d0*btwo(ii)+
     -     63.d0*bthree(ii)+72.d0*bfour(ii)-64.d0*bfive(ii))/44.d0
          end do
          call fun(adata, b7temp,  xmh, xsh, tsid+h,  bseven)



       do ii=1,6
       x(ii)=z(ii)+(h*(11.d0*bone(ii)+81.d0*bthree(ii)+81.d0*bfour(ii)-
     -    32.d0*bfive(ii)-32.d0*bsix(ii)+11.d0*bseven(ii)))/120.d0
       end do






          call fung(adata, z, v, xm, xs, tsid, bone)


          do ii=1,6
          btemp(ii)=v(ii)+(h/3.d0)*bone(ii)
          end do
        call fung(adata, b2temp, btemp , xmh3, xsh3, tsid+h/3.d0,  btwo)


          do ii=1,6
          btemp(ii)=v(ii)+(2.d0*h/3.d0)*btwo(ii)
          end do
       call fung(adata, b3temp, btemp, xm2h3, xs2h3,
     -   tsid+2.d0*h/3.d0, bthree)


          do ii=1,6
          btemp(ii)=v(ii)+h*(bone(ii)+4.d0*btwo(ii)-bthree(ii))/12.d0
          end do
          call fung(adata, b4temp, btemp, xmh3, xsh3,
     -      tsid+h/3.d0,  bfour)



          do ii=1,6
          btemp(ii)=v(ii)+h*(-bone(ii)+18.d0*btwo(ii)-
     -     3.d0*bthree(ii)-6.d0*bfour(ii))/16.d0
          end do
       call fung(adata, b5temp, btemp, xmh2, xsh2, tsid+h/2.d0,  bfive)




          do ii=1,6
          btemp(ii)=v(ii)+h*(9.d0*btwo(ii)-
     -     3.d0*bthree(ii)-6.d0*bfour(ii)+4.d0*bfive(ii))/8.d0
          end do
          call fung(adata, b6temp, btemp, xmh2, xsh2, tsid+h/2.d0, bsix)


          do ii=1,6
          btemp(ii)=v(ii)+h*(9.d0*bone(ii)-36.d0*btwo(ii)+
     -     63.d0*bthree(ii)+72.d0*bfour(ii)-64.d0*bfive(ii))/44.d0
          end do
          call fung(adata, b7temp, btemp, xmh, xsh, tsid+h,  bseven)




       do ii=1,6
       v(ii)=v(ii)+(h*(11.d0*bone(ii)+81.d0*bthree(ii)+81.d0*bfour(ii)-
     -    32.d0*bfive(ii)-32.d0*bsix(ii)+11.d0*bseven(ii)))/120.d0
       end do


       return
       end


********************************************************************
        subroutine fun(adata, x, xm, xs,  tsid, f)
	implicit double precision(a-h,o-z)
	dimension adata(25), x(6),  f(6), xm(4), xs(8)

        pi2=8.d0*datan(1.d0)
        degree= pi2/360.d0



         amiumo=adata(1)
         amius=adata(2)

         a2erth=adata(3)
         aearth=adata(4)

         C20=adata(5)
         C22=adata(6)
         S22=adata(7)
         C30=adata(8)
         C31=adata(9)
         S31=adata(10)
         C32=adata(11)
         S32=adata(12)
         C33=adata(13)
         S33=adata(14)
         
         C40=adata(20)

         are2ms=adata(15)
         crpra2=adata(16)
         cont=adata(17)
         gamma=adata(18)
         slight=adata(19)


        tsidb=280.4606d0*pi2/360.d0+tsid
*       tsidb=tsid

	t2sid=2.d0*tsidb

	t3sid=3.d0*tsidb

	cosnew=C22*cos(t2sid)-S22*sin(t2sid)
        sinnew=C22*sin(t2sid)+S22*cos(t2sid)

        cosw=C31*cos(tsidb)-S31*sin(tsidb)
        sinw=C31*sin(tsidb)+S31*cos(tsidb)

        coew=C32*cos(t2sid)-S32*sin(t2sid)
        siew=C32*sin(t2sid)+S32*cos(t2sid)


        cos3w=C33*cos(t3sid)-S33*sin(t3sid)
        sin3w=C33*sin(t3sid)+S33*cos(t3sid)

         x1=x(1)
         x2=x(2)
         x3=x(3)



        r=dsqrt(x(1)**2+x(2)**2+x(3)**2)
        rdm=dsqrt((x(1)-xm(1))**2+(x(2)-xm(2))**2+(x(3)-xm(3))**2)
        rds=dsqrt((x(1)-xs(1))**2+(x(2)-xs(2))**2+(x(3)-xs(3))**2)
        rm=xm(4)
        rs=xs(4)
        rsdot=xs(8)

        xdovc=x(4)/slight
        ydovc=x(5)/slight
        zdovc=x(6)/slight
        rhdovc=(dsqrt(x(4)**2+x(5)**2+x(6)**2))/slight


        f(1)=x(4)
        f(2)=x(5)
        f(3)=x(6)


         f(4)=  -(x1/r**3) + (a2erth*
     -     (6.d0*cosnew*x1 + 6.d0*sinnew*x2 +
     - (15.d0*x1*(-2.d0*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/r**2 +
     -   C20*((3.d0*x1)/2.d0 - (15.d0*x1*x3**2)/(2.d0*r**2))))/r**5
     c      -amiumo*((x(1)-xm(1))/(exp(3.d0*dlog(rdm)))
     c         +xm(1)/(exp(3.d0*dlog(rm))))
     c      -amius*((x(1)-xs(1))/(exp(3.d0*dlog(rds)))
     c         +xs(1)/(exp(3.d0*dlog(rs))))
     c     +cont*crpra2*are2ms*((x(1)-xs(1))/(exp(3.d0*dlog(rds))))
     c     -crpra2*are2ms*gamma*((x(1)-xs(1))/(exp(3.d0*dlog(rds))))*
     c      (rhdovc-rsdot)
     c     -crpra2*are2ms*gamma*((xdovc-xs(5))/(exp(2.d0*dlog(rds))))


        f(5)=-(x2/r**3) + (a2erth*
     -  (6.d0*sinnew*x1 - 6.d0*cosnew*x2 +
     -  (15.d0*x2*(-2.d0*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/r**2 +
     -       C20*((3.d0*x2)/2.d0 - (15.d0*x2*x3**2)/(2.d0*r**2))))/r**5
     c      -amiumo*((x(2)-xm(2))/(exp(3.d0*dlog(rdm)))
     c         +xm(2)/(exp(3.d0*dlog(rm))))
     c      -amius*((x(2)-xs(2))/(exp(3.d0*dlog(rds)))
     c         +xs(2)/(exp(3.d0*dlog(rs))))
     c     +cont*crpra2*are2ms*((x(2)-xs(2))/(exp(3.d0*dlog(rds))))
     c     -crpra2*are2ms*gamma*((x(2)-xs(2))/(exp(3.d0*dlog(rds))))*
     c      (rhdovc-rsdot)
     c     -crpra2*are2ms*gamma*((ydovc-xs(6))/(exp(2.d0*dlog(rds))))

        f(6)=    -(x3/r**3) + (a2erth*
     - ((15.d0*(-2.d0*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2))*x3)/r**2 +
     -       C20*((9.d0*x3)/2.d0 - (15.d0*x3**3)/(2.d0*r**2))))/r**5
     c      -amiumo*((x(3)-xm(3))/(exp(3.d0*dlog(rdm)))
     c         +xm(3)/(exp(3.d0*dlog(rm))))
     c      -amius*((x(3)-xs(3))/(exp(3.d0*dlog(rds)))
     c         +xs(3)/(exp(3.d0*dlog(rs))))
     c     +cont*crpra2*are2ms*((x(3)-xs(3))/(exp(3.d0*dlog(rds))))
     c     -crpra2*are2ms*gamma*((x(3)-xs(3))/(exp(3.d0*dlog(rds))))*
     c      (rhdovc-rsdot)
     c     -crpra2*are2ms*gamma*((zdovc-xs(7))/(exp(2.d0*dlog(rds))))

       f(4)=f(4)+
     -  (-3.d0*aearth**3*(-5.d0*sinw*x1*x2*(r**2 - 7.d0*x3**2) +
     -      cosw*(r**4 + 35.d0*x1**2*x3**2 -
     -         5.d0*r**2*(x1**2 + x3**2))))/(2.d0*r**9)


       f(5)=f(5)+(3.d0*aearth**3*(5.d0*cosw*x1*x2*(r**2 - 7.d0*x3**2) -
     -      sinw*(r**4 + 35.d0*x2**2*x3**2 -
     -         5.d0*r**2*(x2**2 + x3**2))))/(2.d0*r**9)

       f(6)=f(6)+   (15.d0*aearth**3*(cosw*x1 + sinw*x2)*x3*
     -    (3.d0*r**2 - 7.d0*x3**2))/(2.d0*r**9)


         f(4)=f(4)+(15.d0*aearth**3*(2.d0*siew*(r**2 - 7.d0*x1**2)*x2 +
     -      coew*x1*(2.d0*r**2 - 7.d0*x1**2 + 7.d0*x2**2))*x3)/r**9

        f(5)=f(5)+  (15.d0*aearth**3*(2.d0*siew*x1*(r**2 - 7.d0*x2**2) +
     -      coew*x2*(-2.d0*r**2 - 7.d0*x1**2 + 7.d0*x2**2))*x3)/r**9

          f(6)=f(6)+ (15.d0*aearth**3*(2.d0*siew*x1*x2 +
     -    coew*(x1**2 - x2**2))*(r**2 - 7.d0*x3**2))/r**9


        f(4)=f(4)+ (15*aearth**3*
     -    (cos3w*(-7*x1**2*(x1**2 - 3*x2**2) +
     -         3*r**2*(x1**2 - x2**2)) +
     -      sin3w*x1*x2*
     -       (6*r**2 + 7*(-3*x1**2 + x2**2))))/r**9


       f(5)=f(5)+   (15*aearth**3*
     -    (cos3w*x1*x2*
     -       (-6*r**2 - 7*x1**2 + 21*x2**2) +
     -      sin3w*(3*r**2*(x1**2 - x2**2) +
     -         7*x2**2*(-3*x1**2 + x2**2))))/r**9


       f(6)=f(6)+  (-105*aearth**3*
     -    (cos3w*x1**3 + 3*sin3w*x1**2*x2 -
     -      3*cos3w*x1*x2**2 - sin3w*x2**3)*x3)/r**9

        f(4)=f(4)+ (5*aearth**3*C30*x1*x3*(3*r**2 - 7*x3**2))/(2.*r**9)

        f(5)=f(5)+ (5*aearth**3*C30*x2*x3*(3*r**2 - 7*x3**2))/(2.*r**9)

        f(6)=f(6)-(aearth**3*C30*(3*r**4 -
     -   30*r**2*x3**2 + 35*x3**4))/(2.*r**9)


       f(4)=f(4)+(-15*aearth**4*C40*x1*(r**4 -
     -  14*r**2*x3**2 + 21*x3**4))/
     -  (8.*r**11)

        f(5)=f(5)+(-15*aearth**4*C40*x2*(r**4 -
     -   14*r**2*x3**2 + 21*x3**4))/
     -  (8.*r**11)

        f(6)=f(6)+ (-5*aearth**4*C40*(15*r**4*x3 -
     -   70*r**2*x3**3 + 63*x3**5))/
     -  (8.*r**11)
        
        return
        end
*********************************************************************

*********************************************************************

       subroutine fung(adata, x, v, xm, xs, tsid, g)
       implicit double precision(a-h, o-z)
       dimension adata(25), x(6), v(6),  g(6), xs(8), xm(4)


         pi2=8.d0*datan(1.d0)
        degree= pi2/360.d0


          
         amiumo=adata(1)
         amius=adata(2)

         a2erth=adata(3)
         aearth=adata(4)

         C20=adata(5)
         C22=adata(6)
         S22=adata(7)
         C30=adata(8)
         C31=adata(9)
         S31=adata(10)
         C32=adata(11)
         S32=adata(12)
         C33=adata(13)
         S33=adata(14)
         
         C40=adata(20)

         are2ms=adata(15)
         crpra2=adata(16)
         cont=adata(17)
         gamma=adata(18)
         slight=adata(19)



        tsidb=280.4606d0*pi2/360.d0+tsid


	t2sid=2.d0*tsidb

	t3sid=3.d0*tsidb

	cosnew=C22*cos(t2sid)-S22*sin(t2sid)
        sinnew=C22*sin(t2sid)+S22*cos(t2sid)

        cosw=C31*cos(tsidb)-S31*sin(tsidb)
        sinw=C31*sin(tsidb)+S31*cos(tsidb)

        coew=C32*cos(t2sid)-S32*sin(t2sid)
        siew=C32*sin(t2sid)+S32*cos(t2sid)


        cos3w=C33*cos(t3sid)-S33*sin(t3sid)
        sin3w=C33*sin(t3sid)+S33*cos(t3sid)




         x1=x(1)
         x2=x(2)
         x3=x(3)



        rr=dsqrt(x(1)**2+x(2)**2+x(3)**2)
        rdm=dsqrt((x(1)-xm(1))**2+(x(2)-xm(2))**2+(x(3)-xm(3))**2)
        rds=dsqrt((x(1)-xs(1))**2+(x(2)-xs(2))**2+(x(3)-xs(3))**2)
        rm=xm(4)
        rs=xs(4)

        rsdot=xs(8)

        xdovc=x(4)/slight
        ydovc=x(5)/slight
        zdovc=x(6)/slight
        rhdovc=(dsqrt(x(4)**2+x(5)**2+x(6)**2))/slight



                 df4x= -rr**(-3) + (3*x1**2)/rr**5 -
     -  (5*a2erth*x1*(6*cosnew*x1 + 6*sinnew*x2 +
     -   (15*x1*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -    rr**2 + C20*((3*x1)/2. - (15*x1*x3**2)/(2.*rr**2))))/rr**7
     -    + (a2erth*(6*cosnew +
     -    (15*x1*(-2*cosnew*x1 - 2*sinnew*x2))/rr**2 +
     -    (15*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/rr**2 -
     -    (30*x1**2*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -     rr**4 + C20*(1.5 - (15*x3**2)/(2.*rr**2) +
     -     (15*x1**2*x3**2)/rr**4)))/rr**5
     c  +(amiumo*(3.d0*(x(1)-xm(1))**2-rdm**2))/rdm**5
     c  +((amius-crpra2*are2ms)*(3.d0*(x(1)-xs(1))**2
     c   -rds**2))/rds**5



       df4y= (3*x1*x2)/rr**5 + (a2erth*
     -   (6*sinnew + (15*x1*(-2*sinnew*x1 + 2*cosnew*x2))/rr**2 -
     -   (30*x1*x2*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**4 + (15*C20*x1*x2*x3**2)/rr**4))/rr**5 -
     -  (5*a2erth*x2*(6*cosnew*x1 + 6*sinnew*x2 +
     -   (15*x1*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**2 + C20*((3*x1)/2. - (15*x1*x3**2)/(2.*rr**2))))/rr**7
     c  + (amiumo*3.d0*(x(1)-xm(1))*(x(2)-xm(2)))/rdm**5
     c  + ((amius-crpra2*are2ms)*3.d0*(x(1)-xs(1))*(x(2)
     c  -xs(2)))/rds**5



       df4z=(3*x1*x3)/rr**5 - (5*a2erth*x3*
     -  (6*cosnew*x1 + 6*sinnew*x2 +
     -  (15*x1*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**2 + C20*((3*x1)/2. - (15*x1*x3**2)/(2.*rr**2))))/rr**7
     -  + (a2erth*((-30*x1*
     -   (-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2))*x3)/rr**4 +
     -  C20*((-15*x1*x3)/rr**2 + (15*x1*x3**3)/rr**4)))/rr**5
     c  + (amiumo*3.d0*(x(1)-xm(1))*(x(3)-xm(3)))/rdm**5
     c  + ((amius-crpra2*are2ms)*3.d0*(x(1)-xs(1))*(x(3)
     c  -xs(3)))/rds**5




       df5y=  -rr**(-3) + (3*x2**2)/rr**5 -
     -  (5*a2erth*x2*(6*sinnew*x1 - 6*cosnew*x2 +
     -   (15*x2*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**2 + C20*((3*x2)/2. - (15*x2*x3**2)/(2.*rr**2))))/rr**7
     -   + (a2erth*(-6*cosnew +
     -   (15*x2*(-2*sinnew*x1 + 2*cosnew*x2))/rr**2 +
     -   (15*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/rr**2 -
     -   (30*x2**2*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**4 + C20*(1.5 - (15*x3**2)/(2.*rr**2) +
     -   (15*x2**2*x3**2)/rr**4)))/rr**5
     c  +(amiumo*(3.d0*(x(2)-xm(2))**2-rdm**2))/rdm**5
     c  +((amius-crpra2*are2ms)*(3.d0*(x(2)
     c  -xs(2))**2-rds**2))/rds**5



       df5z=(3*x2*x3)/rr**5 - (5*a2erth*x3*
     -   (6*sinnew*x1 - 6*cosnew*x2 +
     -   (15*x2*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/
     -   rr**2 + C20*((3*x2)/2. - (15*x2*x3**2)/(2.*rr**2))))/rr**7
     -    + (a2erth*((-30*x2*
     -   (-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2))*x3)/rr**4 +
     -   C20*((-15*x2*x3)/rr**2 + (15*x2*x3**3)/rr**4)))/rr**5
     c  + (amiumo*3.d0*(x(2)-xm(2))*(x(3)-xm(3)))/rdm**5
     c  + ((amius-crpra2*are2ms)*3.d0*(x(2)
     c  -xs(2))*(x(3)-xs(3)))/rds**5

       df6z= -rr**(-3) + (3*x3**2)/rr**5 -
     -  (5*a2erth*x3*((15*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2))*
     -   x3)/rr**2 + C20*((9*x3)/2. - (15*x3**3)/(2.*rr**2))))/
     -   rr**7 + (a2erth*((15*
     -   (-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2)))/rr**2 -
     -   (30*(-2*sinnew*x1*x2 + cosnew*(-x1**2 + x2**2))*x3**2)/
     -   rr**4 + C20*(4.5 - (45*x3**2)/(2.*rr**2) +
     -   (15*x3**4)/rr**4)))/rr**5
     c  +(amiumo*(3.d0*(x(3)-xm(3))**2-rdm**2))/rdm**5
     c  +((amius-crpra2*are2ms)*(3.d0*(x(3)
     c  -xs(3))**2-rds**2))/rds**5


       df4x=df4x+ (15*aearth**3*C30*x1**2*x3)/(rr**2)**4.5 +
     -  (5*aearth**3*C30*x3*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**4.5) -
     -  (45*aearth**3*C30*x1**2*x3*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**5.5)

       df4y=df4y+ (15*aearth**3*C30*x1*x2*x3)/(rr**2)**4.5 -
     -  (45*aearth**3*C30*x1*x2*x3*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**5.5)

       df4z=df4z+ (-20*aearth**3*C30*x1*x3**2)/(rr**2)**4.5 +
     -  (5*aearth**3*C30*x1*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**4.5) -
     -  (45*aearth**3*C30*x1*x3**2*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**5.5)

       df5y=df5y+(15*aearth**3*C30*x2**2*x3)/(rr**2)**4.5 +
     -  (5*aearth**3*C30*x3*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**4.5) -
     -  (45*aearth**3*C30*x2**2*x3*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**5.5)

       df5z=df5z+(-20*aearth**3*C30*x2*x3**2)/(rr**2)**4.5 +
     -  (5*aearth**3*C30*x2*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**4.5) -
     -  (45*aearth**3*C30*x2*x3**2*(3*rr**2 - 7*x3**2))/
     -   (2.*(rr**2)**5.5)

       df6z=df6z -(aearth**3*C30*(-48*rr**2*x3 + 80*x3**3))/
     -   (2.*(rr**2)**4.5) +
     -  (9*aearth**3*C30*x3*(3*rr**4 - 30*rr**2*x3**2 + 35*x3**4))/
     -   (2.*(rr**2)**5.5)


        df4x=df4x+ (27*aearth**3*x1*
     -     (-5*sinw*x1*x2*(x1**2 + x2**2 - 6*x3**2) +
     -       cosw*(rr**4 + 35*x1**2*x3**2 -
     -          5*rr**2*(x1**2 + x3**2))))/(2.*rr**11) -
     -  (3*aearth**3*(-10*sinw*x1**2*x2 -
     -       5*sinw*x2*(x1**2 + x2**2 - 6*x3**2) +
     -       cosw*(-6*rr**2*x1 + 70*x1*x3**2 -
     -          10*x1*(x1**2 + x3**2))))/(2.*rr**9)



          df4x=df4x+(15*aearth**3*(-10*coew*x1**2 - 24*siew*x1*x2 +
     -       coew*(2*rr**2 - 7*x1**2 + 7*x2**2))*x3)/rr**9 -
     -  (135*aearth**3*x1*x3*(coew*x1*(2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x2*(-6*x1**2 + x2**2 + x3**2)))/rr**11


         df4x=df4x+ (15*aearth**3*(-30*sin3w*x1**2*x2 +
     -       cos3w*(6*rr**2*x1 - 14*x1**3 -
     -          14*x1*(x1**2 - 3*x2**2) + 6*x1*(x1**2 - x2**2)) +
     -       sin3w*x2*(6*rr**2 + 7*(-3*x1**2 + x2**2))))/rr**9 -
     -  (135*aearth**3*x1*(cos3w*
     -        (-7*x1**2*(x1**2 - 3*x2**2) +
     -          3*rr**2*(x1**2 - x2**2)) +
     -       sin3w*x1*x2*(6*rr**2 + 7*(-3*x1**2 + x2**2))))/rr**11



        df4y=df4y+  (27*aearth**3*x2*
     -     (-5*sinw*x1*x2*(x1**2 + x2**2 - 6*x3**2) +
     -       cosw*(rr**4 + 35*x1**2*x3**2 -
     -          5*rr**2*(x1**2 + x3**2))))/(2.*rr**11) -
     -  (3*aearth**3*(-10*sinw*x1*x2**2 -
     -       5*sinw*x1*(x1**2 + x2**2 - 6*x3**2) +
     -       cosw*(4*rr**2*x2 - 10*x2*(x1**2 + x3**2))))/
     -   (2.*rr**9)


        df4y=df4y+  (15*aearth**3*x3*(18*coew*x1*x2 + 4*siew*x2**2 +
     -       2*siew*(-6*x1**2 + x2**2 + x3**2)))/rr**9 -
     -  (135*aearth**3*x2*x3*(coew*x1*(2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x2*(-6*x1**2 + x2**2 + x3**2)))/rr**11


        df4y=df4y+(15*aearth**3*(26*sin3w*x1*x2**2 +
     -       cos3w*(-6*rr**2*x2 + 42*x1**2*x2 +
     -          6*x2*(x1**2 - x2**2)) +
     -       sin3w*x1*(6*rr**2 + 7*(-3*x1**2 + x2**2))))/rr**9 -
     -  (135*aearth**3*x2*(cos3w*
     -        (-7*x1**2*(x1**2 - 3*x2**2) +
     -          3*rr**2*(x1**2 - x2**2)) +
     -       sin3w*x1*x2*(6*rr**2 + 7*(-3*x1**2 + x2**2))))/rr**11




         df4z=df4z+  (27*aearth**3*x3*
     -     (-5*sinw*x1*x2*(x1**2 + x2**2 - 6*x3**2) +
     -       cosw*(rr**4 + 35*x1**2*x3**2 -
     -          5*rr**2*(x1**2 + x3**2))))/(2.*rr**11) -
     -  (3*aearth**3*(60*sinw*x1*x2*x3 +
     -       cosw*(-6*rr**2*x3 + 70*x1**2*x3 -
     -          10*x3*(x1**2 + x3**2))))/(2.*rr**9)



       df4z=df4z+(15*aearth**3*x3*(4*coew*x1*x3 + 4*siew*x2*x3))/rr**9+
     -  (15*aearth**3*(coew*x1*(2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x2*(-6*x1**2 + x2**2 + x3**2)))/rr**9 -
     -  (135*aearth**3*x3**2*(coew*x1*(2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x2*(-6*x1**2 + x2**2 + x3**2)))/rr**11



          df4z=df4z+   (-135*aearth**3*(cos3w*
     -        (-7*x1**2*(x1**2 - 3*x2**2) +
     -          3*rr**2*(x1**2 - x2**2)) +
     -       sin3w*x1*x2*(6*rr**2 + 7*(-3*x1**2 + x2**2)))*x3)/
     -   rr**11 + (15*aearth**3*
     -     (12*sin3w*x1*x2*x3 + 6*cos3w*(x1**2 - x2**2)*x3))/rr**9



        df5y=df5y+  (-27*aearth**3*x2*
     -     (5*cosw*x1*x2*(x1**2 + x2**2 - 6*x3**2) -
     -       sinw*(rr**4 + 35*x2**2*x3**2 -
     -          5*rr**2*(x2**2 + x3**2))))/(2.*rr**11) +
     -  (3*aearth**3*(10*cosw*x1*x2**2 +
     -       5*cosw*x1*(x1**2 + x2**2 - 6*x3**2) -
     -       sinw*(-6*rr**2*x2 + 70*x2*x3**2 -
     -          10*x2*(x2**2 + x3**2))))/(2.*rr**9)


           df5y=df5y+ (15*aearth**3*(-24*siew*x1*x2 + 10*coew*x2**2 +
     -       coew*(-2*rr**2 - 7*x1**2 + 7*x2**2))*x3)/rr**9 -
     -  (135*aearth**3*x2*x3*(coew*x2*(-2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x1*(x1**2 - 6*x2**2 + x3**2)))/rr**11


          df5y=df5y+   (15*aearth**3*(30*cos3w*x1*x2**2 +
     -       cos3w*x1*(-6*rr**2 - 7*x1**2 + 21*x2**2) +
     -       sin3w*(-6*rr**2*x2 + 14*x2**3 +
     -          6*x2*(x1**2 - x2**2) + 14*x2*(-3*x1**2 + x2**2))))/
     -   rr**9 - (135*aearth**3*x2*
     -     (cos3w*x1*x2*(-6*rr**2 - 7*x1**2 + 21*x2**2) +
     -       sin3w*(3*rr**2*(x1**2 - x2**2) +
     -          7*x2**2*(-3*x1**2 + x2**2))))/rr**11



          df5z=df5z+  (-27*aearth**3*x3*
     -     (5*cosw*x1*x2*(x1**2 + x2**2 - 6*x3**2) -
     -       sinw*(rr**4 + 35*x2**2*x3**2 -
     -          5*rr**2*(x2**2 + x3**2))))/(2.*rr**11) +
     -  (3*aearth**3*(-60*cosw*x1*x2*x3 -
     -       sinw*(-6*rr**2*x3 + 70*x2**2*x3 -
     -          10*x3*(x2**2 + x3**2))))/(2.*rr**9)


       df5z=df5z+(15*aearth**3*x3*(4*siew*x1*x3 - 4*coew*x2*x3))/rr**9+
     -  (15*aearth**3*(coew*x2*(-2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x1*(x1**2 - 6*x2**2 + x3**2)))/rr**9 -
     -  (135*aearth**3*x3**2*(coew*x2*(-2*rr**2 - 7*x1**2 + 7*x2**2) +
     -       2*siew*x1*(x1**2 - 6*x2**2 + x3**2)))/rr**11


          df5z=df5z+   (-135*aearth**3*(cos3w*x1*x2*
     -        (-6*rr**2 - 7*x1**2 + 21*x2**2) +
     -       sin3w*(3*rr**2*(x1**2 - x2**2) +
     -          7*x2**2*(-3*x1**2 + x2**2)))*x3)/rr**11 +
     -  (15*aearth**3*(-12*cos3w*x1*x2*x3 +
     -       6*sin3w*(x1**2 - x2**2)*x3))/rr**9



        df6z=df6z+ (-60*aearth**3*(cosw*x1 + sinw*x2)*x3**2)/rr**9 +
     -  (15*aearth**3*(cosw*x1 + sinw*x2)*
     -     (3*rr**2 - 7*x3**2))/(2.*rr**9) -
     -  (135*aearth**3*(cosw*x1 + sinw*x2)*x3**2*
     -     (3*rr**2 - 7*x3**2))/(2.*rr**11)


         df6z=df6z+ (-180*aearth**3*(2*siew*x1*x2 +
     -   coew*(x1**2 - x2**2))*x3)/rr**9 -
     -  (135*aearth**3*(2*siew*x1*x2 + coew*(x1**2 - x2**2))*x3*
     -     (x1**2 + x2**2 - 6*x3**2))/rr**11



          df6z=df6z+ (-105*aearth**3*(cos3w*x1**3 + 3*sin3w*x1**2*x2 -
     -       3*cos3w*x1*x2**2 - sin3w*x2**3))/rr**9 +
     -  (945*aearth**3*(cos3w*x1**3 + 3*sin3w*x1**2*x2 -
     -       3*cos3w*x1*x2**2 - sin3w*x2**3)*x3**2)/rr**11




        d4x=df4x +
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*(3.d0*(x(1)-xs(1))**2
     c   -rds**2)/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(1)-xs(1))*((xdovc-
     c    xs(5))/(exp(4.d0*dlog(rds))))

       d4y=df4y +
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(1)-xs(1))*(x(2)
     c  -xs(2))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(2)-xs(2))*((xdovc-
     c    xs(5))/(exp(4.d0*dlog(rds))))

        d4z=  df4z +
     c  (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(1)-xs(1))*(x(3)
     c  -xs(3))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(3)-xs(3))*((xdovc-
     c    xs(5))/(exp(4.d0*dlog(rds))))

      d4xdo= -crpra2*are2ms*gamma*((x(1)-xs(1))/(exp(3.d0*dlog(rds))))*
     c      (xdovc/rhdovc)
     c     -crpra2*are2ms*gamma*((1.d0)/(exp(2.d0*dlog(rds))))

      d4ydo= -crpra2*are2ms*gamma*((x(1)-xs(1))/(exp(3.d0*dlog(rds))))*
     c      (ydovc/rhdovc)

      d4zdo= -crpra2*are2ms*gamma*((x(1)-xs(1))/(exp(3.d0*dlog(rds))))*
     c      (zdovc/rhdovc)


       d5x=df4y+
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(1)-xs(1))*(x(2)
     c  -xs(2))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(1)-xs(1))*((ydovc-
     c    xs(6))/(exp(4.d0*dlog(rds))))


          d5y=df5y +
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*(3.d0*(x(2)-xs(2))**2
     c   -rds**2)/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(2)-xs(2))*((ydovc-
     c    xs(6))/(exp(4.d0*dlog(rds))))

        d5z= df5z+
     c  (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(2)-xs(2))*(x(3)
     c  -xs(3))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(3)-xs(3))*((ydovc-
     c    xs(6))/(exp(4.d0*dlog(rds))))

       d5xdo= -crpra2*are2ms*gamma*((x(2)-xs(2))/(exp(3.d0*dlog(rds))))*
     c      (xdovc/rhdovc)


      d5ydo= -crpra2*are2ms*gamma*((x(2)-xs(2))/(exp(3.d0*dlog(rds))))*
     c      (ydovc/rhdovc)
     c     -crpra2*are2ms*gamma*((1.d0)/(exp(2.d0*dlog(rds))))

      d5zdo= -crpra2*are2ms*gamma*((x(2)-xs(2))/(exp(3.d0*dlog(rds))))*
     c      (zdovc/rhdovc)






         d6x=df4z+
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(1)-xs(1))*(x(3)
     c  -xs(3))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(1)-xs(1))*((zdovc-
     c    xs(7))/(exp(4.d0*dlog(rds))))


          d6y=df5z +
     c  (crpra2*are2ms*gamma)*(rhdovc-rsdot)*3.d0*(x(2)-xs(2))*(x(3)
     c  -xs(3))/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(3)-xs(3))*((zdovc-
     c    xs(7))/(exp(4.d0*dlog(rds))))


          d6z= df6z+
     c   (crpra2*are2ms*gamma)*(rhdovc-rsdot)*(3.d0*(x(3)-xs(3))**2
     c   -rds**2)/rds**5
     c   +crpra2*are2ms*gamma*2.d0*(x(3)-xs(3))*((zdovc-
     c    xs(7))/(exp(4.d0*dlog(rds))))




       d6xdo= -crpra2*are2ms*gamma*((x(3)-xs(3))/(exp(3.d0*dlog(rds))))*
     c      (xdovc/rhdovc)


      d6ydo= -crpra2*are2ms*gamma*((x(3)-xs(3))/(exp(3.d0*dlog(rds))))*
     c      (ydovc/rhdovc)


      d6zdo= -crpra2*are2ms*gamma*((x(3)-xs(3))/(exp(3.d0*dlog(rds))))*
     c      (zdovc/rhdovc)
     c     -crpra2*are2ms*gamma*((1.d0)/(exp(2.d0*dlog(rds))))




         df4x=df4x+ (-15*aearth**4*C40*
     -    (rr**6 - 231*x1**2*x3**4 - 7*rr**4*(x1**2 + 2*x3**2) +
     -      21*rr**2*(6*x1**2*x3**2 + x3**4)))/(8.*rr**13)

       df4y=df4y+(105*aearth**4*C40*x1*x2*
     -    (rr**4 - 18*rr**2*x3**2 + 33*x3**4))/(8.*rr**13)

       df4z=df4z+ (105*aearth**4*C40*x1*
     -    (5*rr**4*x3 - 30*rr**2*x3**3 + 33*x3**5))/(8.*rr**13)

       df5y=df5y+  (-15*aearth**4*C40*
     -    (rr**6 - 231*x2**2*x3**4 - 7*rr**4*(x2**2 + 2*x3**2) +
     -      21*rr**2*(6*x2**2*x3**2 + x3**4)))/(8.*rr**13)
       
       df5z=df5z+  (105*aearth**4*C40*x2*
     -    (5*rr**4*x3 - 30*rr**2*x3**3 + 33*x3**5))/(8.*rr**13)

       df6z=df6z+ (-15*aearth**4*C40*
     -    (5*rr**6 - 105*rr**4*x3**2 + 315*rr**2*x3**4 - 231*x3**6)
     -    )/(8.*rr**13)






       g(1)=v(4)
       g(2)=v(5)
       g(3)=v(6)
       g(4)=d4x*v(1)+ d4y*v(2)+d4z*v(3)+d4xdo*v(4)+
     -   d4ydo*v(5)+d4zdo*v(6)
       g(5)=d5x*v(1)+ d5y*v(2)+d5z*v(3)+d5xdo*v(4)+
     -   d5ydo*v(5)+d5zdo*v(6)
       g(6)=d6x*v(1)+ d6y*v(2)+d6z*v(3)+d6xdo*v(4)+
     -  d6ydo*v(5)+d6zdo*v(6)



       return
       end

*********************************************************************

       subroutine sun(xs, tsid)
       implicit double precision (a-h, o-z)
       dimension xs(8)

       pi2=8.d0*datan(1.d0)
       degree= pi2/360.d0
       eps=23.4392911d0*degree
       ttime=(365.242196d0*tsid)/(36525.d0*366.242196d0*pi2)



       anmns=357.5256d0*degree+35999.049d0*degree*ttime



       slam=282.94d0*degree+anmns+(6892.d0/3600.d0)*degree*sin(anmns)
     c  + (72.d0/3600.d0)*degree*sin(2.d0*anmns)

       rs= (149.619d0-2.499d0*cos(anmns)
     c  -0.021d0*cos(2.d0*anmns))*23.71681950544094d0





       xs(1)=rs*cos(slam)
       xs(2)=rs*sin(slam)*cos(eps)
       xs(3)=rs*sin(slam)*sin(eps)
       xs(4)=rs

       slight=299792.458d0*24.d0*60.d0*60.d0*365.242196d0/
     -  (42164.1696d0*366.242196d0*pi2)

       ttdot=365.242196d0/(36525.6363d0*366.242196d0*pi2)

       anmdot=35999.049d0*degree*ttdot

          sldot=anmdot+(6892.d0/3600.d0)*degree*anmdot*cos(anmns)
     c  + (72.d0/3600.d0)*degree*2.d0*anmdot*cos(2.d0*anmns)

        rsdot= ((2.499d0*anmdot*sin(anmns)
     c  +0.021d0*2.d0*anmdot*sin(2.d0*anmns))*23.71681950544094d0)



       xs(5)=(rsdot*cos(slam)-sldot*rs*sin(slam))/slight

       xs(6)=(rsdot*sin(slam)*cos(eps)+
     -  sldot*rs*cos(slam)*cos(eps))/slight

        xs(7)=(rsdot*sin(slam)*sin(eps)+
     -  sldot*rs*cos(slam)*sin(eps))/slight

       xs(8)=rsdot/slight

       return

       end

*************************************************************************

       subroutine moon(xm, tsid)
       implicit double precision (a-h, o-z)
       dimension xm(4)

       pi2=8.d0*datan(1.d0)
       degree= pi2/360.d0
       eps=23.4392911d0*degree



       ttime=(365.242196d0*tsid)/(36525.d0*366.242196d0*pi2)



       ool0=(218.31617d0+481267.88088d0*ttime
     c  -4.06d0/3600.d0*ttime**2)*degree

       ool=(134.96292d0+477198.86753d0*ttime)*degree
       oolpr=(357.52543d0+35999.04944d0*ttime)*degree
       ff=(93.27283d0+483202.01873d0*ttime)*degree
       d=(297.85027d0+445267.11135d0*ttime)*degree



       oolam=ool0+(1.d0/3600.d0)*degree*(22640.d0*sin(ool)
     c  +769.d0*sin(2.d0*ool)-4586.d0*sin(ool-2.d0*d)
     c  +2370.d0*sin(2.d0*d)-668.d0*sin(oolpr)-412.d0*sin(2.d0*ff)
     c  -212.d0*sin(2.d0*ool-2.d0*d)-206.d0*sin(ool+oolpr-2.d0*d)
     c  +192.d0*sin(ool+2.d0*d)-165.d0*sin(oolpr-2.d0*d)
     c  +148.d0*sin(ool-oolpr)-125.d0*sin(d)
     c  -110.d0*sin(ool+oolpr)-55.d0*sin(2.d0*ff-2.d0*d))



        betmo=  (1.d0/3600.d0)*degree*(18520.d0*sin(ff+oolam-ool0
     c  +(1.d0/3600.d0)*degree*412.d0*sin(2.d0*ff)+
     c (1.d0/3600.d0)*degree*541.d0*sin(oolpr))-
     c  526.d0*sin(ff-2.d0*d)
     c  +44.d0*sin(ool+ff-2.d0*d)-31.d0*sin(-ool+ff-2.d0*d)
     c  -25.d0*sin(-2.d0*ool+ff)-23.d0*sin(oolpr+ff-2.d0*d)
     c  +21.d0*sin(-ool+ff)+11.d0*sin(-oolpr+ff-2.d0*d))


       rm=(385000.d0-20905.d0*cos(ool)-3699.d0*cos(2.d0*d-ool)
     c  -2956.d0*cos(2.d0*d)-570*cos(2.d0*ool)+246.d0*cos(2.d0*ool
     c  -2.d0*d)-205.d0*cos(oolpr-2.d0*d)-171.d0*cos(ool+2.d0*d)
     c  -152.d0*cos(ool+oolpr-2.d0*d))/42164.1696d0




       xm(1)=rm*cos(oolam)*cos(betmo)
       xm(2)=rm*sin(oolam)*cos(betmo)*cos(eps)-
     c    rm*sin(betmo)*sin(eps)
       xm(3)=rm*sin(oolam)*cos(betmo)*sin(eps)+
     c   rm*sin(betmo)*cos(eps)
       xm(4)=rm



       return

       end


*********************************************************************************
        subroutine orbital(x, orb)
       implicit double precision (a-h, o-z)
       dimension x(6), orb(7)
       dimension Vr(3), Vv(3), Vh(3), Ve(3), Vn(3), Vi(3), Vj(3), Vk(3)
       dimension VhovSh(3), VnovSn(3)

       pi2=8.d0*datan(1.d0)
       degree= pi2/360.d0
       amuE=1.d0
       tiny=1.d-8

       Vi(1)=1.d0
       Vi(2)=0.d0
       Vi(3)=0.d0

       Vj(1)=0.d0
       Vj(2)=1.d0
       Vj(3)=0.d0

       Vk(1)=0.d0
       Vk(2)=0.d0
       Vk(3)=1.d0


       Vr(1)=x(1)
       Vr(2)=x(2)
       Vr(3)=x(3)
       Sr=dsqrt(x(1)**2+x(2)**2+x(3)**2)


       Vv(1)=x(4)
       Vv(2)=x(5)
       Vv(3)=x(6)
       Sv=dsqrt(x(4)**2+x(5)**2+x(6)**2)

       call cross(Vr, Vv, Vh)
       Sh=dsqrt(Vh(1)**2+Vh(2)**2+Vh(3)**2)

       call cross(Vv, Vh, Ve)
       Ve(1)=Ve(1)/amuE-Vr(1)/Sr
       Ve(2)=Ve(2)/amuE-Vr(2)/Sr
       Ve(3)=Ve(3)/amuE-Vr(3)/Sr
       Se=dsqrt(Ve(1)**2+Ve(2)**2+Ve(3)**2)
       if (Se .eq. 0.d0) then
       Se=Se+tiny
       end if

       Sa=scprod(Vh,Vh)/(amuE*(1-Se**2))


       VhovSh(1)=Vh(1)/Sh
       VhovSh(2)=Vh(2)/Sh
       VhovSh(3)=Vh(3)/Sh
       Si=acos(scprod(Vk,VhovSh))
       if (Si .eq. 0.d0) then
       Si=Si+tiny
       end if

       call cross(Vk, Vh, Vn)
       Sn=dsqrt(scprod(Vn, Vn))
        if (Sn .eq. 0.d0) then
       Sn=Sn+tiny
       end if


       VnovSn(1)=Vn(1)/Sn
       VnovSn(2)=Vn(2)/Sn
       VnovSn(3)=Vn(3)/Sn
       SbOm=acos(scprod(Vi,VnovSn))
       SbOm=dmod(SbOm, pi2)

       if (scprod(Vn,Vj) .lt. 0.d0) then
       SbOm=pi2 -SbOm
       SbOm=dmod(SbOm, pi2)
       end if


       aux=scprod(Vn,Ve)/(Sn*Se)
       if (aux .gt.1.d0) then
       aux=1.d0
       end if
       if (aux .lt. -1.d0) then
       aux=-1.d0
       end if
       Som=acos(aux)
       Som=dmod(Som, pi2)
       if (scprod(Ve,Vk) .lt. 0.d0) then
       Som=pi2 -Som
       Som=dmod(Som, pi2)
       end if


       aux=scprod(Ve,Vr)/(Se*Sr)
       if (aux .gt.1.d0) then
       aux=1.d0
       end if
       if (aux .lt. -1.d0) then
       aux=-1.d0
       end if
       Sf=acos(aux)
       if (scprod(Vr,Vv) .lt. 0.d0) then
       Sf=pi2 -Sf
       end if

       SEE=acos((Se+cos(Sf))/(1.d0+Se*cos(Sf)))
       if ((pi2/2.d0 .lt. Sf) .and. (Sf .lt. pi2)) then
       SEE=pi2 -SEE
       end if

       SM=SEE-Se*sin(SEE)
       SMM=dsqrt(amuE/(Sa*Sa*Sa))
       St=SM/SMM

       orb(1)=Sa
       orb(2)=Se
       orb(3)=Si/degree
       orb(4)=Som/degree
       orb(5)=SbOm/degree
       orb(6)=SM/degree
       orb(7)=St

       return
       end

************************************************************************
       subroutine cross(x, y, v)
       implicit double precision (a-h ,o-z)
       dimension x(3), y(3), v(3)
       v(1)=x(2)*y(3)-y(2)*x(3)
       v(2)=x(3)*y(1)-y(3)*x(1)
       v(3)=x(1)*y(2)-y(1)*x(2)
       return
       end

***************************************************************************
       function scprod(x,y)
       implicit double precision (a-h ,o-z)
       dimension x(3), y(3)
       scprod=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
       return
       end function

***************************************************************



*********************************************************************

*                              DJULIAN

************************************************************************
*Subroutine  for the computation of the Julian Date. The algorithm is described i
*in C.D. Murray and S.F. Dermott, Solar System Dynamics, pp. 526-528

* Imput:
*      i=year, m=month, n=day
*      j=hour, k=minute, s=second (in Universal Time)
*Output:
*       dateju= Julian Date
*Local:
*     judate, l, i1, m1, ut: local parameters

      subroutine djulian(i, m, n, j, k, s, dateju)
      implicit double precision (a-h, o-z)
      if (m .le.2) then
         i1=i-1
         m1=m+12
         else
          i1=i
          m1=m
      end if
      if (i .lt. 1582) then
         l=-2
      else if (i .gt. 1582) then
         l=int(dble(i1)/400.d0)- int(dble(i1)/100.d0)
      else if (i. eq. 1582) then
           if (m .lt. 10) then
              l=-2
           else if (m .gt. 10) then
               l=int(dble(i1)/400.d0)- int(dble(i1)/100.d0)
           else if (m .eq. 10) then
                if (n. le. 4) then
                   l=-2
                else if (n.ge.15) then
                   l=int(dble(i1)/400.d0)- int(dble(i1)/100.d0)
                else
                    write(*,*) 'This date does not exist'
                    write(*,*) 'in the Gregorian calendar'
                    go to 10000
                end if
            end if
      end if
      ut=dble(j)/24.d0+dble(k)/(24.d0*60.d0)+dble(s)/(24.d0*60.d0*60.d0)
      if (i1 .ge. 0) then
         judate=int(365.25d0*dble(i1))+int( 30.6001d0*(dble(m1)+1.d0))
     -   +l+1720996+n
              dateju= dble(judate) +0.5d0+ut
      else
       judate=int(365.25d0*dble(i1))-1+int( 30.6001d0*(dble(m1)+1.d0))
     -   +l+1720996+n
        dateju=dble(judate)+0.5d0+ut
      end if
10000 continue
      return
      end






*********************************************************************
*      Parameters of the problem: masses, positions, asteroid parameters
************************************************************************************
         subroutine param(adata, are2, ipr, isrp,ieh, isun, imoon)
         implicit double precision (a-h, o-z)
         dimension adata(25)


*       Output:
*      adata(1), adata(2): Moon and Sun effects
*      adata(3):
*      adata(4)
*      adata(5)
*      adata(6)
*      adata(7)
*      adata(8)
*      adata(9)
*      adata(12)
*      adata(13)


       pi2=8.d0*datan(1.d0)
       degree= pi2/360.d0

*      Sun and Moon effects
*        amiumo=0.01230314690256d0
*        amius=333060.401628213d0

       amius=332946.0443102152d0
       amiumo=0.012300036573984492d0

        amiumo=amiumo*dble(imoon)
        amius=amius*dble(isun)


*      Earth's effects
	a2erth=0.0228823428560738d0
        aearth=0.15126910740820082d0

	C20=-1082.62602d0*exp(-6.d0*dlog(10.d0))

	C22= 1.574615325d0*exp(-6.d0*dlog(10.d0))
	S22=-0.90387d0*exp(-6.d0*dlog(10.d0))

*        C30=2.53241d-6

*	C31=2.19315d-6
*        S31=0.268087d-6

*        C32=0.30904d-6
*        S32=-0.211431d-6

*        C33=0.100583d-6
*        S33=0.197222d-6

*        C40=1.6199d-6


       C20=C20*dble(ieh)
       C22=C22*dble(ieh)
       S22=S22*dble(ieh)
       C30=C30*dble(ieh)
       C31=C31*dble(ieh)
       S31=S31*dble(ieh)
       C32=C32*dble(ieh)
       S32=S32*dble(ieh)
       C33=C33*dble(ieh)
       S33=S33*dble(ieh)
       C40=C40*dble(ieh)

*********************       SRP             ******************

*     For SRP we can use one of the two approaches and consider either  A/m or either beta as parameters.
*     Now, it is considered the second approach.

*      In the first case, in the units used here, I considered the parameter A/m and the constant
*      crpra2, where
*       crpra2=Cr*Pr*a_s*a_s, with Cr the reflectivity coefficient,
*       Pr the radiation pressure for an object located at a_s=1 au (see Valk et al., Adv. Space Research, 2009).


*     In the second formulation, I used the parameter beta (denoted by are2ms) and
*     the constant \mu_s (denoted by  crpra2).


*      First approach:

*	crpra2=256.022356929d0
	
	crpra2=256.09468404128717d0
	
        are2ms=are2



*      Second approach(\mu_s=crpra2, are2ms=beta):

*      crpra2= 333060.401628213d0
*       are2ms=0.01

*      The parameter cont takes the values 0 and 1. The value 0 means that
*      SRP is killed and 1 means that SRP is taked into account.


       cont=dble(isrp)
***************************************************************************


********************** PR drag ******************************************
*      GAMMA is a constant describing the Poynting-Robertson drag.
*      Delete the comment in front of gamma to take this effect into account.

       gamma=1.d0*dble(ipr)

**************************************************************************



C      SLIGHT is the speed of light in our units

       slight=299792.458d0*24.d0*60.d0*60.d0*365.242196d0/
     -  (42164.1696d0*366.242196d0*pi2)




         adata(1)= amiumo
         adata(2)= amius

         adata(3)= a2erth
         adata(4)=aearth

         adata(5)=C20
         adata(6)=C22
         adata(7)=S22
         adata(8)=C30
         adata(9)=C31
         adata(10)=S31
         adata(11)=C32
         adata(12)=S32
         adata(13)=C33
         adata(14)=S33

         adata(15)= are2ms
         adata(16)= crpra2
         adata(17)= cont
         adata(18)= gamma
         adata(19)= slight

         adata(20)= C40
         adata(21)= 0.d0
         adata(22)= 0.d0
         adata(23)= 0.d0
         adata(24)= 0.d0
         adata(25)= 0.d0

         return
         end
