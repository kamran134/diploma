c
c     DDOPRI5 - вариант с двойной точностью подпрограммы
c     DOPRI5 из Э.Хайрера, С.Нерсетта и Г.Ваннера
c     "Решение обыкновенных дифференциальных уравнений.
c     Нежесткие задачи", М., Мир, 1990, 512 с.,
c     с. 452-454.
c
c     Явный метод Рунге-Кутта (формулы Дорманда-Принса)
c     5 (4) порядка. 
c
c     Рекомендуемый диапазон точности ? - ?.
c     При EPS=1.0d-10 для гармонического осциллятора на интервале
c     интегрирования [ 0.0, 90.0 ] глобальная ошибка равна 1.0d-9
c     и требуется почти 3000 шагов.
c
c     Отлажено на примере гармонического осциллятра.
c
c     Внимание! В алгоритме стоит проверка на число сделанных
c     шагов. Максимально допустимое число шагов NMAX установлено
c     на 3000. При превышении этого числа - аварийный выход с
c     диагностической печатью.
c
      implicit real*8(a-h,o-z)
      external fcn
      real*8 y(4),eps
      common/stat/nfcn,nstep,naccpt,nrejct
      eps=1.d-11
c      write(*,*)' a='
c      read(*,*) a
      x=0.
      y(1)=0.0d0
      y(2)=0.0d0
      y(3)=1
      y(4)=1
      xend=1.0d0
      hmax=1.0d0
      h=0.5d0
      call ddopri5(4,fcn,x,y,xend,eps,hmax,h)
      print *,nfcn,nstep,naccpt,nrejct
      print *,' '
      stop ' '
      end

      subroutine fcn(n,x,y,f)
      implicit real*8(a-h,o-z)
      real*8 y(n),f(n)
C 1=x, 2=dx, 3=lambda
      f(1)=y(2);
      f(2)=-y(1)+y(4);
      f(3)=-y(4);
      f(4)=y(3);
      return
      end

      subroutine solous(nrpnts,x,y,n)
      implicit real*8(a-h,o-z)
      real*8 y(n)
C      common/sol_1/twopi
C 1=x, 2=y, 3=u, 4=v, 5=px, 6=py, 7=pu, 8=pv
      print *,x,y(1),y(2),y(3),y(4)
      return
      end

      subroutine ddopri5(n,fcn,x,y,xend,eps,hmax,h)
      implicit real*8(a-h,o-z)
      real*8 k1(10),k2(10),k3(10),k4(10),k5(10),y1(10),y(n)
      logical reject
      common/stat/nfcn,nstep,naccpt,nrejct
      data nmax/3000/,uround/2.2205d-16/
      posneg=dsign(1.d0,xend-x)
      hmax=dabs(hmax)
      h=dmin1(dmax1(1.d-4,dabs(h)),hmax)
      h=dsign(h,posneg)
      eps=dmax1(eps,7.d0*uround)
      reject=.false.
      naccpt=0
      nrejct=0
      nfcn=1
      nstep=0
      call solous(naccpt+1,x,y,n)
      call fcn(n,x,y,k1)
1     if(nstep.gt.nmax.or.x+.1d0*h.eq.x)goto 79
      if((x-xend)*posneg+uround.gt.0.d0)return
      if((x+h-xend)*posneg.gt.0.d0)h=xend-x
      nstep=nstep+1
      do 22 i=1,n
22    y1(i)=y(i)+h*.2d0*k1(i)
      call fcn(n,x+.2d0*h,y1,k2)
      do 23 i=1,n
23    y1(i)=y(i)+h*((3.d0/40.d0)*k1(i)+(9.d0/40.d0)*k2(i))
      call fcn(n,x+.3d0*h,y1,k3)
      do 24 i=1,n
24    y1(i)=y(i)+h*((44.d0/45.d0)*k1(i)-(56.d0/15.d0)*k2(i)+
     +          (32.d0/9.d0)*k3(i))
      call fcn(n,x+.8d0*h,y1,k4)
      do 25 i=1,n
25    y1(i)=y(i)+h*((19372.d0/6561.d0)*k1(i)-
     &  (25360.d0/2187.d0)*k2(i)+
     +  (64448.d0/6561.d0)*k3(i)-(212.d0/729.d0)*k4(i))
      call fcn(n,x+(8.d0/9.d0)*h,y1,k5)
      do 26 i=1,n
26    y1(i)=y(i)+h*((9017.d0/3168.d0)*k1(i)-(355.d0/33.d0)*k2(i)+
     +    (46732.d0/5247.d0)*k3(i)+(49.d0/176.d0)*k4(i)-
     &    (5103.d0/18656.d0)*k5(i))
      xph=x+h
      call fcn(n,xph,y1,k2)
      do 27 i=1,n
27    y1(i)=y(i)+h*((35.d0/384.d0)*k1(i)+(500.d0/1113.d0)*k3(i)+
     +    (125.d0/192.d0)*k4(i)-(2187.d0/6784.d0)*k5(i)+
     +    (11.d0/84.d0)*k2(i))
      do 61 i=1,n
61    k2(i)=(71.d0/57600.d0)*k1(i)-(71.d0/16695.d0)*k3(i)+
     +    (71.d0/1920.d0)*k4(i)-(17253.d0/339200.d0)*k5(i)+
     +    (22.d0/525.d0)*k2(i)
      call fcn(n,xph,y1,k3)
      do 28 i=1,n
28    k4(i)=(k2(i)-(1.d0/40.d0)*k3(i))*h
      nfcn=nfcn+6
      err=0.
      do 41 i=1,n
      denom=dmax1(1.d-5,dabs(y1(i)),dabs(y(i)),2.d0*uround/eps)
41    err=err+(k4(i)/denom)**2
      err=dsqrt(err/dble(n))
      fac=dmax1(.1d0,dmin1(5.d0,(err/eps)**(0.2d0)/.9d0))
      hnew=h/fac
      if(err.le.eps)then
              naccpt=naccpt+1
              do 44 i=1,n
              k1(i)=k3(i)
44            y(i)=y1(i)
              x=xph
              call solous(naccpt+1,x,y,n)
              if(dabs(hnew).gt.hmax)hnew=posneg*hmax
              if(reject)hnew=posneg*dmin1(dabs(hnew),dabs(h))
              reject=.false.
      else
              reject=.true.
              if(naccpt.ge.1)nrejct=nrejct+1
      end if
      h=hnew
      goto 1
79    write(6,979)x
979   format(' Exit of DDOPRI5 at x=',e11.4)
      return
      end
