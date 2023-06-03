************************************************************************
*                      " KdVB 2022 "                       13.12.2022  *
*                                                                      *
* Расчёт явным методом МакКормака уравнения kdvb с учётом "пелены"     *
*    шириной [x1, x2], x2=x1+hx                                        *
* Сквозной счёт.                                                       *
* Метод коррекции потоков (метод FCT).                                 *
*                                                                      *
*     Intrinsic function: MAX, MIN, DEXP, DSQRT, IDSIGN, IDINT         *
*                                                                      *
*     Обозначения:                                                     *
* cNV - number of variant                                              *
* i_shape=0 - прямоугольник                                            *
*   h_max - высота, [-jx, jx] - ширина (кол-во узлов)                  *                                 *
* i_shape=1 - солитон		                                               *
* i_shape=2 - ступенька	                                               *
* [x1, x1+hx] - ширина "пелены", x1 - начало (x1>j0), hx - ширина      *
* ta - step over time                                                  *
* dx - шаг по координате x                                             *
* t_stop - максимальное время вычисления                               *
* istep - number of steps				                                       *
* i_corr=0 - коррекция потоков                                         *
************************************************************************
      IMPLICIT REAL*8 (a-h,o-z)

      ALLOCATABLE u0(:), us(:), u1(:),
     1            atx(:), atx2(:), btx(:), btx2(:),
     2            cu_d(:), cu_ad(:), u_d(:), u_ad(:), d_u(:),
     3            S(:), h(:), fj(:),fj2(:)

c      CHARACTER*7 cNV(50)
      CHARACTER cNV*4, cSer*2, cSer0*2, fil0*7, fil*7

      DIMENSION open_time(50)

      NAMELIST /solo/j2,delta_j0,t_stop,iout,ta,dx,i_open0,
     1               alf1,alf2,bet1,bet2,gam,ck,x1,hx
      NAMELIST /correct/i_corr,cu_d0,cu_d1,cu_ad0,cu_ad1
      NAMELIST /change/open_time, cNV, n_Ser0
      NAMELIST /shape/i_shape,h_max,jx,eps

      OPEN(20,FILE='..\input.dat',STATUS='old')
      READ(20,solo)
      READ(20,correct)
      READ(20,shape)
      READ(20,change)

      OPEN(21,FILE='..\DATA.IN\'//'i'//cNV//'.dat')

      WRITE(21,solo)
      WRITE(21,correct)
      WRITE(21,shape)
      WRITE(21,change)
      CLOSE(21)

      WRITE(cSER0,'(i2.2)') n_SER0   ! перевод целого в текст

      OPEN(35,FILE='..\DATA.DAT\'//'g'//cNV//'-'//cSER0//'.dat')
      OPEN(45,FILE='..\DATA.OUT\'//'o'//cNV//'.dat')

      WRITE(*,3) cNV

    3 FORMAT(/'    Open file = ',a4/)

      WRITE(45,3) cNV

      j0=IDINT(delta_j0*j2/2)
      IF(j0.GE.j2-2) STOP 'Process terminated, delta_j0 is invalid'
      IF(jx.GE.j0-2) STOP 'Process terminated, jx is invalid'

      m1=IDINT(x1/dx)+j0
      x2=x1+hx
      m2=IDINT(x2/dx)+j0
      IF((m1-j0)*dx.LT.x1) m1=m1+1
      IF((m2-j0)*dx.GT.x2) m2=m2-1

      ALLOCATE(u0(j2))
      ALLOCATE(us(j2))
      ALLOCATE(u1(j2))
      ALLOCATE(h(j2))
      ALLOCATE(atx(j2))
      ALLOCATE(atx2(j2))
      ALLOCATE(btx(j2))
      ALLOCATE(btx2(j2))

      ALLOCATE(S(j2))
      ALLOCATE(fj(j2))
      ALLOCATE(fj2(j2))

      ALLOCATE(cu_d(j2))
      ALLOCATE(cu_ad(j2))
      ALLOCATE(u_d(j2))
      ALLOCATE(u_ad(j2))
      ALLOCATE(d_u(j2))

      x_min=(1-j0)*dx
      x_max=(j2-j0)*dx

*     filling of arrays

      DO j=1,j2
        u0(j)=0.0
        us(j)=0.0
        u1(j)=0.0
        h(j)=0.0
        S(j)=0.0
        fj(j)=0.0
        fj2(j)=0.0
      END DO

      j9=j2-1
      j7=j2-2
      j6=j2-3
      dt2=ta/2
      dt4=ta/4

**********     Начальное распределение      **************

      SELECT CASE(i_shape)

        CASE(0)             ! вид прямоугольника (i_shape=0)
          DO j=j0-jx,j0+jx
            h(j)=h_max
          END DO

        CASE(1)             ! тип солитона (i_shape=1)
          IF(bet1.EQ.0) THEN
            STOP 'Initial value beta_1 for shape soliton equal 0!'
          END IF
          del=DSQRT(eps/(12.0*bet1))
          DO j=3,j7
            x=2*del*(j-j0)*dx
            IF(x.GE.0) THEN
              h(j)=4.0*eps*DEXP(-x)/(1.0+DEXP(-x))**2
            ELSE
              h(j)=4.0*eps*DEXP(x)/(1.0+DEXP(x))**2
            END IF
          END DO
          DO j=m1,j7
            h(j)=0.0
          END DO

        CASE(2)          ! скачок (i_shape=2)
          DO j=1,j0
            h(j)=1.0
          END DO
          DO j=j0+1,j2
           h(j)=0
          END DO
          us(1)=1.0
          us(2)=1.0

        CASE(10)       ! Без распределения
          CONTINUE

      END SELECT

      IF(t_stop.EQ.0.0) THEN       ! вывод начального распределения
        DO j=1,j2
          x=(j-j0)*dx
          WRITE(35,410) x,h(j)
        END DO
        WRITE(*,*)'Received the initial distribution shape = ',i_shape
        STOP 'Initial values is done'
      END IF

**     Запись начального распределения     ***

      DO j= 3,j7
        u0(j)=h(j)
        u1(j)=h(j)
      END DO

*     Начальная "масса"

      V0=0.
      DO j=3,j7
        V0=V0+u1(j)
      END DO
      V0=dx*V0

      WRITE(*,5) V0,x1,x2
    5 FORMAT(/' Initial mass = ',1pe10.4,/'  The width of the shroud'
     1       '  x1 = ',1pd10.4,'  x2 = '1pd10.4,/)
      WRITE(45,5) V0,x1,x2

*     Постоянные

      cS=gam*ta
      cS2=gam*ta/2
      edt=DEXP(-ck*ta)
      edt2=DEXP(-ck*ta/2)
      tx=ta/(2*dx)
      tx2=ta/(4*dx)
      at=ta/dx**2
      at2=ta/(2*dx**2)
      bt=ta/(2*dx**3)
      bt2=ta/(4*dx**3)

      DO j=1,j2
        IF(j.LT.m1.OR.j.GT.m2) THEN
          atx(j)=alf1*at
          atx2(j)=alf1*at2
          btx(j)=bet1*bt
          btx2(j)=bet1*bt2
        ELSE
          atx(j)=alf2*at
          atx2(j)=alf2*at2
          btx(j)=bet2*bt
          btx2(j)=bet2*bt2
        END IF
      END DO

      istep=0
      iout0=0
      i_open=1
      n_SER=n_SER0

      mac=-1     ! переключатель для МакКормака

*********************************************************
*
*       Расчёт по МакКормаку
*
*********************************************************

   10 istep=istep+1
      time=istep*ta
      iout0=iout0+1

      mac=-mac
      IF(mac.LT.0) THEN

******                                                       ***********
*               Первый шаг. Разности "назад"    mac=-1                       *
******                                                       ***********

        DO j=3,j7
          jl=j-1
          jl2=j-2
          jr=j+1
          jr2=j+2

          us(j)=u0(j)-tx*(u0(j)**2-u0(jl)**2)-
     1                atx(j)*(u0(jr)-2*u0(j)+u0(jl))-
     2                btx(j)*(u0(jr2)-2*u0(jr)+2*u0(jl)-u0(jl2))+
     3                fj(j)
        END DO

******                                                       ***********
*                Второй шаг. Разности "вперёд"    mac=-1                     *
******                                                       ***********

        DO j=m1,m2    ! Вычисление источника в промежуточной точке
          S2=S(j)+dt4*(edt2*u0(j)+us(j))
          fj2(j)=-cS2*(us(j)-ck*S2)
        END DO

        DO j=3,j7
          jl=j-1
          jl2=j-2
          jr=j+1
          jr2=j+2

          u1(j)=0.5*(u0(j)+us(j))-
     1          tx2*(us(jr)**2-us(j)**2)-
     2          atx2(j)*(us(jr)-2*us(j)+us(jl))-
     3          btx2(j)*(us(jr2)-2*us(jr)+2*us(jl)-us(jl2))+
     4          fj2(j)
        END DO

        GO TO 23      ! Обход случая mac=+1
      END IF

******                                                       ***********
*               Первый шаг. Разности "вперёд"    mac=+1                *
******                                                       ***********

      DO j=3,j7

        jl=j-1
        jl2=j-2
        jr=j+1
        jr2=j+2

        us(j)=u0(j)-tx*(u0(jr)**2-u0(j)**2)-
     1              atx(j)*(u0(jr)-2*u0(j)+u0(jl))-
     2              btx(j)*(u0(jr2)-2*u0(jr)+2*u0(jl)-u0(jl2))+
     3              fj(j)
      END DO

******                                                       ***********
*                Второй шаг. Разности "назад"      mac=+1              *
******                                                       ***********

      DO j=m1,m2    ! Вычисление источника в промежуточной точке
        S2=S(j)+dt4*(edt2*u0(j)+us(j))
        fj2(j)=-cS2*(us(j)-ck*S2)
      END DO

      DO j=3,j7
        jl=j-1
        jl2=j-2
        jr=j+1
        jr2=j+2

        u1(j)=0.5*(u0(j)+us(j))-
     1        tx2*(us(j)**2-us(jl)**2)-
     2        atx2(j)*(us(jr)-2*us(j)+us(jl))-
     3        btx2(j)*(us(jr2)-2*us(jr)+2*us(jl)-us(jl2))+
     4        fj2(j)
      END DO

   23 CONTINUE     ! Переход после случая для МакКормака mac=-1

      IF(i_corr.EQ.1) GOTO 27   ! обход коррекции  =================

****************    КОРРЕКЦИЯ ПОТОКОВ    **********************
*
*       Обозначения
* ur - полусумма скоростей
* cu_d - диффузионные коэффициенты
* cu_ad - антидиффузионные коэффициенты
*
************************************************************

*     Вычисление коэффициентов

      DO j=2,j9
        ur=u0(j)+u0(j+1)
        txx=(tx*ur)**2
        cu_d(j)=cu_d0+cu_d1*txx
        cu_ad(j)=cu_ad0+cu_ad1*txx
      END DO

      DO j=2,j9
        jr=j+1
        jl=j-1

*       1. Вычисление диффузионных потоков

        u_d(j)=cu_d(j)*(u0(jr)-u0(j))

*       2. Вычисление антидиффузионных потоков

        u_ad(j)=cu_ad(j)*(u1(jr)-u1(j))

      END DO

*     3. Численная диффузия решения

      DO j=2,j9
        u1(j)=u1(j)+u_d(j)-u_d(j-1)
      END DO

*     4. Вычисление первых разностей

      DO j=1,j9
        d_u(j)=u1(j+1)-u1(j)
      END DO

*     5. Ограничение антидиффузионных потоков.

      DO j=3,j7
        SI=DSIGN(1.0d0,u_ad(j))
        u_ad(j)=SI*MAX(0.,MIN(SI*d_u(j-1),ABS(u_ad(j)),SI*d_u(j+1)))
      END DO

*     6. Окончательное решение.

      DO j=3,j7
        u1(j)=u1(j)-u_ad(j)+u_ad(j-1)
      END DO

      IF(i_shape.EQ.2) THEN
        u1(2)=1.0
        u1(j9)=0.0
      ELSE IF(i_shape.EQ.3) THEN
        u1(2)=2.0
        u1(j9)=0.0
      END IF

   27 continue    !  Точка обхода коррекции потоков

****     вычисление интеграла S   ******

      DO j=m1,m2
        S(j)=S(j)+dt2*(edt*u0(j)+u1(j))
        fj(j)=-cS*(u1(j)-ck*S(j))
      END DO

*     Перезапись результатов

      DO j=1,j2
        u0(j)=u1(j)
      END DO

  410 FORMAT(1pd12.3,e12.3)

      IF(iout0.EQ.iout) THEN
        iout0=0
        u_min=u1(1)
        u_max=u1(1)

        DO j=1,j2
          u_max=MAX(u1(j),u_max)
          u_min=MIN(u1(j),u_min)
        END DO

        WRITE(*,70) istep,time,u_min,u_max,u1(m1),u1(m2)
        WRITE(45,70) istep,time,u_min,u_max,u1(m1),u1(m2)
      END IF

   70 FORMAT('step=',I8,' time=',f8.2,
     1       ' Min=',1Pd9.2,' Max=',d9.2,
     2       ' u_x1=',d9.2,' u_x2=',d9.2)

      IF(time.GE.open_time(i_open).AND.i_open.LT.i_open0) THEN

        DO j=1,j2      ! Запись результатов в файл
          x=(j-j0)*dx
          WRITE(35,410) x,u1(j)
        END DO

        IF(time.GE.t_stop) GO TO 80   !  окончание счёта

        i_open=i_open+1
        n_Ser=n_Ser+1
        IF(n_Ser.GT.99) STOP 'overflow n_Ser'

        WRITE(cSer,'(i2.2)') n_Ser   ! перевод целого в текст

        CLOSE(35)
        OPEN(35,FILE='..\DATA.DAT\'//'g'//cNV//'-'//cSER//'.dat')

         fil0=cNV//'-'//cSer0
         fil=cNV//'-'//cSer

        WRITE(*,75) fil0,fil
        WRITE(45,75) fil0,fil
        cSer0=cSer

      END IF
   75 FORMAT(/'      Close file = ',a7,'   Open file = ',a7/)

      GO TO 10

   80 CONTINUE

      DEALLOCATE(u0)
      DEALLOCATE(us)
      DEALLOCATE(u1)
      DEALLOCATE(atx)
      DEALLOCATE(atx2)
      DEALLOCATE(btx)
      DEALLOCATE(btx2)
      DEALLOCATE(h)
      DEALLOCATE(S)
      DEALLOCATE(fj)
      DEALLOCATE(fj2)
      DEALLOCATE(cu_d)
      DEALLOCATE(cu_ad)
      DEALLOCATE(u_d)
      DEALLOCATE(u_ad)
      DEALLOCATE(d_u)

      CLOSE(35)
      CLOSE(40)
      CLOSE(45)

      END
