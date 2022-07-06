IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION X(1:80000),Y(1:80000),Z(1:80000)
DIMENSION C1IX(1:80000),C1IY(1:80000),C1IZ(1:80000)
DIMENSION C4IX(1:80000),C4IY(1:80000),C4IZ(1:80000)
DIMENSION O51IX(1:80000),O51IY(1:80000),O51IZ(1:80000)
DIMENSION O52IX(1:80000),O52IY(1:80000),O52IZ(1:80000)
DIMENSION O51MX(1:80000),O51MY(1:80000),O51MZ(1:80000)
DIMENSION O52MX(1:80000),O52MY(1:80000),O52MZ(1:80000)
DIMENSION O51FX(1:80000),O51FY(1:80000),O51FZ(1:80000)
DIMENSION O52FX(1:80000),O52FY(1:80000),O52FZ(1:80000)
DIMENSION C1MX(1:80000),C1MY(1:80000),C1MZ(1:80000)
DIMENSION C4MX(1:80000),C4MY(1:80000),C4MZ(1:80000)
DIMENSION C1FX(1:80000),C1FY(1:80000),C1FZ(1:80000)
DIMENSION C4FX(1:80000),C4FY(1:80000),C4FZ(1:80000)
DIMENSION X2(1:80000),Y2(1:80000),Z2(1:80000)
DIMENSION DIIMIN(1:80000),DIMMIN(1:80000), DIFMIN(1:80000)
DIMENSION ICMIN(1:80000)
DIMENSION LABEL(1:80000)
DIMENSION LABEL2(1:80000)
DIMENSION ICIMIN(1:80000)
DIMENSION ICMMIN(1:80000)
DIMENSION ICFMIN(1:80000)
DIMENSION IPASS(1:80000)
CHARACTER*40 OPK1,OPK2,OPK3
ntime=1000
ICOUNT=0
open(unit=10, file="lammpstrj.xyz")
C NUMERO DE STEPS DE EQUILIBRACION: AHORA 1000      
do iter=1,1000
read (10,*) NATOMS
read (10,*)
  do i=1,NATOMS
   read (10,*) 
  enddo
enddo
WRITE(6,*) "COUNT ","C ATOM CLOSE", "CO2 MOLEC", "ITER","GO/BACK"
iter=iter-1

c PRIMERA LECTURA
read (10,*) NATOMS
read (10,*) OPK1,OPK2,OPK3
c      write(6,*) OPK1,OPK2,OPK3
K=1
KK=1
KK1=1
KK2=1
ISW=0
do i=1,NATOMS
read (10,*) LABEL(i),X(i),Y(i),Z(i)
IF (LABEL(i).EQ.1) THEN
    C1IX(K)=X(i)
    C1IY(K)=Y(i)
    C1IZ(K)=Z(i)
    K=K+1
endif
IF (LABEL(i).EQ.4) THEN
    C4IX(KK)=X(i)
    C4IY(KK)=Y(i)
    C4IZ(KK)=Z(i)
    KK=KK+1
endif
IF (LABEL(i).EQ.5) THEN
   IF (ISW.EQ.0) THEN
    ISW=1
    O51IX(KK1)=X(i)
    O51IY(KK1)=Y(i)
    O51IZ(KK1)=Z(i)
    KK1=KK1+1
   ELSE
    O52IX(KK2)=X(i)
    O52IY(KK2)=Y(i)
    O52IZ(KK2)=Z(i)
    KK2=KK2+1
    ISW=0
endif
endif
enddo

c      write(6,*) "CARBONOS"
c      do i=1,K-1
c       write(6,*) i,C1IZ(i)
c      enddo
c      write(6,*) "CO2"
c      do i=1,KK-1
c       write(6,*) i,O4IZ(i)
c      enddo

c      write(6,*) iter
c      write(6,*) LABEL(4096),X(4096),Y(4096),Z(4096)


do i=1,KK-1
 DISTM=20000.0

 do j=1,K-1
  DX=(C1IX(j)-C4IX(i))**2
  DY=(C1IY(j)-C4IY(i))**2
  DZ=(C1IZ(j)-C4IZ(i))**2
  DIST=DSQRT(DX+DY+DZ)
  IF (DIST.LT.DISTM) THEN
    DISTM=DIST
    DIIMIN(i)=DISTM
    ICIMIN(i)=j
  ENDIF
  enddo
  enddo

C SEGUNDA LECTURA
   
read (10,*) 
read (10,*) OPK1,OPK2,OPK3
c      write(6,*) OPK1,OPK2,OPK3
K=1
KK=1
KK1=1
KK2=1
ISW=0
do i=1,NATOMS
read (10,*) LABEL(i),X(i),Y(i),Z(i)
IF (LABEL(i).EQ.1) THEN
    C1MX(K)=X(i)
    C1MY(K)=Y(i)
    C1MZ(K)=Z(i)
    K=K+1
endif
IF (LABEL(i).EQ.4) THEN
    C4MX(KK)=X(i)
    C4MY(KK)=Y(i)
    C4MZ(KK)=Z(i)
    KK=KK+1
endif
IF (LABEL(i).EQ.5) THEN
   IF (ISW.EQ.0) THEN
    ISW=1
    O51MX(KK1)=X(i)
    O51MY(KK1)=Y(i)
    O51MZ(KK1)=Z(i)
    KK1=KK1+1
   ELSE
    O52MX(KK2)=X(i)
    O52MY(KK2)=Y(i)
    O52MZ(KK2)=Z(i)
    KK2=KK2+1
    ISW=0
endif
endif
enddo

c      write(6,*) "CARBONOS"
c      do i=1,K-1
c       write(6,*) i,C1FZ(i)
c      enddo
c      write(6,*) "CO2"
c      do i=1,KK-1
c       write(6,*) i,O4FZ(i)
c      enddo

c      write(6,*) iter
c      write(6,*) LABEL(4096),X(4096),Y(4096),Z(4096)


do i=1,KK-1
 DISTM=20000.0

 do j=1,K-1
  DX=(C1MX(j)-C4MX(i))**2
  DY=(C1MY(j)-C4MY(i))**2
  DZ=(C1MZ(j)-C4MZ(i))**2
  DIST=DSQRT(DX+DY+DZ)
  IF (DIST.LT.DISTM) THEN
    DISTM=DIST
    DIMMIN(i)=DISTM
    ICMMIN(i)=j
  ENDIF
  enddo
  enddo

c      write(6,*) "+cercanos"
c      do i=1,KK-1
c       write(6,*) ICIMIN(i),ICFMIN(i),DIIMIN(i),DIFMIN(i)
c      enddo

C TERCERA LECTURA

read (10,*) 
read (10,*) OPK1,OPK2,OPK3
c      write(6,*) OPK1,OPK2,OPK3
K=1
KK=1
KK1=1
KK2=1
KK=1
ISW=0
do i=1,NATOMS
read (10,*) LABEL(i),X(i),Y(i),Z(i)
IF (LABEL(i).EQ.1) THEN
    C1FX(K)=X(i)
    C1FY(K)=Y(i)
    C1FZ(K)=Z(i)
    K=K+1
endif
IF (LABEL(i).EQ.4) THEN
    C4FX(KK)=X(i)
    C4FY(KK)=Y(i)
    C4FZ(KK)=Z(i)
    KK=KK+1
endif
    
IF (LABEL(i).EQ.5) THEN
   IF (ISW.EQ.0) THEN
    ISW=1
    O51FX(KK1)=X(i)
    O51FY(KK1)=Y(i)
    O51FZ(KK1)=Z(i)
    KK1=KK1+1
   ELSE
    O52FX(KK2)=X(i)
    O52FY(KK2)=Y(i)
    O52FZ(KK2)=Z(i)
    KK2=KK2+1
    ISW=0
endif
endif
enddo

c      write(6,*) "CARBONOS"
c      do i=1,K-1
c       write(6,*) i,C1FZ(i)
c      enddo
c      write(6,*) "CO2"
c      do i=1,KK-1
c       write(6,*) i,O4FZ(i)
c      enddo

c      write(6,*) iter
c      write(6,*) LABEL(4096),X(4096),Y(4096),Z(4096)


do i=1,KK-1
 DISTM=20000.0

 do j=1,K-1
  DX=(C1FX(j)-C4FX(i))**2
  DY=(C1FY(j)-C4FY(i))**2
  DZ=(C1FZ(j)-C4FZ(i))**2
  DIST=DSQRT(DX+DY+DZ)
  IF (DIST.LT.DISTM) THEN
    DISTM=DIST
    DIFMIN(i)=DISTM
    ICFMIN(i)=j
  ENDIF
  enddo
  enddo


c      write(6,*) "+cercanos"
c      do i=1,KK-1
c       write(6,*) ICIMIN(i),ICMMIN(i),ICFMIN(i)
c      enddo



C COMPROBACION

  IOQ=0

  do i=1,80000
  IPASS(i)=0
  enddo
C STEPS DE PRODUCCION!
  do while (IOQ.LT.4998)  

DO i=1,KK-1

IF (IPASS(i).EQ.0) THEN


IF ((((O51MZ(i).LT.C1MZ(ICMMIN(i)))).AND.
$   ((O52MZ(i).LT.C1MZ(ICMMIN(i)))).AND.
$   ((C4MZ(i).LT.C1MZ(ICMMIN(i))))) .AND. 
$   ((O51IZ(i).GT.C1IZ(ICIMIN(i))).OR.
$   (O52IZ(i).GT.C1IZ(ICIMIN(i))).OR.
$   (C4IZ(i).GT.C1IZ(ICIMIN(i))))) THEN

ICOUNT=ICOUNT+1
c      WRITE(6,*) "PASA: ",i

IPASS(i)=i

WRITE(6,*) ICOUNT,ICIMIN(i),i, IOQ+1001, "+"
 ENDIF
ELSE

IF ((((O51MZ(i).GT.C1MZ(ICMMIN(i)))).AND.
$   ((O52MZ(i).GT.C1MZ(ICMMIN(i)))).AND.
$   ((C4MZ(i).GT.C1MZ(ICMMIN(i))))) .AND. 
$   ((O51IZ(i).LT.C1IZ(ICIMIN(i))).OR.
$   (O52IZ(i).LT.C1IZ(ICIMIN(i))).OR.
$   (C4IZ(i).LT.C1IZ(ICIMIN(i))))) THEN

ICOUNT=ICOUNT-1
IPASS(i)=0
c      WRITE(6,*) "VUELVE ",i
WRITE(6,*)  ICOUNT,ICIMIN(i),i, IOQ+1001, "-"
ENDIF


ENDIF


ENDDO


DO i=1,KK-1
DIIMIN(i)=DIMMIN(i)
ICIMIN(i)=ICMMIN(i)
DIMMIN(i)=DIFMIN(i)
ICMMIN(i)=ICFMIN(i)
O51IZ(i)=O51MZ(i)
O51IX(i)=O51MX(i)
O51IY(i)=O51MY(i)
O52IZ(i)=O52MZ(i)
O52IX(i)=O52MX(i)
O52IY(i)=O52MY(i)
O51MZ(i)=O51FZ(i)
O51MX(i)=O51FX(i)
O51MY(i)=O51FY(i)
O52MZ(i)=O52FZ(i)
O52MX(i)=O52FX(i)
O52MY(i)=O52FY(i)
C4IX(i)=C4MX(i)
C4IY(i)=C4MY(i)
C4IZ(i)=C4MZ(i)
C4MX(i)=C4FX(i)
C4MY(i)=C4FY(i)
C4MZ(i)=C4FZ(i)
ENDDO
DO i=1,K-1
C1IZ(i)=C1MZ(i)
C1IY(i)=C1MY(i)
C1IX(i)=C1MX(i)
C1MZ(i)=C1FZ(i)
C1MY(i)=C1FY(i)
C1MX(i)=C1FX(i)
ENDDO

C NEW LECTURA

read (10,*) 
read (10,*) OPK1,OPK2,OPK3
c      write(6,*) OPK1,OPK2,OPK3
K=1
KK=1
KK1=1
KK2=1
ISW=0
do i=1,NATOMS
read (10,*) LABEL(i),X(i),Y(i),Z(i)
IF (LABEL(i).EQ.1) THEN
    C1FX(K)=X(i)
    C1FY(K)=Y(i)
    C1FZ(K)=Z(i)
    K=K+1
endif
IF (LABEL(i).EQ.4) THEN
    C4FX(KK)=X(i)
    C4FY(KK)=Y(i)
    C4FZ(KK)=Z(i)
    KK=KK+1
endif
IF (LABEL(i).EQ.5) THEN
   IF (ISW.EQ.0) THEN
    ISW=1
    O51FX(KK1)=X(i)
    O51FY(KK1)=Y(i)
    O51FZ(KK1)=Z(i)
    KK1=KK1+1
   ELSE
    O52FX(KK2)=X(i)
    O52FY(KK2)=Y(i)
    O52FZ(KK2)=Z(i)
    KK2=KK2+1
    ISW=0
endif
endif
enddo

c      write(6,*) "CARBONOS"
c      do i=1,K-1
c       write(6,*) i,C1FZ(i)
c      enddo
c      write(6,*) "CO2"
c      do i=1,KK-1
c       write(6,*) i,O4FZ(i)
c      enddo

c      write(6,*) iter
c      write(6,*) LABEL(4096),X(4096),Y(4096),Z(4096)


do i=1,KK-1
 DISTM=20000.0

 do j=1,K-1
  DX=(C1FX(j)-C4FX(i))**2
  DY=(C1FY(j)-C4FY(i))**2
  DZ=(C1FZ(j)-C4FZ(i))**2
  DIST=DSQRT(DX+DY+DZ)
  IF (DIST.LT.DISTM) THEN
    DISTM=DIST
    DIFMIN(i)=DISTM
    ICFMIN(i)=j
  ENDIF
  enddo
  enddo

c      write(6,*) "+cercanos"
c      do i=1,KK-1
c       write(6,*) ICIMIN(i),ICMMIN(i),ICFMIN(i)
c      enddo

IOQ=IOQ+1

ENDDO

WRITE(6,*) "TOTAL PERMEATED: ",ICOUNT

close(10)
end



