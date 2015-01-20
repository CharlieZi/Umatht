      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),
     1 COORD(*)

      
      TEMP_0 = STATEV(2)
      T_K = TEMP_0+ 273.0


      IF ((KSTEP.EQ.1).AND.(KINC.EQ.1)) THEN
            STATEV(1) = 0.0001
      ENDIF

      ALPHA = STATEV(1)


      DALPHA = 67198.43*EXP(-67342.89/8.314/(T_K))*(ALPHA**0.2857)*(1.0-ALPHA)**1.2373*DTIME
      ALPHA = ALPHA + DALPHA



      STATEV(1) = ALPHA
      FIELD(1) = STATEV(1) 



      RETURN
      END





      subroutine uexpan(expan,dexpandt,temp,time,dtime,predef,dpred,
     $     statev,cmname,nstatv,noel)
c
      include 'aba_param.inc'
c
      character*80 cmname
c
      dimension expan(*),dexpandt(*),temp(2),time(2),predef(*),
     $     dpred(*),statev(nstatv)
c

      alpha_chemi = 0
      IF (STATEV(1)  .ge. 0.5) THEN
      	alpha_chemi = 3.9d-05
      ENDIF


      STATEV(2) = temp(1)

      STATEV(3) = alpha
      STATEV(4) = temp(2)


      alpha_1 = 0
      alpha_2 = -2.63d-05
      alpha_3 = -2.63d-05 + alpha_chemi


	expan(1) = alpha_1*temp(2)
	expan(2) = alpha_2*temp(2)
	expan(3) = alpha_3*temp(2)

c
      return
      end