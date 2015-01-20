c
c
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,
     $     dtemp,dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,
     $     props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)
c
      include 'aba_param.inc'
c
      character*80 cmname
c
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     $     dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),time(2),
     $     predef(1),dpred(1),props(nprops),coords(3)
c


C unitSwitch, must Check before run

      unitSwitch = 'mm'


C percentage of fiber

      v_f = 0.55


C  国际单位制 m\kg

      rho = 1180
      H_r = 191000
      kx = 2.0
      kyz =2.0




      IF ((kstep.EQ.1).AND.(kinc.eq.1)) THEN
            statev(1) = 0.0001
      ENDIF



C 热力学计算
      T_k = temp+ 273.0
      alpha = statev(1)
      specht=(-0.0078+0.00112*T_K-0.00274*alpha+1.2e-6*T_K*T_K-0.0516*alpha*alpha)*1000.0		

C must be after specht, and before dq
      IF (unitSwitch .EQ. 'mm') THEN
            rho = rho*1e12
            H_r = H_r/1e6
            kx = kx
            kyz = kyz
            specht = specht/1e6
      ENDIF


      dudt = specht
      du = dudt*dtemp


C 固化度系数计算
      dalpha = 67198.43*exp(-67342.89/8.314/(T_k))*(alpha**0.2857)*(1.0-alpha)**1.2373*dtime
      alpha = alpha + dalpha

      dq = (1-v_f)*H_r*dalpha/dtime
      q = dq*dtime

      u = u+du - q


c      kx=-1.875+0.00955*T_K-0.232*alpha-5.672e-6*T_K*T_K+0.00725*alpha*alpha
c      kyz=-1.337+0.00654*T_K-0.266*alpha-2.236e-6*T_K*T_K+0.00777*alpha*alpha



C  热流  input flux = -[k]*{dtemdx}
      do i=1, ntgrd
         flux(i) = -kx*dtemdx(i) 
      end do

C 热传导矩阵
      do i=1, ntgrd
         dfdg(i,i) = kx
      end do
      dfdg(2,3) = kyz
c





C output Flags

      statev(1) = alpha
      statev(2) = dalpha
      statev(3) = dtime
      statev(4) = specht
      statev(5) = q
      statev(6) = kx
      statev(7) = kyz




      return
      end
      
