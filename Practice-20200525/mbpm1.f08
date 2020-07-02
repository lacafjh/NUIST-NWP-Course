!正压原始方程模式
!1985 10 29 by  Shen tongli
!m=20 为x方向格点数，n=16 为y方向格点数，d为网格距，rm为放大系数
!f为地转参数，w为工作数组，cla，clo分别为区域中心纬度和经度
!dt为时间步长，s为平滑系数
!ua，ub，uc分别为n-1，n，n+1时间层的x方向风速
!va，vb，vc分别为n-1，n，n+1时间层的y方向风速
!za，zb，zc分别为n-1，n，n+1时间层的位势高度
!na用于控制12小时的预报；nb用于记录时间积分步数；nt2=72用于判别
!是否积分12小时，是否该做内点平滑；nt4=6用于判定是否该做边界平滑；
!nt5用于判定是否该做时间平滑。
!zo是为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度。
program shen2

      integer                 :: m = 20, n = 16, nt2 = 72, nt4 = 6, nt5 = 36
      real, parameter         :: d = 300000.0, cla = 51.0, clo = 118.0, dt = 600.0, &
                                    zo = 2500.0, s = 0.5, c1 = dt / 2.0, c2 = dt * 2.0
      real, dimension(m:n)    :: ua, va, za, ub, vb, zb, uc, vc, zc, rm, f, w
      
      !为便于Grads做图而建立的位势高度场数据文件h.grd（包括初始场和预报场）
      open(unit=10, file='h.grd', form='unformatted', iostat=ios, status="new")
      if ( ios /= 0 ) stop "Error opening file 'h.grd'"
      !计算放大系数和地转参数，并写入数据文件中
      call cmf(rm,f,d,cla,m,n)
      open(unit=11, file='rm.dat', iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 'rm.dat'"
      write(unit=11, fmt=101, iostat=ios, advance='NO') rm
      if ( ios /= 0 ) stop "Write error in file unit 11"
      close(unit=11, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 11"
 101  format(20f10.5)
      open(unit=12, file='f.dat', iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 'f.dat'"
      write(unit=12, fmt=102, iostat=ios, advance='NO') f
      if ( ios /= 0 ) stop "Write error in file unit 12"
      close(unit=12, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 12"
 102  format(20e15.5)

      !读入初始资料场 
      open(unit=13, file='za.dat', iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file 'za.dat'"
      read(unit=13, fmt=103, iostat=istat) za
      if ( istat /= 0 ) stop "Read error in file unit 13"
      close(unit=13, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 13"
 103  format(20f6.0)

      !为便于做图,将初始数据资料写入二进制数据文件h.grd中
      write(unit=10, iostat=ios, advance='NO') ((za(i,j),i=1,m),j=1,n)     
      if ( ios /= 0 ) stop "Write error in file unit 10"

      !加入地转风子程序后此处需要修改
      !计算地转风初值   
      call cgw(ua, va, za, rm, f, d, m, n)
      !open(unit=14, file='ua.dat', iostat=ios, status="new", action="write")
      !if ( ios /= 0 ) stop "Error opening file 'ua.dat'"
      !write(unit=14, fmt=104, iostat=ios, advance='NO') ua
      !if ( ios /= 0 ) stop "Write error in file unit 14"
      !close(unit=14, iostat=ios, status="keep")
      !if ( ios /= 0 ) stop "Error closing file unit 14"
      !open(unit=15, file='va.dat', iostat=ios, status="new", action="write")
      !if ( ios /= 0 ) stop "Error opening file 'va.dat'"
      !write(unit=15, fmt=104, iostat=ios, advance='NO') va
      !if ( ios /= 0 ) stop "Write error in file unit 15"
      !close(unit=15, iostat=ios, status="keep")
      !if ( ios /= 0 ) stop "Error closing file unit 15"
      open(unit=14, file='ua.dat', iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file 'ua.dat'"
      read(unit=14, fmt=104, iostat=istat) ua
      if ( istat /= 0 ) stop "Read error in file unit 14"
      close(unit=14, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 14"
      open(unit=15, file='va.dat', iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file 'va.dat'"
      read(unit=15, fmt=104, iostat=istat) va
      if ( istat /= 0 ) stop "Read error in file unit 15"
      close(unit=15, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 15"
 104  format(20f10.5)

      !边值传送子程序
      call tbv(ub, vb, zb, ua, va, za, m, n)
      call tbv(uc, vc, zc, ua, va, za, m, n)

      !开始预报  
      do 10 na=1,2
         nb=0
         !欧拉后差积分1小时
         do 20 nn=1,6
            call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,dt,zo,m,n)
            call ti(ua,va,za,ub,vb,zb,ua,va,za,rm,f,d,dt,zo,m,n)
            nb=nb+1
  20     continue

         !边界平滑子程序
         call ssbp(za,w,s,m,n)
         call ssbp(ua,w,s,m,n)
         call ssbp(va,w,s,m,n)

         !前差积分半步
         call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,c1,zo,m,n)
         !中央差积分半步
         call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
         nb=nb+1
         !数组传送子程序
         call ta(ub,vb,zb,uc,vc,zc,m,n)
         !中央差积分一步,共积分11小时
         do 30 nn=1,66
         call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,c2,zo,m,n)
         nb=nb+1
         !打印积分步数,na大循环步,nb小循环步
         call pv(na,nb)
         !判断是否积分12小时
         if(nb == nt2) go to 80
         !判断是否做边界平滑
         if(nb/nt4*nt4 == nb) go to 40
         go to 50
  40     call ssbp(zc,w,s,m,n)
         call ssbp(uc,w,s,m,n)
         call ssbp(vc,w,s,m,n)

         !判断是否做时间平滑
  50     if(nb == nt5) go to 60
         if(nb == nt5+1) go to 60
         go to 70
         !时间平滑子程序
  60     call ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,m,n)

         !数组传送,为下一轮积分做准备
  70     call ta(ua,va,za,ub,vb,zb,m,n)
         call ta(ub,vb,zb,uc,vc,zc,m,n)
  30     continue

         !区域内点平滑
  80     call ssip(zc,w,s,m,n,2)
         call ssip(uc,w,s,m,n,2)
         call ssip(vc,w,s,m,n,2)

         !打印积分步数
         call pv(na,nb)

         !数组传送,为后12小时的积分做准备
  10     call ta(ua,va,za,uc,vc,zc,m,n)

      !存放预报结果
      open(6,file='zc.dat',status='new')
      write(6,103) zc
      close(6)
	   write(10) ((zc(i,j),i=1,m),j=1,n)    
	   open(7,file='uc.dat',status='new')
      write(7,104) uc
      close(7)
      open(8,file='vc.dat',status='new')
      write(8,104) vc
      close(8)
      stop

end program shen2

!computing map factors and coriolis parameter
!rk为圆锥常数,rlq为兰勃特投影映像平面上赤道到北极点的距离,a为地球半径
!sita为标准余纬,psx为区域中心余纬,r为模式中心到北极的距离
module LaCAfJH_NWP
contains
   subroutine cmf(rm, f, d, cla, m, n)
      implicit none
      integer, intent(inout)  :: m, n
      real, intent(inout)     :: d, cla
      real, dimension(m:n), intent(inout) :: rm, f
      integer  :: i, j, xi, yj, xi0, yj0
      real, parameter   :: rk = 0.7156, rlq = 11423370.0, a = 6371000.0,            &
               conv = 57.29578, w1 = 2.0/rk, sita = 30.0/conv, psx = (90.0-cla)/conv
      real              :: cel0, cel, r, rl, w2, sinl

      !计算模式中心到北极的距离 r
      cel0  = a*sin(sita)/rk
      cel   = (tan(psx/2.0))/(tan(sita/2.0))
      r     = cel0*cel**rk

      !确定网格坐标原点在地图坐标系中的位置
      xi0 = -(m-1)/2.0
      yj0 = r/d+(n-1)/2.0

      !求各格点至北极点的距离 rl, (xj, yi)为模式各格点在地图坐标系中的位置
      do i=1,m
         do j=1,n
            xi = xi0 + (i - 1)
            yj = yj0 - (j - 1)
            rl = sqrt(xi**2 + yj**2)*d

            !求放大系数rm和柯氏参数f
            w2       =(rl/rlq)**w1
            sinl     = (1.0 - w2)/(1.0 + w2)
            rm(i, j) = rk*rl/(a*sqrt(1.0 - sinl**2))
            f(i, j)  = 1.4584e-4*sinl
         end do
      end do
      return
   end subroutine cmf

   !computing geostrophic winds
   !请同学编写地转风初值的子程序！！！应用书中（4.134）式，要求:能用中央差的地方一定要用中央差!
   subroutine cgw(ua, va, za, rm, f, d, m, n)
      real ,dimension(m:n) :: ua, va, za, rm, f
      do i = 2, m-1
         do j = 2, n-1
            ua(i, j) = -rm(i, j)*9.8*(za(i, j + 1) - za(i, j - 1))/2/d/f(i, j)
            va(i, j) = rm(i, j)*9.8*(za(i + 1, j) - za(i - 1, j))/2/d/f(i, j)
         end do
      end do
      return
   end subroutine cgw

   !time integrations
   subroutine ti(ua, va, za, ub, vb, zb, uc, vc, zc, rm, f, d, dt, zo, m, n)
      real, dimension(m:n) :: ua, va, za, ub, vb, zb, uc, vc, zc, rm, f
      c=0.25/d
      m1=m-1
      n1=n-1
      do i=2,m1
         do j=2,n1
            e=-c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(ub(i+1,j)-ub(i,j))    &
               +(ub(i,j)+ub(i-1,j))*(ub(i,j)-ub(i-1,j))              &
               +(vb(I,j-1)+vb(i,j))*(ub(i,j)-ub(i,j-1))              &
               +(vb(I,j)+vb(i,j+1))*(ub(i,j+1)-ub(i,j))              &
               +19.6*(zb(i+1,j)-zb(i-1,j)))+f(i,j)*vb(i,j)
            uc(i,j)=ua(i,j)+e*dt
            g=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(vb(i+1,j)-vb(i,j))    &
                     +(ub(I,j)+ub(i-1,j))*(vb(i,j)-vb(i-1,j))        &
                     +(vb(I,j-1)+vb(i,j))*(vb(i,j)-vb(i,j-1))        &
                     +(vb(I,j)+vb(i,j+1))*(vb(i,j+1)-vb(i,j))        &
                     +19.6*(zb(i,j+1)-zb(i,j-1)))-f(i,j)*ub(i,j)
            vc(i,j)=va(i,j)+g*dt
         end do
      end do
      do i=2,m1
         do j=2,n1
            h=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(zb(i+1,j)-zb(i,j))    &
                  +(ub(I,j)+ub(i-1,j))*(zb(i,j)-zb(i-1,j))           &
                  +(vb(I,j-1)+vb(i,j))*(zb(i,j)-zb(i,j-1))           &
                  +(vb(I,j)+vb(i,j+1))*(zb(i,j+1)-zb(i,j))           &
               +2.0*(zb(i,j)-zo)*(ub(i+1,j)-ub(i-1,j)+vb(i,j+1)-vb(i,j-1)))
            zc(i,j)=za(i,j)+h*dt
         end do
      end do  
      return
   end subroutine ti

   !time smoothimg
   subroutine ts(ua, ub, uc, va, vb, vc, za, zb, zc, s, m, n)
      real, dimension(m:n) ::ua, va, za, ub, vb, zb, uc, vc, zc
       m1=m-1
         n1=n-1
         do i=2,m1
            do j=2,n1
               ub(i,j)=ub(i,j)+s*(ua(i,j)+uc(i,j)-2.0*ub(i,j))/2.0
               vb(i,j)=vb(i,j)+s*(va(i,j)+vc(i,j)-2.0*vb(i,j))/2.0
               zb(i,j)=zb(i,j)+s*(za(i,j)+zc(i,j)-2.0*zb(i,j))/2.0
            end do
         end do
      return
   end subroutine ts

   !space smoothing for internal points 区域内5点平滑(正逆平滑)
   !请同学编写区域内5点平滑(正逆平滑)的子程序！！！应用书中（4.126）式
   !注：此程序必须设计成开关形式，保证既可选做正逆平滑，又可选做正平滑   l=1为只执行正平滑，l=2为执行正逆平滑.
   subroutine ssip(a, w, s, m,n, l)
      real, dimension(m:n) :: a, w
      !计算五点平滑（正平滑）
      if (l == 1) then
         do i = 2, m - 1
            do j = 2, n - 1
               w(i,j) = a(i,j) + s*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1) - 4*a(i,j))/4.0
            end do
         end do
         do i = 2, m - 1
            do j = 2, n - 1
               a(i, j) = w(i, j)
            end do
         end do
      return
      !计算五点平滑（正逆平滑）
      else
         do i = 2, m - 1
            do j = 2, n - 1
               w(i,j) = a(i,j) + s*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1) - 4*a(i,j))/4.0
            end do
         end do
         do i = 2, m - 1
            do j = 2, n - 1
               a(i, j) = w(i, j)
            end do
         end do
         do i = 2,m - 1
            do j = 2, n - 1
               w(i,j) = a(i,j) - s*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1) - 4*a(i,j))/4.0
            end do
         end do
         do i = 2,m - 1
            do j = 2, n - 1
               a(i, j) = w(i, j)
            end do
         end do
      end if
      return
   end subroutine ssip

   !space smoothing for boundary points 边界九点平滑
   subroutine ssbp(a,w,s,m,n)
      real, dimension(m:n) :: a, w
      m1=m-1
      m3=m-3
      n1=n-1
      n2=n-2
      n3=n-3
      do i=2,m1
         do j=2,n1,n3
            w(i,j)=a(i,j)+0.5*s*(1.0-s)*(a(i-1,j)                       &
               +a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))+0.25*s*s*        &
               (a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
         end do
      end do
      do i=2,m1,m3
         do j=3,n2
            w(i,j)=a(i,j)+0.5*s*(1.0-s)*(a(i-1,j)                       &
               +a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))+0.25*s*s*        &
               (a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
         end do
      end do
      do i=2,m1
         do j=2,n1,n3
            a(i,j)=w(i,j)
         end do
      end do
      do i=2,m1,m3
         do j=3,n2
            a(i,j)=w(i,j)
         end do
      end do
      return
   end subroutine ssbp

   !transmiting arrays  数组传送
   subroutine ta(ua, va, za, ub, vb, zb, m, n)
      real, dimension(m:n) :: ua, va, za, ub, vb, zb
      do i=1,m
         do j=1,n
            ua(i,j)=ub(i,j)
            va(i,j)=vb(i,j)
            za(i,j)=zb(i,j)
         end do
      end do
   return
   end subroutine ta

   !transmiting boundary valaus  赋固定边界值
   subroutine tbv(ua, va, za, ub, vb, zb, m, n)
      real, dimension(m:n) :: ua, va, za, u, vb, zb
      m1=m-1
      n1=n-1
      do  i=1,m
         do j=1,n,n1
            ua(i,j)=ub(i,j)
            va(i,j)=vb(i,j)
            za(i,j)=zb(i,j)
         end do
      end do
      do i=1,m,m1
         do j=1,n
            ua(i,j)=ub(i,j)
            va(i,j)=vb(i,j)
            za(i,j)=zb(i,j)
         end do
      end do
      return
   end subroutine tbv

   !printing variables  打印积分步数
   subroutine pv(na,nb)
      write(*,100)na,nb
100   format(//////5x,3hna=,i3,5x,3hnb=,i2/)
   return
   end subroutine pv
end module LaCAfJH_NWP
