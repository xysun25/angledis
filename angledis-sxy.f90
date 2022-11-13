module atomdef   ! 封装atomdef变量
    implicit none
    type atom    ! 自定义一个叫atom的数据类型,后面可以用来声明变量的数据类型为atom
        character(len=10)::element   ! 记录元素名称
        real::x,y,z,m,e,sig,eps   ! 记录原子坐标与 L-J 势中截断距离和能量
        end type   ! 自定义数据类型结束
    
    type mol   ! 自定义一个叫mol 的数据类型
        character(len=10)::numberofmol   ! 记录分子数目
        real::x,y,z   ! 记录分子坐标
        end type   ! 自定义数据类型结束
    
end module 

! 主程序
program main
    use atomdef   ! 调用atomdef模块
    implicit none
    character(len=20)::filename,mmol1,mmol2,mmol3,trjfile   ! 记录文件名称
    integer::num, numberofAIM
    integer i,j,k,u,l,p,q,r,trj,ls,u1,u2,u3,ff1,ff2,lnn,np,nn,mm,up1,upan,uu1,uu2,uu3,nsite,tt,qt
    integer numberoftrj,o1,h1,oo,ohh,ha
    integer num1,num2,num3,numberofAIM1,numberofAIM2,numberofAIM3
    integer numberofatom1,numberofatom2,numberofatom3,numberofatoms_2
    type(atom),allocatable::atoms_1(:),atoms_2(:),atoms(:,:)   ! 声明两个atom类型的变量，是两个可变一维数组，和一个atom类型的变量：可变二维数组
    real X,Y,Z,M,Xcom,Ycom,Zcom,x00,y00,dr,high,oh,ohx,ohy,ohz,ohl,oohx,oohy,oohz,oohhx,oohhy,oohhz,oohhl,ohhx,ohhy,ohhz,ool
    real a0,b0,c0,a,b,c,ax,by,cz,interv,theta,bin,angle,c1,c2,thetah,an,bn,cn,anx,any,anz,bnx,bny,bnz,cnx,cny,cnz
    integer numberofinterv
    integer::numberofions,numberofatoms
    real::o(2),h(2)
    real,allocatable::Hang(:),angdis(:,:,:)   ! 声明一个可变一维数组，一个可变三维数组
    
    ! 要写in文件作为输入，里面包含下面声明和操作文件的内容
    read (*,*) filename   ! in文件，包含下面文件名称和变量名称
    write(*,*) 'the input file is:', filename
    open(100,file=filename,status="old")
    read(100,*) trjfile,numberoftrj
    write(*,*) 'trjfile and numberof trj',trjfile,numberoftrj
    read (100,*) mmol1,mmol2,mmol3
    write(*,*) 'the mmol files are:', mmol1,mmol2,mmol3
    read (100,*) num1,num2,num3
    write(*,*) 'the number of molecules are:', num1,num2,num3
    read(100,*) a0,b0,c0,ax,by,cz
    read(100,*) interv, numberofinterv, bin
    read(100,*) u1, u2, u3
    read(100,*) uu1, uu2, uu3
    
    o=[4,12]
    h=[16,18]
    o1=19
    h1=20
    
    ! 处理EmimOH.mmol 文件
    open(11,file=mmol1,status="old")   ! 打开已存在的mmol1文件，指定文件代码为11
    read(11,*)   ! 读这个文件
    read(11,*) numberofAIM1    ! 20
    allocate(atoms_1(numberofAIM1))   ! allocate 配置atoms_1一维数组的内存空间，后面要释放的
    ! 循环
    do i=1,numberofAIM1   ! i起始数值，numberofAIM1终止数值，增量默认为1;循环20次
        ! 使用atom类型变量atoms_1中element等元素的方法：变量和元素用%间隔
read(11,*)atoms_1(i)%element,atoms_1(i)%x,atoms_1(i)%y,atoms_1(i)%z,atoms_1(i)%m,atoms_1(i)%e,atoms_1(i)%sig,atoms_1(i)%eps        
    end do
    rewind(11)
    
    ! 处理NTf2.mmol 文件
    open(12,file=mmol2,status="old")
    read(12,*)
    read(12,*) numberofAIM2   ! 14
    allocate (atoms_2(numberofAIM2))
    do j=1,numberofAIM2
read(12,*)atoms_2(j)%element,atoms_2(j)%x,atoms_2(j)%y,atoms_2(j)%z,atoms_2(j)%m,atoms_2(j)%e,atoms_2(j)%sig,atoms_2(j)%eps
    end do
    rewind(12)
    
    ! 处理mxene.mmol 文件
    open(13,file=mmol3,status="old")
    read(13,*) numberofAIM3   ! 14
    close(13)
    
    numberofatom3=num3*numberofAIM3   ! 600*14=8400
    num=num1*numberofAIM1+num2*numberofAIM2+num3*numberofAIM3   ! 500*20+500*15+600*14=25900
    numberofAIM=numberofAIM1+numberofAIM2   ! 20+15=35
    nsite=int(180/bin)   ! 180/0.5=360
    mm=num1*num3   ! 500*600=300000
    ha=2/numberofinterv   !2/100(140)(70)
     
    ! 打开轨迹文件
    open(10,file=trjfile,status="old")
    allocate(angdis(4,numberofinterv,nsite))   ! 配置三维可变数组angdis的内存angdis(4,100,360)
    allocate(Hang(nsite))   ! 配置一维数组Hang的内存Hang(360)
    
    open(16,file="cation_angl.dat",status="unknown")
    open(17,file="cation_ang2.dat",status="unknown")
    open(18,file="anion_angl.dat",status="unknown")
    open(19,file="anion_ang2.dat",status="unknown")
    open(20,file="hydrogen_ang.dat",status="unknown")
    
    ff1=0
    ff2=0
    upan=0
    up1=0
    np=100
    nn=0
    ! 在轨迹文件里循环
    do trj=1,numberoftrj   ! trj起始数值，numberoftrj终止数值（帧个数）
        allocate(atoms(num,numberofAIM))   ! atoms是自定义数据类型atom的一个二维数组atoms(25900,35)
        read(10,*)
        read(10,*)
        do k=1,num1   ! 循环500次
            do u=1,numberofAIM1   ! 循环20次
                read(10,*) atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
            end do
        end do 
        do k=num1+1,num2+num1
		    do u=1,numberofAIM2
			    read (10,*)atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
		    end do
	    end do

	    do k=num2+num1+1,num2+num1+num3
		    do r=1,numberofAIM3
			    read(10,*) atoms(k,r)%element,atoms(k,r)%x,atoms(k,r)%y,atoms(k,r)%z
		    end do
        end do

        ! 调整穿过边界的阳离子
        ! 转移阳离子
        do k=1,num1    ! num1=500 循环500次 一共有500个阳离子分子
            tt=0
            ! 求得x边界
            x00=atoms(k,1)%x   ! x00等于每个阳离子中的第一个原子的x坐标，atoms是一个自定义数据类型atom的变量，二维数组，使用变量atoms中x这个元素
            do u=2,numberofAIM1   ! numberofAIM1=20 循环20次 一个分子有20个原子
                dr=abs(atoms(k,u)%x-x00)   ! 阳离子中其他19个原子x坐标减去第一个原子的x坐标
                if (dr>(0-a0)) then   ! 如果一个阳离子中的原子穿越了边界，到了另一边
                    tt=1
                end if 
            end do 
            if (tt==1) then
                do u=1,numberofAIM1   ! 20，循环20次，对一个阳离子的20个原子操作
                    if (atoms(k,u)%x>0.0) then 
                        atoms(k,u)%x=atoms(k,u)%x+a0*2
                    end if 
                end do 
            end if
            
            ! 求得y边界
            y00 = atoms(k,1)%y
            tt = 0
            do u=2,numberofAIM1
                dr = abs(atoms(k,u)%y-y00)
                if(dr >(0-b0)) then
                    tt =1 
                end if
            end do
            if (tt == 1) then
                do u=1,numberofAIM1
                    if(atoms(k,u)%y >0.0) then
                        atoms(k,u)%y = atoms(k,u)%y + b0*2
                    end if
                end do
            end if
        end do 


        ! 调整穿过边界的阳离子
        ! 调整阴离子
        do k=num1+1,num2+num1
            tt = 0
            ! 求得x边界
            x00 = atoms(k,1)%x
            do u=2,numberofAIM2
                dr = abs(atoms(k,u)%x-x00)
                if(dr >(0-a0)) then
                    tt =1    ! 如果一个阳离子中的原子穿越了边界，到了另一边
                end if
            end do
            if (tt ==1) then
                do u=1,numberofAIM2
                    if(atoms(k,u)%x > 0.0) then
                        atoms(k,u)%x = atoms(k,u)%x + a0*2
                    end if
                end do
            end if
            ! 求得x边界
            y00 = atoms(k,1)%y
            tt = 0
            do u=2,numberofAIM2
                dr = abs(atoms(k,u)%y-y00)
                if(dr >(0-b0)) then
                    tt =1     ! 如果一个阳离子中的原子穿越了边界，到了另一边
                end if
            end do
            if (tt == 1) then
                do u=1,numberofAIM2
                    if(atoms(k,u)%y >0.0) then
                        atoms(k,u)%y = atoms(k,u)%y + b0*2
                    end if
                 end do
            end if
    end do
    
    ! write(*,*)"coordinate of atoms transfers have completed" 
    
    
    
    ! 求阳离子离子质心
    do k=1,num1   ! 500个阳离子，循环500次
        X=0.0
        Y=0.0
        Z=0.0
        M=0.0   ! 均设置初始值为0.0，M是质量和
        do u=1,numberofAIM1   ! 20,阳离子中有20个原子，循环20次 ,找到阳离子的质量中心
            X=X+atoms(k,u)%x*atoms_1(u)%m
            Y=Y+atoms(k,u)%y*atoms_1(u)%m
		    Z=Z+atoms(k,u)%z*atoms_1(u)%m
		    M=M+atoms_1(u)%m
        end do 
        Xcom=X/M
	    Ycom=Y/M
	    Zcom=Z/M
        do ls=1,numberofinterv   ! numberofinterv(层间间隔层的数目)=70/100/140,循环70/10/140次
            b=by+(ls-1)*interv   ! interv=0.1,在1.0层间距的体系，by=-5，b=-5+（2-1）*0.1=-4.9，如ls=2
            if(Ycom>=b.and.Ycom<(b+interv))then   ! 层间距1.0的体系，b=4.9,b+interv=5.0
                if (Zcom<cz) then   ! 层间距1.0体系，zc=-52
                    ff1=ff1+1
                    ! N 原子的（x,y,z） 阳离子中第一个原子 u1=1
                    anx=atoms(k,u1)%x    ! N原子的x坐标
                    any=atoms(k,u1)%y    ! N原子的y坐标
                    anz=atoms(k,u1)%z    ! N原子的z坐标
                    ! C 原子的（x,y,z） 阳离子中第二个原子 u2=2
                    bnx=atoms(k,u2)%x
                    bnx=atoms(k,u2)%y
                    bnz=atoms(k,u2)%z
                    ! 第二个N 原子的（x,y,z）,也就是阳离子中的第三个原子 u3=3
                    cnx=atoms(k,u3)%x
                    cnx=atoms(k,u3)%y
                    cnz=atoms(k,u3)%z
                    ! 三点确定一个平面（环），求阳离子平面（环）的法向量
                    an = (bny-any)*(cnz-anz)-(cny-any)*(bnz-anz)    ! 法向量
				    bn = (bnz-anz)*(cnx-anx)-(cnz-anz)*(bnx-anx)
				    cn = (bnx-anx)*(cny-any)-(cnx-anx)*(bny-any)
				
				    theta = bn/sqrt((an**2 + bn**2 + cn**2))
				    theta = ACOS(theta)
				    theta = theta*180.0/3.1415926
                    write(16,"(I5,2x,F12.6)") ls,theta   ! 写入cation_ang1.dat，阳离子在负极
                elseif(Zcom>(0-cz))then
				    ff2=ff2+1
                    anx = atoms(k,u1)%x  ! A (x,y,z)
				    any = atoms(k,u1)%y
				    anz = atoms(k,u1)%z
				    bnx = atoms(k,u2)%x  ! B (x,y,z)
				    bny = atoms(k,u2)%y
				    bnz = atoms(k,u2)%z
				    cnx = atoms(k,u3)%x  ! C (x,y,z)
				    cny = atoms(k,u3)%y
				    cnz = atoms(k,u3)%z
   
    				an = (bny-any)*(cnz-anz)-(cny-any)*(bnz-anz) ! normal vector 
	    			bn = (bnz-anz)*(cnx-anx)-(cnz-anz)*(bnx-anx)
		    		cn = (bnx-anx)*(cny-any)-(cnx-anx)*(bny-any)
				
			    	theta = bn/sqrt((an**2 + bn**2 + cn**2))
				    theta = ACOS(theta)
				    theta = theta*180.0/3.1415926
			    	write(17,"(I5,2x,F12.6)") ls,theta   ! 写入cation_ang2.dat阳离子在正极
			    end if
		    end if
        end do
        if(abs(Ycom)<=7.and.Zcom<=cz)then
            ohx = atoms(k,o1)%x-atoms(k,h1)%x
            ohy = atoms(k,o1)%y-atoms(k,h1)%y
            ohz = atoms(k,o1)%z-atoms(k,h1)%z
            oh = sqrt((ohx**2 + ohy**2 + ohz**2))
            do j=num2+num1+1,num2+num1+num3/2
                do i=1, 2
                    oo = o(i)
                    ohh = h(i)
                    oohx = atoms(j,oo)%x-atoms(k,h1)%x
                    oohy = atoms(j,oo)%y-atoms(k,h1)%y
                    oohz = atoms(j,oo)%z-atoms(k,h1)%z
                    ohhx = atoms(j,ohh)%x-atoms(k,o1)%x
                    ohhy = atoms(j,ohh)%y-atoms(k,o1)%y
                    ohhz = atoms(j,ohh)%z-atoms(k,o1)%z
                    oohhx = atoms(j,oo)%x-atoms(j,ohh)%x
                    oohhy = atoms(j,oo)%y-atoms(j,ohh)%y
                    oohhz = atoms(j,oo)%z-atoms(j,ohh)%z
                    ool = sqrt((oohx**2 + oohy**2 + oohz**2))
                    ohl = sqrt((ohhx**2 + ohhy**2 + ohhz**2))
                    oohhl = sqrt((oohhx**2 + oohhy**2 + oohhz**2))
                    an = (ohx*oohx+ohy*oohy+ohz*oohz)/oh/ool
                    bn = (ohhx*oohhx+ohhy*oohhy+ohhz*oohhz)/ohl/oohhl
                    if(ool<4)then
                        nn = nn+1
                        thetah = (ohx*oohx+ohy*oohy+ohz*oohz)/oh/ool
                        thetah = ACOS(thetah)
                        thetah = thetah*180.0/3.1415926
                        write(20,"(2F12.6)") ool,thetah
                    elseif(ohl<4)then
                        nn = nn+1
                        thetah = (ohhx*oohhx+ohhy*oohhy+ohhz*oohhz)/ohl/oohhl
                        thetah = ACOS(thetah)
                        thetah = thetah*180.0/3.1415926
                        write(20,"(2F12.6)") ohl,thetah
                    end if
                end do
            end do
        end if
    end do 

    ! 求阴离子质心
    do k=num1+1,num2+num1
	    X=0.0
        Y=0.0
        Z=0.0
        M=0.0
            do u=1,numberofAIM2
		        X=X+atoms(k,u)%x*atoms_2(u)%m
		        Y=Y+atoms(k,u)%y*atoms_2(u)%m
		        Z=Z+atoms(k,u)%z*atoms_2(u)%m
		        M=M+atoms_2(u)%m
	        end do 
	        Xcom=X/M
	        Ycom=Y/M
	        Zcom=Z/M
	        do ls=1,numberofinterv
		        b=by+(ls-1)*interv
		        if(Ycom>=b.and.Ycom<(b+interv))then
			        if(Zcom<cz)then
				        upan=upan+1
                        anx = atoms(k,uu1)%x  ! A (x,y,z)
                        any = atoms(k,uu1)%y
                        anz = atoms(k,uu1)%z
                        bnx = atoms(k,uu2)%x  ! B (x,y,z)
                        bny = atoms(k,uu2)%y
                        bnz = atoms(k,uu2)%z

                        an = anx-bnx
                        bn = any-bny
                        cn = any-bny
				        theta = bn/sqrt((an**2 + bn**2 + cn**2))
			            theta = ACOS(theta)
				        theta = theta*180.0/3.1415926
				        write(18,"(I5,2x,F12.6)") ls,theta
                    elseif(Zcom>(0-cz))then
				        up1=up1+1
				        anx = atoms(k,uu1)%x  ! A (x,y,z)
                        any = atoms(k,uu1)%y
                        anz = atoms(k,uu1)%z
                        bnx = atoms(k,uu2)%x  ! B (x,y,z)
                        bny = atoms(k,uu2)%y
                        bnz = atoms(k,uu2)%z

                        an = anx-bnx
                        bn = any-bny
                        cn = any-bny
				        theta = bn/sqrt((an**2 + bn**2 + cn**2))
				        theta = ACOS(theta)
				        theta = theta*180.0/3.1415926
                        write(19,"(I5,2x,F12.6)") ls,theta
			        end if
		        end if
	        end do
        end do

        deallocate(atoms)
        end do  
        rewind(16)
        rewind(17)
        rewind(18)
        rewind(19)
        rewind(20) 
        
        ! 开始统计
        do k=1, ff1   
            read(16,*)ls, theta
            dr=abs(theta-90)/bin
            u=floor(dr)
            angdis(1,ls,u)=angdis(1,ls,u)+1
        end do
        do k=1,ff2
            read(17,*)ls, theta
            dr=abs(theta-90)/bin
            u=floor(dr)
            angdis(2,ls,u)=angdis(2,ls,u)+1
        end do
        do k=1,upan
            read(18,*) ls,theta
            dr=abs(theta-90)/bin
            u=floor(dr)
            angdis(3,ls,u)=angdis(3,ls,u)+1
        end do
        do k=1,up1
            read(19,*)ls, theta
            dr=abs(theta-90)/bin
            u=floor(dr)
            angdis(4,ls,u)=angdis(4,ls,u)+1
        end do
        !do k=1,nn
         !   read(20,*)ls, theta
          !  dr=theta/bin
           ! u=floor(dr)
           ! Hang(u)=Hang(u)+1
        !end do


        !average
        open(15,file="averageangle",status="unknown")
        open(14,file="Hangle",status="unknown")
        do ls=1,numberofinterv
	        do k=1,nsite/2
write(15,"(F8.3,2x,F8.3,2x,4F12.8)")(ls*interv-7),(90-k*bin),angdis(1,ls,k)/numberoftrj,angdis(2,ls,k)/numberoftrj,&
    angdis(3,ls,k)/numberoftrj,angdis(4,ls,k)/numberoftrj
	        end do
	        !write(15,"(I4,2x,4F10.4)")ls,angledis(1,ls),angledis(2,ls),angledis(3,ls),angledis(4,ls)
         end do
      !  do i=1,nsite/2
!write(14,"(F8.3,2x,5F12.8)")i*bin,Hang(i)/numberoftrj,sum(angdis(1,:,i))/numberoftrj,sum(angdis(2,:,i))/numberoftrj,&
    !sum(angdis(3,:,i))/numberoftrj,sum(angdis(4,:,i))/numberoftrj
       ! end do
       ! do i=1+nsite/2,nsite
!write(14,"(F8.3,2x,F12.8)")i*bin,Hang(i)/numberoftrj
 !       end do
deallocate(angdis)
!deallocate(Hang)


stop
end program
        