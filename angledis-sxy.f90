module atomdef   ! ��װatomdef����
    implicit none
    type atom    ! �Զ���һ����atom����������,�����������������������������Ϊatom
        character(len=10)::element   ! ��¼Ԫ������
        real::x,y,z,m,e,sig,eps   ! ��¼ԭ�������� L-J ���нضϾ��������
        end type   ! �Զ����������ͽ���
    
    type mol   ! �Զ���һ����mol ����������
        character(len=10)::numberofmol   ! ��¼������Ŀ
        real::x,y,z   ! ��¼��������
        end type   ! �Զ����������ͽ���
    
end module 

! ������
program main
    use atomdef   ! ����atomdefģ��
    implicit none
    character(len=20)::filename,mmol1,mmol2,mmol3,trjfile   ! ��¼�ļ�����
    integer::num, numberofAIM
    integer i,j,k,u,l,p,q,r,trj,ls,u1,u2,u3,ff1,ff2,lnn,np,nn,mm,up1,upan,uu1,uu2,uu3,nsite,tt,qt
    integer numberoftrj,o1,h1,oo,ohh,ha
    integer num1,num2,num3,numberofAIM1,numberofAIM2,numberofAIM3
    integer numberofatom1,numberofatom2,numberofatom3,numberofatoms_2
    type(atom),allocatable::atoms_1(:),atoms_2(:),atoms(:,:)   ! ��������atom���͵ı������������ɱ�һά���飬��һ��atom���͵ı������ɱ��ά����
    real X,Y,Z,M,Xcom,Ycom,Zcom,x00,y00,dr,high,oh,ohx,ohy,ohz,ohl,oohx,oohy,oohz,oohhx,oohhy,oohhz,oohhl,ohhx,ohhy,ohhz,ool
    real a0,b0,c0,a,b,c,ax,by,cz,interv,theta,bin,angle,c1,c2,thetah,an,bn,cn,anx,any,anz,bnx,bny,bnz,cnx,cny,cnz
    integer numberofinterv
    integer::numberofions,numberofatoms
    real::o(2),h(2)
    real,allocatable::Hang(:),angdis(:,:,:)   ! ����һ���ɱ�һά���飬һ���ɱ���ά����
    
    ! Ҫдin�ļ���Ϊ���룬����������������Ͳ����ļ�������
    read (*,*) filename   ! in�ļ������������ļ����ƺͱ�������
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
    
    ! ����EmimOH.mmol �ļ�
    open(11,file=mmol1,status="old")   ! ���Ѵ��ڵ�mmol1�ļ���ָ���ļ�����Ϊ11
    read(11,*)   ! ������ļ�
    read(11,*) numberofAIM1    ! 20
    allocate(atoms_1(numberofAIM1))   ! allocate ����atoms_1һά������ڴ�ռ䣬����Ҫ�ͷŵ�
    ! ѭ��
    do i=1,numberofAIM1   ! i��ʼ��ֵ��numberofAIM1��ֹ��ֵ������Ĭ��Ϊ1;ѭ��20��
        ! ʹ��atom���ͱ���atoms_1��element��Ԫ�صķ�����������Ԫ����%���
read(11,*)atoms_1(i)%element,atoms_1(i)%x,atoms_1(i)%y,atoms_1(i)%z,atoms_1(i)%m,atoms_1(i)%e,atoms_1(i)%sig,atoms_1(i)%eps        
    end do
    rewind(11)
    
    ! ����NTf2.mmol �ļ�
    open(12,file=mmol2,status="old")
    read(12,*)
    read(12,*) numberofAIM2   ! 14
    allocate (atoms_2(numberofAIM2))
    do j=1,numberofAIM2
read(12,*)atoms_2(j)%element,atoms_2(j)%x,atoms_2(j)%y,atoms_2(j)%z,atoms_2(j)%m,atoms_2(j)%e,atoms_2(j)%sig,atoms_2(j)%eps
    end do
    rewind(12)
    
    ! ����mxene.mmol �ļ�
    open(13,file=mmol3,status="old")
    read(13,*) numberofAIM3   ! 14
    close(13)
    
    numberofatom3=num3*numberofAIM3   ! 600*14=8400
    num=num1*numberofAIM1+num2*numberofAIM2+num3*numberofAIM3   ! 500*20+500*15+600*14=25900
    numberofAIM=numberofAIM1+numberofAIM2   ! 20+15=35
    nsite=int(180/bin)   ! 180/0.5=360
    mm=num1*num3   ! 500*600=300000
    ha=2/numberofinterv   !2/100(140)(70)
     
    ! �򿪹켣�ļ�
    open(10,file=trjfile,status="old")
    allocate(angdis(4,numberofinterv,nsite))   ! ������ά�ɱ�����angdis���ڴ�angdis(4,100,360)
    allocate(Hang(nsite))   ! ����һά����Hang���ڴ�Hang(360)
    
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
    ! �ڹ켣�ļ���ѭ��
    do trj=1,numberoftrj   ! trj��ʼ��ֵ��numberoftrj��ֹ��ֵ��֡������
        allocate(atoms(num,numberofAIM))   ! atoms���Զ�����������atom��һ����ά����atoms(25900,35)
        read(10,*)
        read(10,*)
        do k=1,num1   ! ѭ��500��
            do u=1,numberofAIM1   ! ѭ��20��
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

        ! ���������߽��������
        ! ת��������
        do k=1,num1    ! num1=500 ѭ��500�� һ����500�������ӷ���
            tt=0
            ! ���x�߽�
            x00=atoms(k,1)%x   ! x00����ÿ���������еĵ�һ��ԭ�ӵ�x���꣬atoms��һ���Զ�����������atom�ı�������ά���飬ʹ�ñ���atoms��x���Ԫ��
            do u=2,numberofAIM1   ! numberofAIM1=20 ѭ��20�� һ��������20��ԭ��
                dr=abs(atoms(k,u)%x-x00)   ! ������������19��ԭ��x�����ȥ��һ��ԭ�ӵ�x����
                if (dr>(0-a0)) then   ! ���һ���������е�ԭ�Ӵ�Խ�˱߽磬������һ��
                    tt=1
                end if 
            end do 
            if (tt==1) then
                do u=1,numberofAIM1   ! 20��ѭ��20�Σ���һ�������ӵ�20��ԭ�Ӳ���
                    if (atoms(k,u)%x>0.0) then 
                        atoms(k,u)%x=atoms(k,u)%x+a0*2
                    end if 
                end do 
            end if
            
            ! ���y�߽�
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


        ! ���������߽��������
        ! ����������
        do k=num1+1,num2+num1
            tt = 0
            ! ���x�߽�
            x00 = atoms(k,1)%x
            do u=2,numberofAIM2
                dr = abs(atoms(k,u)%x-x00)
                if(dr >(0-a0)) then
                    tt =1    ! ���һ���������е�ԭ�Ӵ�Խ�˱߽磬������һ��
                end if
            end do
            if (tt ==1) then
                do u=1,numberofAIM2
                    if(atoms(k,u)%x > 0.0) then
                        atoms(k,u)%x = atoms(k,u)%x + a0*2
                    end if
                end do
            end if
            ! ���x�߽�
            y00 = atoms(k,1)%y
            tt = 0
            do u=2,numberofAIM2
                dr = abs(atoms(k,u)%y-y00)
                if(dr >(0-b0)) then
                    tt =1     ! ���һ���������е�ԭ�Ӵ�Խ�˱߽磬������һ��
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
    
    
    
    ! ����������������
    do k=1,num1   ! 500�������ӣ�ѭ��500��
        X=0.0
        Y=0.0
        Z=0.0
        M=0.0   ! �����ó�ʼֵΪ0.0��M��������
        do u=1,numberofAIM1   ! 20,����������20��ԭ�ӣ�ѭ��20�� ,�ҵ������ӵ���������
            X=X+atoms(k,u)%x*atoms_1(u)%m
            Y=Y+atoms(k,u)%y*atoms_1(u)%m
		    Z=Z+atoms(k,u)%z*atoms_1(u)%m
		    M=M+atoms_1(u)%m
        end do 
        Xcom=X/M
	    Ycom=Y/M
	    Zcom=Z/M
        do ls=1,numberofinterv   ! numberofinterv(����������Ŀ)=70/100/140,ѭ��70/10/140��
            b=by+(ls-1)*interv   ! interv=0.1,��1.0�������ϵ��by=-5��b=-5+��2-1��*0.1=-4.9����ls=2
            if(Ycom>=b.and.Ycom<(b+interv))then   ! ����1.0����ϵ��b=4.9,b+interv=5.0
                if (Zcom<cz) then   ! ����1.0��ϵ��zc=-52
                    ff1=ff1+1
                    ! N ԭ�ӵģ�x,y,z�� �������е�һ��ԭ�� u1=1
                    anx=atoms(k,u1)%x    ! Nԭ�ӵ�x����
                    any=atoms(k,u1)%y    ! Nԭ�ӵ�y����
                    anz=atoms(k,u1)%z    ! Nԭ�ӵ�z����
                    ! C ԭ�ӵģ�x,y,z�� �������еڶ���ԭ�� u2=2
                    bnx=atoms(k,u2)%x
                    bnx=atoms(k,u2)%y
                    bnz=atoms(k,u2)%z
                    ! �ڶ���N ԭ�ӵģ�x,y,z��,Ҳ�����������еĵ�����ԭ�� u3=3
                    cnx=atoms(k,u3)%x
                    cnx=atoms(k,u3)%y
                    cnz=atoms(k,u3)%z
                    ! ����ȷ��һ��ƽ�棨��������������ƽ�棨�����ķ�����
                    an = (bny-any)*(cnz-anz)-(cny-any)*(bnz-anz)    ! ������
				    bn = (bnz-anz)*(cnx-anx)-(cnz-anz)*(bnx-anx)
				    cn = (bnx-anx)*(cny-any)-(cnx-anx)*(bny-any)
				
				    theta = bn/sqrt((an**2 + bn**2 + cn**2))
				    theta = ACOS(theta)
				    theta = theta*180.0/3.1415926
                    write(16,"(I5,2x,F12.6)") ls,theta   ! д��cation_ang1.dat���������ڸ���
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
			    	write(17,"(I5,2x,F12.6)") ls,theta   ! д��cation_ang2.dat������������
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

    ! ������������
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
        
        ! ��ʼͳ��
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
        