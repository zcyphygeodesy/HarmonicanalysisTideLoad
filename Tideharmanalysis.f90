      subroutine Tideharmanalysis(dtmfl,tidefl,knd,dk,itd)
      implicit none
	character*800::dtmfl,tidefl,infl,outfl,chfl,prfl
	character*800::line,str(80),nstr
	integer i,j,nlon,nlat,sn,kln,kk,nk,astat(12),maxn,nd,m,n
	integer knd,nlon0,nlat0,pm
	real*8 dk,itd,rec(800),hd(6),hd0(6),tmp,doodson
	real*8::GRS(6),ae,rst1(4),rst2(4),rst(4),CGrdPntD2,BLH(3),bq,dp(80)
      real*8::a0,b0,era,erb,xx(4),x2(4),x0(2),x1(2),LL(4)
	real*8,allocatable::td1(:,:),td10(:,:),td11(:,:),dtm(:,:)
	real*8,allocatable::td2(:,:),td20(:,:),td21(:,:),td00(:,:)
	real*8,allocatable::cilm(:,:,:),cilm0(:,:,:)
	real*8,allocatable::ciln(:,:,:),ciln0(:,:,:)
	integer::status=0!
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0;ae=GRS(2)
      !输入迭代终止条件：残差标准差为原格值标准的dk倍
      !input the iteration condition parameters
      !Iteration termination condition: The standard deviation of the residual grid value is less than dk of the standard deviation
      !of the original grid value, or the difference of the residual standard deviation of the previous step iteration relative to 
      !the current step iteration is less than itd of the standard deviation of the original grid values.
      dk=1.d-3; itd=1.d-4
      !knd=0 landwater, =1 sea level variation =-1 surface atmosphere variation
      knd=1!knd=0陆地水，=1海平面变化，=-1地面大气压
	BLH(3)=0.d0;nk=0 !nk为程序批量构造的负荷球谐系数模型文件数
      !打开陆海地形球坐标格网模型，用于分离陆海区域knd=0 or =1
      !The spatial resolution of the land-sea terrain grid should not be lower than that of the surface load grid.
      open(unit=8,file=dtmfl,status="old",iostat=status)
      if(status/=0)goto 904
      read(8,'(a)') line
      call PickRecord(line,kln,rec,sn)
      if(sn<6)then
          close(8);goto 904
      endif
	hd0(1:6)=rec(1:6)
	nlat0=nint((hd0(4)-hd0(3))/hd0(6))
	nlon0=nint((hd0(2)-hd0(1))/hd0(5))
	hd0(5)=(hd0(2)-hd0(1))/dble(nlon0)
	hd0(6)=(hd0(4)-hd0(3))/dble(nlat0)
	allocate(dtm(nlat0,nlon0), stat=astat(1))
	if (sum(astat(1:1)) /= 0) then
         close(8);goto 904
	endif
	do i=1,nlat0
	  read(8,*,end=705)(dtm(i,j),j=1,nlon0)
      enddo
705   close(8)
      open(unit=58,file=tidefl,status="old",iostat=status)
      if(status/=0) goto 902 
      open(unit=20,file='tdloadharmodel.dat',status="replace")
      write(20,'(a)')trim('Global tidal load normalized spherical harmonic coefficient model in cm/hpa.')
      write(20,'(a)')trim('Created by ETideLoad, ZHANG Chuanyin, Chinese Academy of Surveying and Mapping.')
      write(20,'(a8,3a5,4a14,a11,a10,a11,a10)')'Doodson','name','n','m','Csin+','Ccos+','Csin-','Ccos-','C+','eps+','C-','eps-'
      open(unit=40,file='tideloadOne.dat',status="replace")
      write(40,'(F13.9,F13.2)')GRS(1)*1.d-14,GRS(2)
      write(40,'(a5,a8,8a16)')'name','Doodson','C10+','C10-','C11+','C11-','S11+','S11-'
	do while(.not.eof(58))
        read(58,'(a)') infl
        open(unit=8,file=infl,status="old",iostat=status)
        if(status/=0)goto 901
        read(8,'(a)') line
        call PickRecord(line,kln,rec,sn)
        if(sn<8)goto 901
	  hd(1:6)=rec(1:6);doodson=rec(7)*1.d-3
        call PickRecstr(line,kln,str,sn)
        write(nstr,*)str(8)
	  nlat=nint((hd(4)-hd(3))/hd(6))
	  nlon=nint((hd(2)-hd(1))/hd(5))
	  hd(5)=(hd(2)-hd(1))/dble(nlon)
	  hd(6)=(hd(4)-hd(3))/dble(nlat)
        maxn=nlat
        pm=nint(2.d0/hd(6))!统计时移去两极2度范围
	  if (maxn>2160) then
           close(8);goto 901
        endif
        if(nk>0) goto 1111
	  allocate(td1(nlat,nlon), stat=astat(1))
	  allocate(td10(nlat,nlon), stat=astat(2))
	  allocate(td11(nlat,nlon), stat=astat(3))
	  allocate(td2(nlat,nlon), stat=astat(4))
	  allocate(td20(nlat,nlon), stat=astat(5))
	  allocate(td21(nlat,nlon), stat=astat(6))
	  allocate(cilm(2,maxn+1,maxn+1), stat=astat(7))
	  allocate(cilm0(2,maxn+1,maxn+1), stat=astat(8))
	  allocate(ciln(2,maxn+1,maxn+1), stat=astat(9))
	  allocate(ciln0(2,maxn+1,maxn+1), stat=astat(10))
	  allocate(td00(nlat-2*pm,nlon), stat=astat(11))
	  if (sum(astat(1:11)) /= 0) then
           close(8);goto 901
	  endif
1111    td1=0.d0;td2=0.d0
	  do i=1,nlat
	    read(8,*,end=905)(td1(i,j),j=1,nlon)
          do j=1,nlon
            if(td1(i,j)>9900.0)td1(i,j)=0.d0
            tmp=CGrdPntD2(BLH(2),BLH(1),dtm,nlat0,nlon0,hd0)
            if(knd==1)then
              if(tmp>0.d0)td1(i,j)=0.d0!!zero on land
            endif
          enddo
        enddo
	  do i=1,nlat
	    read(8,*,end=905)(td2(i,j),j=1,nlon)
          do j=1,nlon
            if(td2(i,j)>9900.0)td2(i,j)=0.d0
            tmp=CGrdPntD2(BLH(2),BLH(1),dtm,nlat0,nlon0,hd0)
            if(knd==1)then
              if(tmp>0.d0)td2(i,j)=0.d0!!zero on land
            endif
         enddo
	  enddo
905	  close(8)
	  write(prfl,*)trim("testdata\"//nstr),trim(str(7)),"pro.ini"
	  write(outfl,*)trim("testdata\"//nstr),trim(str(7)),"cs.dat"
	  write(chfl,*)trim("testdata\"//nstr),trim(str(7)),"rnt.dat"
        open(unit=10,file=outfl,status="replace")
        open(unit=12,file=chfl,status="replace")
        open(unit=14,file=prfl,status="replace")
        write(10,*)'in-phase amplitude spherical harmonic coefficient model'
        write(12,101)(hd(i),i=1,6)
        write(14,*)'Iterative residual statistics of in-phase amplitude'
        call StatGrid(td1,nlat,nlon,rst1)
        bq=rst1(2);td10=td1;cilm=0;nd=0
        dp(1)=rst1(2);kk=1
        write(14,102)nd,(rst1(j),j=1,4)
        !采用一维FFT方法计算分潮负荷球谐系数
        !calculate tidal load spherical harmonic coefficients by 1-D FFT method.
908     call SphHarmExpandFFT(td10,hd,cilm0,2,maxn,nlat,nlon,GRS)
        !采用一维FFT方法计算分潮异常幅值
        !spherical harmonic synthesis of tidal constituent by 1-D FFT method.
        call SphHarmPointFFT(td11,cilm0,2,nlat,nlon,maxn,hd,GRS)
        td11=td10-td11 !residual tidal constituent
        td00=td11(pm+1:nlat-pm,1:nlon)
        call StatGrid(td00,nlat-2*pm,nlon,rst)
        write(14,102)nd+1,(rst(j),j=1,4)
        bq=dp(kk);kk=kk+1;dp(kk)=rst(2)
        if(rst(2)>rst1(2)*dk.and.bq-rst(2)>rst1(2)*itd.and.dp(kk)<dp(kk-1).and.nd<40.and.kk<42)then  !!!!!迭代条件
           cilm=cilm+cilm0;td10=td11;nd=nd+1;goto 908
        endif
        if(dp(kk)<dp(kk-1))then
          cilm=cilm+cilm0;td10=td11
        endif
        td00=td10(pm+1:nlat-pm,1:nlon)
        call StatGrid(td00,nlat-2*pm,nlon,rst)
        a0=cilm(1,1,1)*ae;era=rst(2)/rst1(2)*1.d2!zero-degree term (cm/hPa) and relative error %
        write(10,'(F13.9,F13.2,F12.4,F8.3)')GRS(1)*1.d-14,GRS(2),a0,era
 	  do n=1,maxn
	    do m=0,n
	      write(10,'(2I6,2ES24.16)')n,m,cilm(1,n+1,m+1),cilm(2,n+1,m+1)
	    enddo
	  enddo
	  do i=1,nlat
	    write(12,'(15F12.4)')(td10(i,j),j=1,nlon)
	  enddo
        write(10,*)'cross-phase amplitude spherical harmonic coefficient model'
        write(14,*)'Iterative residual statistics of cross-phase amplitude'
        call StatGrid(td2,nlat,nlon,rst1)
        bq=rst1(2);td20=td2;ciln=0;nd=0
        dp(1)=rst1(2);kk=1
        write(14,102)nd,(rst1(j),j=1,4)
        !采用一维FFT方法计算分潮负荷球谐系数
        !calculate tidal load spherical harmonic coefficients by 1-D FFT method.
708     call SphHarmExpandFFT(td20,hd,ciln0,2,maxn,nlat,nlon,GRS)
        !采用一维FFT方法计算分潮异常幅值
        !spherical harmonic synthesis of tidal constituent by 1-D FFT method.
        call SphHarmPointFFT(td21,ciln0,2,nlat,nlon,maxn,hd,GRS)
        td21=td20-td21 !residual tidal constituent
        td00=td21(pm+1:nlat-pm,1:nlon)
        call StatGrid(td00,nlat-2*pm,nlon,rst)
        write(14,102)nd+1,(rst(j),j=1,4)
        bq=dp(kk);kk=kk+1;dp(kk)=rst(2)
        if(rst(2)>rst1(2)*dk.and.bq-rst(2)>rst1(2)*itd.and.dp(kk)<dp(kk-1).and.nd<40)then  
           ciln=ciln+ciln0;td20=td21;nd=nd+1;goto 708
        endif
        if(dp(kk)<dp(kk-1))then
          ciln=ciln+ciln0;td20=td21
        endif
        td00=td20(pm+1:nlat-pm,1:nlon)
        call StatGrid(td00,nlat-2*pm,nlon,rst)
        b0=ciln(1,1,1)*ae;erb=rst(2)/rst1(2)*1.d2!zero-degree term (cm/hPa)，relative error %
        write(10,'(F13.9,F13.2,F12.4,F8.3)')GRS(1)*1.d-14,GRS(2),b0,erb!
 	  do n=1,maxn
	    do m=0,n
	      write(10,'(2I6,2ES24.16)')n,m,ciln(1,n+1,m+1),ciln(2,n+1,m+1)
	    enddo
	  enddo
	  do i=1,nlat
	    write(12,'(15F12.4)')(td20(i,j),j=1,nlon)
	  enddo
        close(10)
        close(12)
        close(14)
 	  do n=1,maxn
	    do m=0,n
            LL(1)=cilm(1,n+1,m+1);LL(2)=cilm(2,n+1,m+1)
            LL(3)=ciln(1,n+1,m+1);LL(4)=ciln(2,n+1,m+1)
            LL(1:4)=LL(1:4)*GRS(2);xx(1)=(LL(1)-LL(4))/2.d0;xx(2)=(LL(2)+LL(3))/2.d0
            xx(3)=LL(4)+xx(1);xx(4)=LL(3)-xx(2)
            x0(1)=xx(4);x0(2)=xx(3);call ReimtoHg(x0,x1)
            x2(3:4)=x1(1:2)
            x0(1)=xx(2);x0(2)=xx(1);call ReimtoHg(x0,x1)
            x2(1:2)=x1(1:2)
	      write(20,103)doodson,trim(nstr),n,m,(xx(i),i=1,4),(x2(i),i=1,4)
	    enddo
	  enddo
	  write(40,104)trim(nstr),doodson,cilm(1,2,1),ciln(1,2,1),cilm(1,2,2),ciln(1,2,2),cilm(2,2,2),ciln(2,2,2)
        write(*, '(a8,a8,F11.4,F8.3,F12.4,F8.3)')trim(nstr),trim(str(7)),a0,era,b0,erb
        !kk-1迭代次数
        nk=nk+1
901     continue
      enddo
	deallocate(td1,td10,td11,td2,td20,td21,td00,cilm,cilm0,ciln,ciln0)
902   continue
      close(58);close(20);close(40)
904   continue
101   format(4F6.1,2F13.8,F14.1)
102   format(I3,4F12.4)
103   format(f8.3,a5,2I5,4F14.8,F11.6,F10.4,F11.6,F10.4)
104   format(a5,f8.3,8E16.8)
      end
!
!************************************************************************
!
      subroutine ReimtoHg(xx,x1)
!将同相/异相幅值转换为振幅/迟角
!Convert the in-phase and out-of-phase amplitude into the amplitude and phase
	integer::nn,i,j
	real*8::xx(2),x1(2),pi,RAD
!---------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      x1(1)=dsqrt(xx(1)**2+xx(2)**2)
      if(x1(1)<1.d-16)then
        x1(2)=0.d0;goto 901
      endif
      x1(1)=dsqrt(xx(2)**2+xx(1)**2)
      x1(2)=datan(xx(2)/xx(1))/RAD
      if(xx(1)>0.and.xx(2)<0)x1(2)=x1(2)+360.d0
      if(xx(1)<0)x1(2)=x1(2)+180.d0
901   return
      end
!
!******************************************************************
!
      subroutine StatGrid(grid,row,col,rst)
!-------------------------------------------------------------
      implicit none
	integer::row,col,i,j,kk
	real*8::grid(row,col),rst(4),maxt,mint,pv,err,fr
!-----------------------------------------------------------------------
	pv=0.d0;err=0.d0;maxt=-9.d28;mint=9.d28
      fr=dble(row*col);kk=0
      do i=1,row
        do j=1,col
          if(grid(i,j)>9000.d0)goto 1001
          kk=kk+1;pv=pv+grid(i,j)/fr
	    if(maxt<grid(i,j))maxt=grid(i,j)
          if(mint>grid(i,j))mint=grid(i,j)
1001      continue
        enddo
      enddo
      pv=pv/dble(kk)*fr
      do i=1,row
        do j=1,col
          if(grid(i,j)>9000.d0)goto 1002
          err=err+(grid(i,j)-pv)**2/dble(kk)
1002      continue
        enddo
      enddo
      err=dsqrt(err)
	rst(1)=pv; rst(2)=err; rst(3)=mint;rst(4)=maxt
	return
      end
