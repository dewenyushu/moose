        program geteulerangles
        implicit none
        character*15:: head_1,head_2,head_3,head_4,head_5,head_6,head_7
        character*15:: head_8,head_9,head_10,head_11,head_12,head_13
        character*15:: head_14,head_15,head_16,head_17,head_18,head_19
        character*15:: head_20,head_21,head_22,head_23,head_24,head_25
        character*15:: head_26,head_27,head_28,head_29,head_30,head_31
        character*15:: head_32,head_33,head_34,head_35,head_36,head_37
        character*15:: head_38,head_39,head_40,head_41,head_42,head_43
        character*15:: head_44,head_45,head_46,head_47,head_48,head_49
        character*15:: head_50,head_51,head_52,head_53,head_54,head_55
        character*15:: head_56,head_57,head_58,head_59,head_60,head_61
        character*15:: head_62,head_63,head_64,head_65,head_66,head_67
        character*15:: head_68,head_69,head_70,head_71,head_72,head_73
        character*15:: head_74,head_75,head_76,head_77,head_78,head_79
        character*15:: head_80,head_81,head_82,head_83,head_84,head_85
        character*15:: head_86,head_87,head_88,head_89,head_90,head_91
        real,dimension(101,82)::A
        real,dimension(3,3)::rotation
        real::ph,th,tm,norm_first,norm_sec,norm_third
        real::det_fp1,det_fp2,det_fp3,det_fp
        integer::i,j,linesinfile
        linesinfile = 101

        open(unit=18, file='cp_Maniaty_periodic_out.csv' , 
     &  status='old')

        read(18,*) head_1,head_2,head_3,head_4,head_5,head_6,head_7,
     &      head_8,head_9,head_10,head_11,head_12, head_13,head_14,
     &      head_15,head_16,head_17,head_18,head_19,head_20,head_21,
     &      head_22,head_23,head_24,head_26,head_27,head_28,head_29,
     &      head_30,head_31,head_32,head_33,head_34,head_35,head_36,
     &      head_37,head_38,head_39,head_40,head_41,head_42,head_43,
     &      head_44,head_45,head_46,head_47,head_48,head_49,head_50,
     &      head_51,head_52,head_53,head_54,head_55,head_56,head_57,
     &      head_58,head_59,head_60,head_61,head_62,head_63,head_64,
     &      head_65,head_66,head_67,head_68,head_69,head_70,head_71,
     &      head_72,head_73,head_74,head_75,head_76,head_77,head_78,
     &      head_79,head_80,head_81,head_82!,head_83,head_84,head_85,
!     &      head_86,head_87,head_88,head_89,head_90,head_91


       do i=1,linesinfile
       read(18,*)(A(i,j),j=1,82)
!300    format(10F11.6,A8)
       end do

        do i=1,linesinfile

              rotation(1,1)=A(i,47)
              rotation(1,2)=A(i,48) ! note that some indices are transposed due to C++ to fortran conversion
              rotation(1,3)=A(i,49)
              rotation(2,1)=A(i,50)
              rotation(2,2)=A(i,51)
              rotation(2,3)=A(i,52)
              rotation(3,1)=A(i,53)
              rotation(3,2)=A(i,54)
              rotation(3,3)=A(i,55)

!         det_fp1 = A(i,4)*(A(i,8)*A(i,12)-A(i,11)*A(i,9))
!         det_fp2 = -A(i,5)*(A(i,7)*A(i,12)-A(i,10)*A(i,9))
!         det_fp3 = A(i,6)*(A(i,7)*A(i,11)-A(i,10)*A(i,8))
!         det_fp = det_fp1 + det_fp2 + det_fp3
              
              call euler(1,ph,th,tm,rotation)

               write(71,*) 180-tm, th, 180-ph
!               write(72,*) i,det_fp

        end do

        close(18)
        end program

c *****************************************************************************
      subroutine euler (iopt,ph,th,tm,a)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
c *****************************************************************************
      implicit none
      real,dimension(3,3)::a(3,3)
      real::pi,ph,th,tm,sth,sph,cph,cth,stm,ctm
      integer iopt
      pi=4.*atan(1.d0)
 
      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        stm=sin(tm*pi/180.)
        ctm=cos(tm*pi/180.)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end

