        program F2_parametrizations
        implicit real*8 (a-h,o-z)
        dimension q2feld(4),q2feld_lowx(4)
        common / Order / iorder
        
c      
c       gfortran F2_test.f F2_luxlike.f mstw2008structurefunctions.f 
c                alphaS.f mstwpdf.f -o F2.prg
c


        open(250,file='F2_Luxlike_MX_Q2_0_225.dat',status='unknown')
        open(251,file='F2_res_MX_Q2_0_225.dat',status='unknown')
        open(260,file='F2_Luxlike_MX_Q2_1_25.dat',status='unknown')
        open(261,file='F2_res_MX_Q2_1_25.dat',status='unknown')
        open(270,file='F2_Luxlike_MX_Q2_2_5.dat',status='unknown')
        open(271,file='F2_res_MX_Q2_2_5.dat',status='unknown')
        open(280,file='F2_Luxlike_MX_Q2_4_5.dat',status='unknown')
        open(281,file='F2_res_MX_Q2_4_5.dat',status='unknown')
        open(350,file='F2_Luxlike_Q2_0_225.dat',status='unknown')
        open(351,file='F2_res_Q2_0_225.dat',status='unknown')
        open(360,file='F2_Luxlike_Q2_1_25.dat',status='unknown')
        open(361,file='F2_res_Q2_1_25.dat',status='unknown')
        open(370,file='F2_Luxlike_Q2_2_5.dat',status='unknown')
        open(371,file='F2_res_Q2_2_5.dat',status='unknown')
        open(380,file='F2_Luxlike_Q2_4_5.dat',status='unknown')
        open(381,file='F2_res_Q2_4_5.dat',status='unknown')

        open(450,file='lowx_F2_Luxlike_Q2_1_5.dat',status='unknown')
        open(451,file='lowx_F2_res_Q2_1_5.dat',status='unknown')
        open(460,file='lowx_F2_Luxlike_Q2_2_5.dat',status='unknown')
        open(461,file='lowx_F2_res_Q2_2_5.dat',status='unknown')
        open(470,file='lowx_F2_Luxlike_Q2_5_0.dat',status='unknown')
        open(471,file='lowx_F2_res_Q2_5_0.dat',status='unknown')
        open(480,file='lowx_F2_Luxlike_Q2_8_5.dat',status='unknown')
        open(481,file='lowx_F2_res_Q2_8_5.dat',status='unknown')

        q2feld(1) = 0.225d0
        q2feld(2) = 1.25d0
        q2feld(3) = 2.5d0
        q2feld(4) = 4.5d0

        q2feld_lowx(1) = 1.5d0
        q2feld_lowx(2) = 2.5d0
        q2feld_lowx(3) = 5.d0
        q2feld_lowx(4) = 8.5d0


        wmin = 1.1
        wmax = 4.d0
        nw = 300
        dw = (wmax - wmin)/dble(nw)

        w = wmin - dw
        do iw = 1,nw
        w = w + dw
        
        q2val = q2feld(1)
        xval = q2val/(q2val + w**2 - .939**2)
        call F2_fit_luxlike(xval,q2val,F2_1,FL_1)
        call F2_res(xval,q2val,F2_res_1,FL_res_1)
        q2val = q2feld(2)
        xval = q2val/(q2val + w**2 - .939**2)
        call F2_fit_luxlike(xval,q2val,F2_2,FL_2)
        call F2_res(xval,q2val,F2_res_2,FL_res_2)
        q2val = q2feld(3)
        xval = q2val/(q2val + w**2 - .939**2)
        call F2_fit_luxlike(xval,q2val,F2_2,FL_2)
        call F2_res(xval,q2val,F2_res_3,FL_res_3)
        q2val = q2feld(4)
        xval = q2val/(q2val + w**2 - .939**2)
        call F2_fit_luxlike(xval,q2val,F2_2,FL_2)
        call F2_res(xval,q2val,F2_res_4,FL_res_4)

        write(250,*) w, F2_1,FL_1
        write(260,*) w, F2_2,FL_2
        write(270,*) w, F2_3,FL_3
        write(280,*) w, F2_4,FL_4

        write(251,*) w, F2_res_1,FL_res_1
        write(261,*) w, F2_res_1,FL_res_1
        write(271,*) w, F2_res_1,FL_res_1
        write(281,*) w, F2_res_1,FL_res_1

        enddo

        xmin = 0.001d0
        xmax = 1.d0
        nx = 300
        dx = (xmax - xmin)/dble(nx)

        x = xmin - dx
        do ix = 1,nx
        x = x + dx

        q2val = q2feld(1)
        call F2_fit_luxlike(x,q2val,F2_1,FL_1)
        call F2_res(x,q2val,F2_res_1,FL_res_1)
        q2val = q2feld(2)
        call F2_fit_luxlike(x,q2val,F2_2,FL_2)
        call F2_res(x,q2val,F2_res_2,FL_res_2)
        q2val = q2feld(3)
        call F2_fit_luxlike(x,q2val,F2_3,FL_3)
        call F2_res(x,q2val,F2_res_3,FL_res_3)
        q2val = q2feld(4)
        call F2_fit_luxlike(x,q2val,F2_4,FL_4)
        call F2_res(x,q2val,F2_res_4,FL_res_4)
        write(350,*) x, F2_1, FL_1
        write(360,*) x, F2_2, FL_2
        write(370,*) x, F2_3, FL_3
        write(380,*) x, F2_4, FL_4

        write(351,*) x, F2_res_1, FL_res_1
        write(361,*) x, F2_res_2, FL_res_2
        write(371,*) x, F2_res_3, FL_res_3
        write(381,*) x, F2_res_4, FL_res_4


        enddo
 

        ymin = -12.d0
        ymax = 0.d0
        ny = 300
        dy = (ymax - ymin)/dble(ny)

        y = ymin - dy
        do iy = 1,ny
        y = y + dy
        x = dexp(y)

        q2val = q2feld_lowx(1)
        call F2_fit_luxlike(x,q2val,F2_1,FL_1)
        call F2_res(x,q2val,F2_res_1,FL_res_1)
        q2val = q2feld_lowx(2)
        call F2_fit_luxlike(x,q2val,F2_2,FL_2)
        call F2_res(x,q2val,F2_res_2,FL_res_2)
        q2val = q2feld_lowx(3)
        call F2_fit_luxlike(x,q2val,F2_3,FL_3)
        call F2_res(x,q2val,F2_res_3,FL_res_3)
        q2val = q2feld_lowx(4)
        call F2_fit_luxlike(x,q2val,F2_4,FL_4)
        call F2_res(x,q2val,F2_res_4,FL_res_4)
        write(450,*) x, F2_1, FL_1
        write(460,*) x, F2_2, FL_2
        write(470,*) x, F2_2, FL_3
        write(480,*) x, F2_4, FL_4

        write(451,*) x, F2_res_1, FL_res_1
        write(461,*) x, F2_res_2, FL_res_2
        write(471,*) x, F2_res_3, FL_res_3
        write(481,*) x, F2_res_4, FL_res_4



        enddo
 
        end
c       ------------------------------------------------------

