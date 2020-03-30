
! Answers: a) cv at peak position for L=10 is 1965 and for L=7 is 575
! b) Value of energy per spin for L=7 at T=3.8 is -1.95
! c) T_p for L=7,8,10 is 4.3, 4.36, 4.4 respectively.

program ising
    implicit none
    !Define the number of iterations along with equilibration steps
    integer :: i,j,k,L=10,p,a,b,c,d,x,y,niter=1010000,time,mm,nn,oo,N
    ! Ei and Ef refer to initial and final energy after spin flip. dE refers to change in energy
    ! M is the instantaneous magnetization, cv is the specific heat at constant volume, chi is the magnetic susceptibility
    real*8 :: r,h,q,E,M,mag,Ei,Ef,dE,u
    real*8 :: T,J_ising=1.0
    real*8 :: avg_m,avg_E,avg_m_N,avg_E_N,avg_m2,avg_E2,cv,chi
    integer, dimension(:,:,:), allocatable :: spin
    integer :: T_temp,nequil=10000,cnt=0


    allocate(spin(L,L,L))
    E=0.0
    M=0.0
    N=L*L*L

    !Start with a random spin configuration on the lattice
    call random_seed

    do i=1,L
        do j=1,L
            do k=1,L
                call RANDOM_NUMBER(r)


                if(r<0.5)then
                    spin(k,j,i)=-1
                else
                    spin(k,j,i)=+1
                end if

            end do
        end do
    end do


    do i=1,L
        do j=1,L
            do k=1,L
                ! Ensure periodic boundary conditions
                a=i+1;b=i-1;c=j+1;d=j-1;x=k+1;y=k-1
                if(i==1)b=L
                if(i==L)a=1
                if(j==1)d=L
                if(j==L)c=1
                if(k==1)y=L
                if(k==L)x=1
                M=M+spin(i,j,k)
                E=E-J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
            end do
        end do
    end do
    ! Energy is halved as double contributions are taken during the calculation of the energy in the previous step
    E=E*0.5d0
    ! Start the temperature decrement loop. Note, the decrement counter has to be an integer.
    do T_temp=470,380,-2
        T=dfloat(T_temp)/100.0d0
        avg_m=0.0d0; avg_E=0.0d0; avg_m_N=0.0d0; avg_E_N=0.0d0; avg_m2=0.0d0; avg_E2=0.0d0
        ! In the above step, we initialize the average energy, average magnetization, average energy per particle, average magnetization per particle, average energy squared, average magnetization squared
        do time=1,niter
            do mm=1,L
                do nn=1,L
                    do oo=1,L
                        call RANDOM_NUMBER(r); i=int(r*dfloat(L)+1)
                        call RANDOM_NUMBER(r); j=int(r*dfloat(L)+1)
                        call RANDOM_NUMBER(r); k=int(r*dfloat(L)+1)
                        ! Ensure periodic boundary conditions
                        a=i+1;b=i-1;c=j+1;d=j-1;x=k+1;y=k-1
                        if(i==1)b=L
                        if(i==L)a=1
                        if(j==1)d=L
                        if(j==L)c=1
                        if(k==1)y=L
                        if(k==L)x=1
                        M=M+spin(i,j,k)
                        Ei=-J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                        ! Flip spin
                        spin(i,j,k)=-spin(i,j,k)
                        Ef=-J_ising*dfloat(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                        dE=Ef-Ei
                        ! Main step of Metropolis algorithm
                        if(dE<=0.0)then
                            E=E+dE
                            M=M+(2.0*dfloat(spin(i,j,k)))
                        else
                            u=exp(-dE/T)
                            call RANDOM_NUMBER(h)
                            if(h<u)then
                                E=E+dE
                                M=M+(2.0*dfloat(spin(i,j,k)))
                            else
                                spin(i,j,k)=-spin(i,j,k)
                            end if
                        end if
                    end do
                end do
            end do
            ! Start after equilibration
            if(time.gt.nequil)then
                mag=abs(M)/(dfloat(N))
                avg_m=avg_m+mag; avg_E=avg_E+E/(dfloat(N))
                avg_m_N=avg_m_N+abs(M); avg_E_N=avg_E_N+E
                avg_m2=avg_m2+(M*M); avg_E2=avg_E2+(E*E)
            end if


        end do
        avg_m=avg_m/dfloat(niter-nequil);avg_E=avg_E/dfloat(niter-nequil)
        avg_m_N=avg_m_N/dfloat(niter-nequil);avg_E_N=avg_E_N/dfloat(niter-nequil)
        avg_m2=avg_m2/dfloat(niter-nequil);avg_E2=avg_E2/dfloat(niter-nequil)
        cv=(avg_E2-(avg_E_N*avg_E_N))/(T*T)
        chi=(avg_m2-(avg_m_N*avg_m_N))/(T)
        cnt=cnt+1
        ! Print the results
        print*,cnt,avg_m,avg_E,cv,chi

    end do


end program ising






