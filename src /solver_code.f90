
module precision_module
    use iso_fortran_env, only: real64, real32
   
   implicit none
    integer, parameter :: wp=real64
end module precision_module

module all_subroutines
    use precision_module
   
   implicit none
    contains
    
     subroutine cd_cal(n, ny, nx, dx, dy, ucn, vcn, pn, cd, re)
    	use precision_module
        implicit none
        
        integer :: n, ny, nx, i, j
        real(wp) :: re,dx, dy, fp, fv, cd, mu
        real(wp) :: ucn(0:ny+1, 0:nx+1), vcn(0:ny+1, 0:nx+1), pn(0:ny+1, 0:nx+1)
        mu=1.0_wp/re
    	fp=0.0_wp
    	 !$acc parallel loop  reduction(+:fp) present(pn,ucn) copy(fp,fv,cd)
    	do i=NINT(9.5*n)+1,NINT(10.5*n)	
    		fp=fp+(pn(i,11*n)-pn(i,12*n+1))*dy
		end do
		
		fv=0.0_wp    
    	!$acc parallel loop  reduction(+:fv) present(pn,ucn) copy(fp,fv,cd)
    	do j=(11*n)+1,(12*n)	
    		fv=fv + ((2.0_wp*mu*dx/dy)*(ucn(NINT(10.5*n+1),j)+ucn(NINT(9.5*n),j)))
		end do
    	cd=2.0_wp*(fp+fv)
    	
    end subroutine cd_cal
    
    subroutine write_snapshot(n, ny, nx, dx, dy, ucn, vcn, pn)
        use precision_module
        implicit none
        integer :: n, ny, nx, i, j
        real(wp) :: dx, dy
        real(wp) :: ucn(0:ny+1, 0:nx+1), vcn(0:ny+1, 0:nx+1), pn(0:ny+1, 0:nx+1)
        real(wp) :: x, y

        print *, "Writing data to flow_t300.dat..."
        open(unit=10, file='flow_t300.dat', status='replace')
        
        ! Header for easy identification
        write(10, *) 'X, Y, U, V, P'
        
        ! Loop through the domain to write coordinates and variables
        do j = 1, nx
            x = real(j, wp) * dx
            do i = 1, ny
                y = real(i, wp) * dy
                
                ! If inside the cylinder, write zeros to make the plot clean
                if (j >= 11*n+1 .and. j <= 12*n .and. i >= NINT(9.5*n)+1 .and. i <= NINT(10.5*n)) then
                    write(10, '(5(E15.6, 2x))') x, y, 0.0_wp, 0.0_wp, 0.0_wp
                else
                    write(10, '(5(E15.6, 2x))') x, y, ucn(i,j), vcn(i,j), pn(i,j)
                end if
            end do
            ! GNUplot reads 2D grid blocks better with a blank line between columns
            write(10, *) '' 
        end do
        
        close(10)
        print *, "Data export cnulllete."
    end subroutine write_snapshot
    
    
    Subroutine write_vorticity_snapshot(n, ny, nx, dx, dy, frame_num, ucn, vcn,pn)
    use precision_module
    implicit none
    integer :: n, ny, nx, i, j, frame_num
    real(wp) :: dx, dy
    real(wp) :: ucn(0:ny+1, 0:nx+1), vcn(0:ny+1, 0:nx+1), pn(0:ny+1, 0:nx+1)
    real(wp) :: x, y, omega
    character(len=40) :: filename

    write(filename, '(A,I6.6,A)') 'vort_', frame_num, '.dat'
    open(unit=20, file=trim(filename), status='replace')
    write(20, *) 'X, Y, OMEGA ,p'

    do j = 1, nx
        x = real(j, wp) * dx
        do i = 1, ny
            y = real(i, wp) * dy
            if (j >= 11*n+1 .and. j <= 12*n .and. i >= NINT(9.5*n)+1 .and. i <= NINT(10.5*n)) then
                write(20, '(3(E15.6, 2x))') x, y, 0.0_wp ,pn(i,j)
            else
                omega = (vcn(i,j+1)-vcn(i,j-1))/(2.0_wp*dx) - (ucn(i+1,j)-ucn(i-1,j))/(2.0_wp*dy)
                write(20, '(3(E15.6, 2x))') x, y, omega ,pn(i,j)
            end if
        end do
        write(20, *) ''
    end do
    close(20)
    end subroutine write_vorticity_snapshot

    
    subroutine uv_gp_update(n,ny,nx,uco, vco)
        use precision_module
        implicit none

        integer:: n,nx,ny
        real(wp):: uco(0:ny+1,0:nx + 1),vco(0:ny+1,0:nx + 1)
        !$acc kernels present(uco,vco)
        !$null workshare
        uco(:,0) = 2.0_wp - uco(:,1) !left u=1
        vco(:,0) = -vco(:,1) !left v=0
        
        uco(:,nx + 1) = uco(:,nx ) !right  du/dx=0
        vco(:,nx + 1) = vco(:,nx ) !right  du/dx=0

        uco(0,:)=uco(1,:) !bottom 
        vco(0,:)=-vco(1,:) !bottom 

        uco(ny+1,:)=uco(ny,:) !top
        vco(ny+1,:)=-vco(ny,:) !top
        
        !left wall square
        uco(NINT(9.5 * n)+1:NINT(10.5*n),11*n+1)=-uco(NINT(9.5*n)+1:NINT(10.5*n),11*n)
        vco(NINT(9.5*n)+1:NINT(10.5*n),11*n+1)=-vco(NINT(9.5*n)+1:NINT(10.5*n),11*n)
        
        !suqare right wall
        uco(NINT(9.5*n)+1:NINT(10.5*n),12*n)=-uco(NINT(9.5*n)+1:NINT(10.5*n),12*n+1)
        vco(NINT(9.5*n)+1:NINT(10.5*n),12*n)=-vco(NINT(9.5*n)+1:NINT(10.5*n),12*n+1)

        !suqare top wall
        uco(NINT(10.5*n),11*n+1:12*n)=-uco(NINT(10.5*n)+1,11*n+1:12*n)
        vco(NINT(10.5*n),11*n+1:12*n)=-vco(NINT(10.5*n)+1,11*n+1:12*n)
        
        !square bottom wall
        uco(NINT(9.5*n)+1,11*n+1:12*n)=-uco(NINT(9.5*n),11*n+1:12*n)
        vco(NINT(9.5*n)+1,11*n+1:12*n)=-vco(NINT(9.5*n),11*n+1:12*n)
        
        uco(NINT(9.5*n)+1,11*n+1) = 0.5_wp*(-uco(NINT(9.5*n),11*n+1) - uco(NINT(9.5*n)+1,11*n)) ! bottom-left
        vco(NINT(9.5*n)+1,11*n+1) = 0.5_wp*(-vco(NINT(9.5*n),11*n+1) - vco(NINT(9.5*n)+1,11*n))
		
		uco(NINT(9.5*n)+1,12*n) = 0.5_wp*(-uco(NINT(9.5*n),12*n) - uco(NINT(9.5*n)+1,12*n+1)) ! bottom-right
		vco(NINT(9.5*n)+1,12*n) = 0.5_wp*(-vco(NINT(9.5*n),12*n) - vco(NINT(9.5*n)+1,12*n+1))

		uco(NINT(10.5*n),11*n+1) = 0.5_wp*(-uco(NINT(10.5*n)+1,11*n+1) - uco(NINT(10.5*n),11*n)) ! top-left
		vco(NINT(10.5*n),11*n+1) = 0.5_wp*(-vco(NINT(10.5*n)+1,11*n+1) - vco(NINT(10.5*n),11*n))

		uco(NINT(10.5*n),12*n) = 0.5_wp*(-uco(NINT(10.5*n)+1,12*n) - uco(NINT(10.5*n),12*n+1)) ! top-right
		vco(NINT(10.5*n),12*n) = 0.5_wp*(-vco(NINT(10.5*n)+1,12*n) - vco(NINT(10.5*n),12*n+1)) 


          !$acc end kernels 
        
        !$null end workshare
    end subroutine uv_gp_update



    subroutine p_gp_update(n,ny,nx,po)
        use precision_module
        implicit none
        
        integer:: n,nx,ny
        real(wp):: po(0:ny+1,0:nx +1)
        
        !$acc kernels present(po)
        po(:,0) =     po(:,1) !left dp/dx=0
        po(:,nx + 1 ) = -po(:,nx) !right p0
        po(0,:)=po(1,:)!bottom dp/dy=0
        po(ny+1,:)=po(ny,:)!top dp/dy=0


        !left wall square
        po(NINT(9.5*n)+1:NINT(10.5*n),11*n+1)=po(NINT(9.5*n)+1:NINT(10.5*n),11*n)
        !suqare right wall
        po(NINT(9.5*n)+1:NINT(10.5*n),12*n)=po(NINT(9.5*n)+1:NINT(10.5*n),12*n+1)
        !suqare top wall
        po(NINT(10.5*n),11*n+1:12*n)=po(NINT(10.5*n)+1,11*n+1:12*n)
        !square bottom wall
        po(NINT(9.5*n)+1,11*n+1:12*n)=po(NINT(9.5*n),11*n+1:12*n)
        
         po(NINT(9.5*n)+1,11*n+1) = 0.5_wp*(po(NINT(9.5*n),11*n+1) + po(NINT(9.5*n)+1,11*n)) ! bottom-left
		
		po(NINT(9.5*n)+1,12*n) = 0.5_wp*(po(NINT(9.5*n),12*n) + po(NINT(9.5*n)+1,12*n+1)) ! bottom-right

		po(NINT(10.5*n),11*n+1) = 0.5_wp*(po(NINT(10.5*n)+1,11*n+1) + po(NINT(10.5*n),11*n)) ! top-left

		po(NINT(10.5*n),12*n) = 0.5_wp*(po(NINT(10.5*n)+1,12*n) + po(NINT(10.5*n),12*n+1)) ! top-right        
        
        !$acc end kernels
    end subroutine p_gp_update
    
    subroutine wall_velocity_update(n, ny,nx, uws, vws, ucs, vcs)
        use precision_module
        implicit none
        integer:: i,j,n,ny,nx
        real(wp):: uws(0:ny,0:nx), vws(0:ny ,0:nx), ucs(0:ny+1,0:nx + 1), vcs(0:ny+1,0:nx + 1)
        
        !$acc parallel loop  collapse(2) present(uws,ucs, vcs,vws)
       do j=0,nx
            do i=0,ny
                    uws(i,j)=0.5*(ucs(i+1,j+1)+ucs(i+1,j))    
                    vws(i,j)=0.5*(vcs(i,j+1)+vcs(i+1,j+1))
            end do
        end do
         !$null end do
    end subroutine wall_velocity_update
    
    subroutine pressure_poission(n, ny,nx, dx, dt, uws, vws, pn, po, be, ga_p, iter, residue, residue1, rhs_p)
        use precision_module
        implicit none
        integer:: i,j,k,n,nx,ny, iter
        real(wp) :: be, ga_p
        real(wp):: uws(0:ny,0:nx), vws(0:ny,0:nx)
        real(wp):: pn(0:ny+1,0:nx+1), po(0:ny+1,0:nx+1)
        real(wp) :: rhs_p(ny,nx), residue,residue1 , dx,dt
        iter=0
        residue=1.0_wp
        residue1=0.0_wp

        
     
      !$acc parallel loop  collapse(2)
        do j=1,11*n
                do i=1,ny
                
                rhs_p(i,j)=(dx/dt)*(uws(i-1,j)-uws(i-1,j-1)+ be*(vws(i,j-1)-vws(i-1,j-1)))
                    
                end do
            end do
            !$acc parallel loop  collapse(2)
           do j=11*n+1,12*n
                do i= 1,(NINT(9.5*n))
                    
                rhs_p(i,j)=(dx/dt)*(uws(i-1,j)-uws(i-1,j-1)+ be*(vws(i,j-1)-vws(i-1,j-1)))
                                                                
                end do
            end do    
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= (NINT(10.5*n+1)),ny
                
                rhs_p(i,j)=(dx/dt)*(uws(i-1,j)-uws(i-1,j-1)+ be*(vws(i,j-1)-vws(i-1,j-1)))
                    
                end do
            end do
            !$acc parallel loop  collapse(2)
            do j=12*n+1,nx
                do i=1,ny
                    
                rhs_p(i,j)=(dx/dt)*(uws(i-1,j)-uws(i-1,j-1)+ be*(vws(i,j-1)-vws(i-1,j-1)))
                    
                end do
            end do
          
        do while(residue >1e-6 )!.and. abs(residue1-residue)>1e-7)
        !$acc kernels
            po=pn
		!$acc end kernels
            !$acc parallel loop  collapse(2)
            do j=1,11*n
                do i=1,ny
                
                pn(i,j)=(rhs_p(i,j)-be*be*(po(i+1,j)+po(i-1,j))-(po(i,j+1)+po(i,j-1)))/ga_p
                    
                end do
            end do
            iter=iter+1
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= 1,(NINT(9.5*n))
                    
                pn(i,j)=(rhs_p(i,j)-be*be*(po(i+1,j)+po(i-1,j))-(po(i,j+1)+po(i,j-1)))/ga_p
                                                                
                end do
            end do    
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= (NINT(10.5*n+1)),ny
                
                pn(i,j)=(rhs_p(i,j)-be*be*(po(i+1,j)+po(i-1,j))-(po(i,j+1)+po(i,j-1)))/ga_p
                    
                end do
            end do
            !$acc parallel loop  collapse(2)
            do j=12*n+1,nx
                do i=1,ny
                    
                pn(i,j)=(rhs_p(i,j)-be*be*(po(i+1,j)+po(i-1,j))-(po(i,j+1)+po(i,j-1)))/ga_p
                    
                end do
            end do

            !pn(8, 8)=0.0 
            !po(8, 8)=0.0 
            residue=0.0_wp
            !$acc parallel loop  collapse(2) reduction(max:residue)
     		do j=1,nx
     			do i=1,ny
     				residue=max(residue, abs (pn(i,j)-po(i,j)))
				end do
			end do

            call p_gp_update(n,ny,nx,pn)
        end do

        print *,"exited pressure poission after", iter, "iterations with residue=", residue
        !pn(ny/2+1, nx/2+1)=0.0 
        !po(ny/2+1, nx/2+1)=0.0

        end subroutine pressure_poission
        
        
    subroutine momentum_solver(n,ny,nx, uco, vco, ucs, vcs, ucn ,vcn, uwo, vwo,uwn,vwn, dx ,dy ,dt ,re ,be ,al ,ga,iter, residue, residue1, rhs_v, rhs_u, u_temp, v_temp)
        use precision_module
        implicit none
        integer:: n,nx,ny
        real(wp) :: dx, dy, dt, re, be, al, ga
        real(wp):: uco(0:ny+1,0:nx + 1),vco(0:ny+1,0:nx + 1),ucs(0:ny+1,0:nx + 1), &
                   vcs(0:ny+1,0:nx + 1),ucn(0:ny+1,0:nx + 1),vcn(0:ny+1,0:nx + 1)
        real(wp):: uwo(0:ny,0:nx), vwo(0:ny,0:nx), uwn(0:ny,0:nx), vwn(0:ny,0:nx)
        real(wp):: cc,co,dc, residue, residue1, rhs_u(ny,nx), rhs_v(ny,nx), u_temp(0:ny+1,0:nx + 1), v_temp(0:ny+1,0:nx + 1)
        integer :: i,j,k, iter
       !$null single
       residue=1.0_wp
        residue1=1.0_wp
        iter=0
        !$null end single
        
        !$acc parallel loop  collapse(2) private(cc,co,dc)
        do j=1,11*n
              do i=1,ny
                
                cc=(0.5/dx)*(((ucn(i,j)+ucn(i,j+1))*uwn(i-1,j)-(ucn(i,j)+ucn(i,j-1))*uwn(i-1,j-1)) &
                     +((ucn(i,j)+ucn(i+1,j))*vwn(i,j-1)-(ucn(i,j)+ucn(i-1,j))*vwn(i-1,j-1)))
                co=(0.5/dx)*(((uco(i,j)+uco(i,j+1))*uwo(i-1,j)-(uco(i,j)+uco(i,j-1))*uwo(i-1,j-1)) &
                     +((uco(i,j)+uco(i+1,j))*vwo(i,j-1)-(uco(i,j)+uco(i-1,j))*vwo(i-1,j-1)))
                dc=(1.0_wp/re)*((ucn(i,j+1)+ucn(i,j-1)-2.0_wp*ucn(i,j))/(dx*dx) &
                     +(ucn(i-1,j)+ucn(i+1,j)-2.0_wp*ucn(i,j))/(dy*dy))        
                rhs_u(i,j)=(-2.0_wp*al/dt)*ucn(i,j)+al*(3.0_wp*cc-co-dc)
                
                
                cc=(0.5/dx)*(((vcn(i,j)+vcn(i,j+1))*uwn(i-1,j)-(vcn(i,j)+vcn(i,j-1))*uwn(i-1,j-1)) &
                     +((vcn(i,j)+vcn(i+1,j))*vwn(i,j-1)-(vcn(i,j)+vcn(i-1,j))*vwn(i-1,j-1)))
                co=(0.5/dx)*(((vco(i,j)+vco(i,j+1))*uwo(i-1,j)-(vco(i,j)+vco(i,j-1))*uwo(i-1,j-1)) &
                     +((vco(i,j)+vco(i+1,j))*vwo(i,j-1)-(vco(i,j)+vco(i-1,j))*vwo(i-1,j-1)))
                dc=(1.0_wp/re)*((vcn(i,j+1)+vcn(i,j-1)-2.0_wp*vcn(i,j))/(dx*dx) &
                     +(vcn(i-1,j)+vcn(i+1,j)-2.0_wp*vcn(i,j))/(dy*dy))        
                rhs_v(i,j)=(-2.0_wp*al/dt)*vcn(i,j)+al*(3.0_wp*cc-co-dc)
                
                
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2) private(cc,co,dc)
            do j=11*n+1,12*n
                do i= 1,(NINT(9.5*n))
                    
                    cc=(0.5/dx)*(((ucn(i,j)+ucn(i,j+1))*uwn(i-1,j)-(ucn(i,j)+ucn(i,j-1))*uwn(i-1,j-1)) &
                         +((ucn(i,j)+ucn(i+1,j))*vwn(i,j-1)-(ucn(i,j)+ucn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((uco(i,j)+uco(i,j+1))*uwo(i-1,j)-(uco(i,j)+uco(i,j-1))*uwo(i-1,j-1)) &
                         +((uco(i,j)+uco(i+1,j))*vwo(i,j-1)-(uco(i,j)+uco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((ucn(i,j+1)+ucn(i,j-1)-2.0_wp*ucn(i,j))/(dx*dx) &
                         +(ucn(i-1,j)+ucn(i+1,j)-2.0_wp*ucn(i,j))/(dy*dy))        
                    rhs_u(i,j)=(-2.0_wp*al/dt)*ucn(i,j)+al*(3.0_wp*cc-co-dc)
                
                
                
                    cc=(0.5/dx)*(((vcn(i,j)+vcn(i,j+1))*uwn(i-1,j)-(vcn(i,j)+vcn(i,j-1))*uwn(i-1,j-1)) &
                         +((vcn(i,j)+vcn(i+1,j))*vwn(i,j-1)-(vcn(i,j)+vcn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((vco(i,j)+vco(i,j+1))*uwo(i-1,j)-(vco(i,j)+vco(i,j-1))*uwo(i-1,j-1)) &
                         +((vco(i,j)+vco(i+1,j))*vwo(i,j-1)-(vco(i,j)+vco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((vcn(i,j+1)+vcn(i,j-1)-2.0_wp*vcn(i,j))/(dx*dx) &
                         +(vcn(i-1,j)+vcn(i+1,j)-2.0_wp*vcn(i,j))/(dy*dy))        
                    rhs_v(i,j)=(-2.0_wp*al/dt)*vcn(i,j)+al*(3.0_wp*cc-co-dc)
                
                end do
            end do    
            !$null end do nowait
            !$acc parallel loop  collapse(2) private(cc,co,dc)
            do j=11*n+1,12*n
                do i= (NINT(10.5*n+1)),ny
                    cc=(0.5/dx)*(((ucn(i,j)+ucn(i,j+1))*uwn(i-1,j)-(ucn(i,j)+ucn(i,j-1))*uwn(i-1,j-1)) &
                         +((ucn(i,j)+ucn(i+1,j))*vwn(i,j-1)-(ucn(i,j)+ucn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((uco(i,j)+uco(i,j+1))*uwo(i-1,j)-(uco(i,j)+uco(i,j-1))*uwo(i-1,j-1)) &
                         +((uco(i,j)+uco(i+1,j))*vwo(i,j-1)-(uco(i,j)+uco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((ucn(i,j+1)+ucn(i,j-1)-2.0_wp*ucn(i,j))/(dx*dx) &
                         +(ucn(i-1,j)+ucn(i+1,j)-2.0_wp*ucn(i,j))/(dy*dy))        
                    rhs_u(i,j)=(-2.0_wp*al/dt)*ucn(i,j)+al*(3.0_wp*cc-co-dc)
                
                
                
                    cc=(0.5/dx)*(((vcn(i,j)+vcn(i,j+1))*uwn(i-1,j)-(vcn(i,j)+vcn(i,j-1))*uwn(i-1,j-1)) &
                         +((vcn(i,j)+vcn(i+1,j))*vwn(i,j-1)-(vcn(i,j)+vcn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((vco(i,j)+vco(i,j+1))*uwo(i-1,j)-(vco(i,j)+vco(i,j-1))*uwo(i-1,j-1)) &
                         +((vco(i,j)+vco(i+1,j))*vwo(i,j-1)-(vco(i,j)+vco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((vcn(i,j+1)+vcn(i,j-1)-2.0_wp*vcn(i,j))/(dx*dx) &
                         +(vcn(i-1,j)+vcn(i+1,j)-2.0_wp*vcn(i,j))/(dy*dy))        
                    rhs_v(i,j)=(-2.0_wp*al/dt)*vcn(i,j)+al*(3.0_wp*cc-co-dc)
                    
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2) private(cc,co,dc)
            do j=12*n+1,nx
                do i=1,ny
                    
                    cc=(0.5/dx)*(((ucn(i,j)+ucn(i,j+1))*uwn(i-1,j)-(ucn(i,j)+ucn(i,j-1))*uwn(i-1,j-1)) &
                         +((ucn(i,j)+ucn(i+1,j))*vwn(i,j-1)-(ucn(i,j)+ucn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((uco(i,j)+uco(i,j+1))*uwo(i-1,j)-(uco(i,j)+uco(i,j-1))*uwo(i-1,j-1)) &
                         +((uco(i,j)+uco(i+1,j))*vwo(i,j-1)-(uco(i,j)+uco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((ucn(i,j+1)+ucn(i,j-1)-2.0_wp*ucn(i,j))/(dx*dx) &
                         +(ucn(i-1,j)+ucn(i+1,j)-2.0_wp*ucn(i,j))/(dy*dy))        
                    rhs_u(i,j)=(-2.0_wp*al/dt)*ucn(i,j)+al*(3.0_wp*cc-co-dc)
                
                
                
                    cc=(0.5/dx)*(((vcn(i,j)+vcn(i,j+1))*uwn(i-1,j)-(vcn(i,j)+vcn(i,j-1))*uwn(i-1,j-1)) &
                         +((vcn(i,j)+vcn(i+1,j))*vwn(i,j-1)-(vcn(i,j)+vcn(i-1,j))*vwn(i-1,j-1)))
                    co=(0.5/dx)*(((vco(i,j)+vco(i,j+1))*uwo(i-1,j)-(vco(i,j)+vco(i,j-1))*uwo(i-1,j-1)) &
                         +((vco(i,j)+vco(i+1,j))*vwo(i,j-1)-(vco(i,j)+vco(i-1,j))*vwo(i-1,j-1)))
                    dc=(1.0_wp/re)*((vcn(i,j+1)+vcn(i,j-1)-2.0_wp*vcn(i,j))/(dx*dx) &
                         +(vcn(i-1,j)+vcn(i+1,j)-2.0_wp*vcn(i,j))/(dy*dy))        
                    rhs_v(i,j)=(-2.0_wp*al/dt)*vcn(i,j)+al*(3.0_wp*cc-co-dc)

                end do
            end do
            !$null end do 
        
       ! do k=1,5000
          do while (residue >1e-6 .or. residue1>1e-6)  
          
          	!$null workshare
          	!$acc kernels
            u_temp=ucs
            v_temp=vcs
            !$acc end kernels
            !$acc parallel loop  collapse(2)
            do j=1,11*n
                do i=1,ny
                

                ucs(i,j)=(rhs_u(i,j)-be*be*(u_temp(i+1,j)+u_temp(i-1,j))-u_temp(i,j+1)-u_temp(i,j-1))/ga
                vcs(i,j)=(rhs_v(i,j)-be*be*(v_temp(i+1,j)+v_temp(i-1,j))-v_temp(i,j+1)-v_temp(i,j-1))/ga
                
                end do
            end do
            !$null end do nowait
                       !$null master
            iter=iter+1
            !$null end master
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= 1,(NINT(9.5*n))
                    
                ucs(i,j)=(rhs_u(i,j)-be*be*(u_temp(i+1,j)+u_temp(i-1,j))-u_temp(i,j+1)-u_temp(i,j-1))/ga
                vcs(i,j)=(rhs_v(i,j)-be*be*(v_temp(i+1,j)+v_temp(i-1,j))-v_temp(i,j+1)-v_temp(i,j-1))/ga
                                
                end do
            end do    
            !$null end do nowait
 
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= (NINT(10.5*n+1)),ny
                
                ucs(i,j)=(rhs_u(i,j)-be*be*(u_temp(i+1,j)+u_temp(i-1,j))-u_temp(i,j+1)-u_temp(i,j-1))/ga
                vcs(i,j)=(rhs_v(i,j)-be*be*(v_temp(i+1,j)+v_temp(i-1,j))-v_temp(i,j+1)-v_temp(i,j-1))/ga
                    
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j=12*n+1,nx
                do i=1,ny
                    
                ucs(i,j)=(rhs_u(i,j)-be*be*(u_temp(i+1,j)+u_temp(i-1,j))-u_temp(i,j+1)-u_temp(i,j-1))/ga
                vcs(i,j)=(rhs_v(i,j)-be*be*(v_temp(i+1,j)+v_temp(i-1,j))-v_temp(i,j+1)-v_temp(i,j-1))/ga
                    
                end do
            end do
            !$null end do

            call uv_gp_update(n,ny,nx,ucs, vcs)

			!$null single
			residue=0.0_wp
			residue1=0.0_wp
			!$null end single
            !$acc parallel loop  collapse(2) reduction(max:residue)
            do j=1,nx
            	do i=1,ny
            residue=max(residue, abs(ucs(i,j)-u_temp(i,j)))
            	end do
            end do
            !$null end do


            !$acc parallel loop  collapse(2) reduction(max:residue1)
            do j=1,nx
            	do i=1,ny
	            residue1=max(residue1,abs(vcs(i,j)-v_temp(i,j)))
            	end do
            end do
            !$null end do

        end do
       !$null master

       print *,"exited momentum after", iter, "iterations with residue=", residue,"residue1=", residue1
       !$null end master
       !$null barrier
    end subroutine momentum_solver


end module all_subroutines

program square_cylinder
    
    use precision_module
    use all_subroutines
    implicit none
    
    integer, parameter :: n=12
    integer :: i,j,k, iter, nstep
    integer, parameter :: ny=20*n
    integer, parameter :: nx=37*n
    
    real(wp)::dx=1.0_wp/n,dy=1.0_wp/n,dt=0.001_wp,t=0.0_wp,re=100.0_wp,be,al,ga,ga_p
    real(wp)::cu,diu,cv, div,cd
    
    real(wp):: uco(0:ny +1,0:nx +1)= 0.0_wp, vco(0:ny+1,0:nx+1) = 0.0_wp, po(0:ny+1,0:nx+1) = 0.0_wp
    real(wp):: ucs(0:ny+1,0:nx +1)= 0.0_wp, vcs(0:ny+1,0:nx+1) =0.0_wp
    real(wp):: ucn(0:ny+1,0:nx +1) =0.0_wp,vcn(0:ny+1,0:nx+1)=0.0_wp,pn(0:ny+1,0:nx+1)=0.0_wp
    
    real(wp):: uwo(0:ny,0:nx)= 0.0_wp, vwo(0:ny,0:nx) = 0.0_wp
    real(wp):: uwn(0:ny,0:nx) =0.0_wp,vwn(0:ny,0:nx)=0.0_wp
    
    real(wp) :: rhs_p(ny,nx), rhs_v(ny,nx), rhs_u(ny,nx),u_temp(0:ny+1,0:nx + 1), v_temp(0:ny+1,0:nx + 1)
    
    real(wp)::residue
    real(wp)::residue1
    real(wp)::outer_residue
    open(unit=30, file='cd_history.dat', status='replace')
    
    nstep=0
    outer_residue=1.0_wp
    residue=1.0_wp
    be=dx/dy
    al=dx*dx*re
    ga=-2.0_wp*al/dt-2.0_wp*(be*be+1.0_wp)
    ga_p=-2.0_wp*(be*be+1.0_wp)
    
    
    !$null parallel
    !$acc data copy(uco,vco,ucs,vcs,ucn,vcn,pn,po,uwo,vwo,uwn,vwn) create(u_temp, v_temp, rhs_u, rhs_v, rhs_p)
    call uv_gp_update(n,ny,nx,uco, vco)
    call uv_gp_update(n,ny,nx,ucn, vcn)
    call p_gp_update(n,ny,nx,pn)
    call p_gp_update(n,ny,nx,po)
!$acc serial present(vcn)
    vcn(NINT(10.5_wp*n)+1,12*n+1)=vcn(NINT(10.5_wp*n)+1,12*n+1)+0.05_wp
    !$acc end serial
    
    
    do while ( outer_residue > 1e-6 .and. t<200.001_wp)
    
        !$null master
        
       print *,"going into momentum"
        !$null end master
        call momentum_solver(n,ny,nx, uco, vco, ucs, vcs, ucn ,vcn, uwo, vwo,uwn,vwn, dx ,dy ,dt ,re ,be ,al ,ga,iter, residue, residue1, rhs_v, rhs_u, u_temp, v_temp)
        !$null master
        
       print *,"outside momentum"
        !$null end master
        !$null barrier
        call wall_velocity_update(n, ny,nx, uwn, vwn, ucs, vcs)
        !$null barrier
        !!$null single
        !print *,"going into pressure"
        !!$null end single
        
        call pressure_poission(n, ny,nx, dx, dt, uwn, vwn, pn, po, be, ga_p,iter, residue, residue1, rhs_p)
        !$null workshare
        !$acc kernels
        uwo=uwn
        vwo=vwn
        uco=ucn
        vco=vcn
        po=pn
        !$acc end kernels
        !$null end workshare
        !corrections step
        !$acc parallel loop  collapse(2)
        do  j=1,11*n
                do i=1,ny
                
                    ucn(i, j)=ucs(i, j)-((dt/(2.0_WP*dx))*(pn(i,j+1)-pn(i,j-1))) 
                    vcn(i, j)=vcs(i, j)-((dt/(2.0_WP*dx))*(pn(i+1,j)-pn(i-1,j))) 
                
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= 1,(NINT(9.5*n))

                    ucn(i, j)=ucs(i, j)-((dt/(2.0_WP*dx))*(pn(i,j+1)-pn(i,j-1))) 
                    vcn(i, j)=vcs(i, j)-((dt/(2.0_WP*dx))*(pn(i+1,j)-pn(i-1,j))) 
   
                end do
            end do    
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j=11*n+1,12*n
                do i= (NINT(10.5*n+1)),ny

                    ucn(i, j)=ucs(i, j)-((dt/(2.0_WP*dx))*(pn(i,j+1)-pn(i,j-1))) 
                    vcn(i, j)=vcs(i, j)-((dt/(2.0_WP*dx))*(pn(i+1,j)-pn(i-1,j))) 
    
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do  j=12*n+1,nx
                do i=1,ny
                    
                    ucn(i, j)=ucs(i, j)-((dt/(2.0_WP*dx))*(pn(i,j+1)-pn(i,j-1))) 
                    vcn(i, j)=vcs(i, j)-((dt/(2.0_WP*dx))*(pn(i+1,j)-pn(i-1,j))) 
                end do
            end do
            !$null end do
            
            !$acc parallel loop 
            do j = 0, nx
                uwn(0,j) = uwo(0,j) - (dt/dx)*(pn(0+1,j+1)-pn(0+1,j))
                   vwn(0,j) = vwo(0,j) - (dt/dy)*(pn(0+1,j+1)-pn(0,j+1))
            end do
            
            !$acc parallel loop 
            do i = 0, ny
                uwn(i,0) = uwo(i,0) - (dt/dx)*(pn(i+1,0+1)-pn(i+1,0))
                vwn(i,0) = vwo(i,0) - (dt/dy)*(pn(i+1,0+1)-pn(i,0+1))
            end do
            !$null end do
            !$acc parallel loop  collapse(2) 
            do j = 1, 11*n
                do i = 1, ny
                uwn(i,j) = uwo(i,j) - (dt/dx)*(pn(i+1,j+1)-pn(i+1,j))
                vwn(i,j) = vwo(i,j) - (dt/dy)*(pn(i+1,j+1)-pn(i,j+1))
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j = 11*n+1, 12*n-1 
                do i = 1, NINT(9.5*n)
                    uwn(i,j) = uwo(i,j) - (dt/dx)*(pn(i+1,j+1)-pn(i+1,j))
                    vwn(i,j) = vwo(i,j) - (dt/dy)*(pn(i+1,j+1)-pn(i,j+1))
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j = 11*n+1, 12*n-1
                do i = NINT(10.5*n), ny 
                    uwn(i,j) = uwo(i,j) - (dt/dx)*(pn(i+1,j+1)-pn(i+1,j))
                    vwn(i,j) = vwo(i,j) - (dt/dy)*(pn(i+1,j+1)-pn(i,j+1))
                end do
            end do
            !$null end do nowait
            !$acc parallel loop  collapse(2)
            do j = 12*n, nx 
                do i = 1, ny
                    uwn(i,j) = uwo(i,j) - (dt/dx)*(pn(i+1,j+1)-pn(i+1,j))
                    vwn(i,j) = vwo(i,j) - (dt/dy)*(pn(i+1,j+1)-pn(i,j+1))
                end do
            end do 
            !$null end do
            call uv_gp_update(n,ny,nx,ucn, vcn)
            !$null  barrier
            !$null single
			outer_residue=0.0_wp
			
			!$null end single
            !$acc parallel loop  collapse(2) reduction(max:outer_residue)
            do j=1,nx
            	do i=1,ny
            outer_residue=max(outer_residue, abs(ucn(i,j)-uco(i,j)))
            	end do
            end do
            !$null end do nowait
            
            !$null master 
             t=t+dt
             call cd_cal(n, ny, nx, dx, dy, ucn, vcn, pn, cd, re)
             write(30,'(2(E15.6,2x))')t,cd
            print*, "at t =", t,"residue  =", outer_residue, "and cd = ",cd
        	!i=null_get_num_threads()
            nstep =nstep+1
            !if(mod(nstep,50)==0) then
            !	call write_vorticity_snapshot(n, ny, nx, dx, dy, nstep, ucn, vcn)
        	!end if
            !$null end master 
            !$null barrier
            !$acc kernels
            !$null workshare
            ucs=ucn
            vcs=vcn
            uwo=uwn
            vwo=vwn
            !$acc end kernels
            !$null end workshare
			!$null  barrier
            !if(t>2.998) then
            !print *, "iam thread number ", null_get_thread_num(),"at t=",t
            !end if
            !!$null barrier
    end do
	!$acc end data    
    print *, "iam thread number ", i,"outside the while"
    !$null end parallel
    close (30)

	call write_vorticity_snapshot(n, ny, nx, dx, dy, i, ucn, vcn,pn)
    call write_snapshot(n, ny, nx, dx, dy, ucn, vcn, pn)
    
end program square_cylinder
