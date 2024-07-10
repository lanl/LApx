module gb_normals_mod
#include "macro_debug.fpp"
     implicit none
   contains
     subroutine initGBproperties(griddata, grain_id_full)
       use all_grid_data_mod, only :all_grid_data
       use global, only : kperiodic_image, all_mighty_grid_global
       use kinds
       use log_file_mod, only : writeToScreen
       implicit none
       type(all_grid_data), intent(inout) :: griddata
       integer, intent(in) :: grain_id_full(:,:,:)
       integer, allocatable :: neibs_full(:,:,:,:), idgb_full(:,:,:)
       real(k_real), allocatable :: tridistance_full(:,:,:)
       integer :: nx,ny,nz
       __SUPPRESS_CLASS_UNUSED__(griddata)
       ! GRAIN BOUNDARY RELATED STAFF, USED ONLY FOR VACANCY DIFFUSION MODEL ATM

       call all_mighty_grid_global%AMGgetGlobalGridDimension(nx,ny,nz)
       allocate(idgb_full(nx,ny,nz), &
                tridistance_full(nx,ny,nz), &
                neibs_full(7,nx,ny,nz))
   
       kperiodic_image=1
       call GBenum(grain_id_full, idgb_full, neibs_full)
       call triline(idgb_full, tridistance_full)
       call computeGBNormalsAAK(tridistance_full, grain_id_full, idgb_full, neibs_full) !-> this one computes normals need to change the name
       deallocate(idgb_full, tridistance_full, neibs_full)
       call writeToScreen('GRAIN BOUNDARY NORMALS AND TRIPLE LINES INITIALIZED')
   
     end subroutine
   
   subroutine GBenum(grain_id_full, idgb_full, neibs_full)
      use global, only : kperiodic_image, idgb, neibs, n_gb_voxels, all_mighty_grid_global
      use kinds
      implicit none
      integer, intent(in) :: grain_id_full(:,:,:)
      integer, intent(out) :: idgb_full(:,:,:), neibs_full(:,:,:,:)
      integer :: ix,iy,iz,c,a,ll,mm,n
      integer :: xn,yn,zn,sList(0:26),checklist(0:26),pList(0:26), neighbors
      integer,allocatable :: counts(:,:,:,:)
      logical :: not_skip
      integer :: npts1, npts2, npts3
      integer :: x_start, x_end, y_start, y_end, z_start, z_end
      integer :: x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank

      call all_mighty_grid_global%AMGgetGlobalGridDimension(npts1, npts2, npts3)
      call all_mighty_grid_global%AMGGetLoopLimitsGlobal(x_start, x_end, y_start, y_end, z_start, z_end)
      call all_mighty_grid_global%AMGGetRankBoxLimit(x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank)

      allocate(counts(20,npts1,npts2,npts3))
   
      counts(1,:,:,:)=0
      counts(2:20,:,:,:)=-2
    
    !!establish the kind of grain boundary point (plane, TJ, QP)!!
      n_gb_voxels = -1
    do  iz = 1,npts3  !klimit(3)
    do  iy = 1,npts2  !klimit(2)
      do  ix = 1,npts1  !klimit(1)
    
         neighbors=0
    
         if (grain_id_full(ix,iy,iz).ge.0) then
            checklist=-2
            sList=-2
            pList=1
            sList(0)=grain_id_full(ix,iy,iz)
            pList(0)=0 !igas(phase_id_full(ix,iy,iz))
            c=1
            do ll=-1,1
            do mm=-1,1
            do n=-1,1
               if ((abs(ll)+abs(mm)+abs(n)).ne.0) then
                  not_skip=.true.
                  if (kperiodic_image.eq.1) then
                     xn = mod(ix+ll-1,npts1)+1
                     if (xn.le.0) xn=xn+npts1
                     yn = mod(iy+mm-1,npts2)+1
                     if (yn.le.0) yn=yn+npts2
                     zn = mod(iz+n-1,npts3)+1
                     if (zn.le.0) zn=zn+npts3
                  else
                     xn = ix+ll
                     yn = iy+mm
                     zn = iz+n
                     if (xn.gt.npts1.or.xn.lt.1.or.yn.gt.npts2.or.yn.lt.1.or.zn.gt.npts3.or.zn.lt.1) then
                        not_skip=.false.
                     end if
                  end if
                  if (not_skip) then
                     if ((abs(ll)+abs(mm)+abs(n)).eq.1) then
                        sList(c)=grain_id_full(xn,yn,zn)
                        pList(c)=0 !igas(phase_id_full(xn,yn,zn))
                        c=c+1
                     end if
                  end if
               end if
            end do
            end do
            end do
            do ll=0,26
               c=0
               if (sList(ll).ge.0.and.pList(ll).eq.0) then
                  do mm=0,26
                     if (sList(ll).ne.checklist(mm)) c=c+1
                  end do
                  if (c.eq.27) then
                     checklist(ll)=sList(ll)
                     neighbors=neighbors+1
                  end if
               end if
            end do
    !cdiff
           if(neighbors.gt.6) neighbors=6
    !cdiff
            counts(1,ix,iy,iz)=neighbors
    !cdiff
    !cdiff      if(neighbors.gt.4) write(*,*) 'neigh =',neighbors,'voxel',ix,iy,iz
    !cdiff
            a=0
            do mm=0,26
               if (checklist(mm).ge.0.and.pList(mm).eq.0) then
                  counts(2+a,ix,iy,iz)=checklist(mm)
                  a=a+1
               end if
            end do
         end if
         ! compute number of GB voxel
         if(counts(1,ix,iy,iz)>1) n_gb_voxels = n_gb_voxels+1
      end do
      end do
      end do
      idgb(:,:,:)=counts(1,x_box_start_rank:x_box_end_rank,y_box_start_rank:y_box_end_rank, z_box_start_rank:z_box_end_rank)
      neibs(:,:,:,:)=counts(2:8,x_box_start_rank:x_box_end_rank,y_box_start_rank:y_box_end_rank, z_box_start_rank:z_box_end_rank)
      idgb_full = counts(1,:,:,:)
      neibs_full =counts(2:8,:,:,:)
      deallocate(counts)
      return
    end

   
   
   subroutine triLine(idgb_full, tridistance_full)
       use global, only :  kperiodic_image, tridistance, all_mighty_grid_global
       use kinds
     integer, intent(in) :: idgb_full(:,:,:)
     real(k_real), intent(out) :: tridistance_full(:,:,:)
     integer :: ix,iy,iz,l,m,n,shell
     integer :: xn,yn,zn
     real(k_real) :: mindist,r_crit
   
     logical :: not_skip,not_tj,surf_not_found
     integer :: npts1, npts2, npts3
     integer :: x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank

     call all_mighty_grid_global%AMGgetGlobalGridDimension(npts1, npts2, npts3)
     call all_mighty_grid_global%AMGGetRankBoxLimit(x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank)


     r_crit=2.5
   
     do iz=1,npts3; do iy=1,npts2; do ix=1,npts1
   
        not_tj=.true.
        surf_not_found=.true.
        if (idgb_full(ix,iy,iz).ge.2) then
        tridistance_full(ix,iy,iz)=900.0
        shell=1
   !!!check if current boundary point is TJ!!!
        if (idgb_full(ix,iy,iz).ge.3) then
           not_tj=.false.
           tridistance_full(ix,iy,iz)=r_crit
        end if
        do while (not_tj.and.surf_not_found)
        do l=-1*shell,shell; do m=-1*shell,shell; do n=-1*shell,shell
           if (abs(l).eq.shell.or.abs(m).eq.shell.or.abs(n).eq.shell) then
           not_skip=.true.
           if (kperiodic_image.eq.1) then
              xn = mod(ix+l-1,npts1)+1
              if (xn.le.0) xn=xn+npts1
              yn = mod(iy+m-1,npts2)+1
              if (yn.le.0) yn=yn+npts2
              zn = mod(iz+n-1,npts3)+1
              if (zn.le.0) zn=zn+npts3
           else
              xn = ix+l
              yn = iy+m
              zn = iz+n
              if (xn.gt.npts1.or.xn.lt.1.or.yn.gt.npts2.or.yn.lt.1.or.zn.gt.npts3.or.zn.lt.1) then
              not_skip=.false.
              end if
           end if
           if (not_skip) then
   !!!   determine the scaling towards triple junction!!!
           if (idgb_full(xn,yn,zn).ge.3) then
           mindist = sqrt( int2real(l**2) + int2real(m**2)+ int2real(n**2))
           if (mindist.lt.tridistance_full(ix,iy,iz)) tridistance_full(ix,iy,iz)=mindist
           if (tridistance_full(ix,iy,iz).lt.r_crit) tridistance_full(ix,iy,iz)=r_crit
           not_tj=.false.
           end if
   
   !!!   determine scaling towards volume edge if non-periodic!!!
           if (kperiodic_image.eq.0) then
           if (xn.gt.npts1.or.xn.lt.1.or.yn.gt.npts2.or.yn.lt.1.or.zn.gt.npts3.or.zn.lt.1) then
           mindist = sqrt( int2real(l**2) + int2real(m**2)+ int2real(n**2))
           if (mindist.lt.tridistance_full(ix,iy,iz)) tridistance_full(ix,iy,iz)=mindist
           if (tridistance_full(ix,iy,iz).lt.r_crit) tridistance_full(ix,iy,iz)=r_crit
           surf_not_found=.false.
           end if
           end if
           end if
           end if
        end do; end do; end do
        !!!if TJ or surface not found, expand shell out one layer and check again!!!
        shell=shell+1
      end do !while
        else
           tridistance_full(ix,iy,iz)=0.0
        end if
     end do; end do;  end do
     tridistance = tridistance_full(x_box_start_rank:x_box_end_rank,y_box_start_rank:y_box_end_rank, z_box_start_rank:z_box_end_rank)
     return
   end
   
   subroutine computeGBNormalsAAK(tridistance_full, grain_id_full, idgb_full, neibs_full)
     use global, only : kperiodic_image, gbnormals, all_mighty_grid_global
     use kinds
     use, intrinsic :: IEEE_ARITHMETIC
     implicit none
     real(k_real), intent(in) :: tridistance_full(:,:,:)
     integer, intent(in) :: grain_id_full(:,:,:), neibs_full(:,:,:,:), idgb_full(:,:,:)
     integer :: myID,myGBID,neighID,nCounts(7),searchDist,r2search
     integer :: ix,iy,iz,l,m,n,xn,yn,zn,nn
     real(k_real) :: cMass(7,3), AllCOM(7,3),NeighCOMtot(3)
     real(k_real) :: temp(3)
     real(k_real), allocatable :: gbnormals_full(:,:,:,:)!, neinormals_full(:,:,:,:,:)
     real(k_real), parameter :: voxel_size(3)=(/1._k_real, 1._k_real, 1._k_real/)
     integer :: npts1, npts2, npts3
     integer :: x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank

     call all_mighty_grid_global%AMGgetGlobalGridDimension(npts1, npts2, npts3)
     call all_mighty_grid_global%AMGGetRankBoxLimit(x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank)

     allocate(gbnormals_full(3,npts1,npts2,npts3))
     !allocate(neinormals_full(3,6,npts1,npts2,npts3))
     do iz = 1,npts3; do iy = 1,npts2; do ix = 1,npts1
       myID=grain_id_full(ix,iy,iz)
       myGBID=idgb_full(ix,iy,iz)
       
       !Do not look for normals if I am not on a grain boundary
       if(myGBID.LE.1) cycle
       
       !Set search radius and zero coms
       searchDist=ceiling(tridistance_full(ix,iy,iz))
       r2search=searchDist**2
       cMass=0._k_real
       nCounts=0
       
       !loop over cells and find center of mass for each neighbor
       do l=-searchDist,searchDist; do m=-searchDist,searchDist; do n=-searchDist,searchDist
         
         !Only include voxels which are within search radius
         if((l*l+m*m+n*n).GT.r2search) cycle
         
         !Account for periodic boundaries
         if (kperiodic_image.eq.1) then
           xn = mod(ix+l-1,npts1)+1
           if (xn.le.0) xn=xn+npts1
           yn = mod(iy+m-1,npts2)+1
           if (yn.le.0) yn=yn+npts2
           zn = mod(iz+n-1,npts3)+1
           if (zn.le.0) zn=zn+npts3
         else
           xn = ix+l         !mod(i+l-1,dim)+1
           yn = iy+m         !mod(j+m-1,dim)+1
           zn = iz+n         !mod(k+n-1,dim)+1
           if (xn.gt.npts1.or.xn.lt.1.or.yn.gt.npts2.or.yn.lt.1.or.zn.gt.npts3.or.zn.lt.1) then
             cycle
           end if
         end if
         
         neighID=-2
         !Figure out which neighbor the voxel belongs to (myID is neighbor #1), and ignore if not a neighbor
         do nn=1,myGBID
           if(neibs_full(nn,ix,iy,iz).EQ.grain_id_full(xn,yn,zn)) neighID=nn
         end do
         if(neighID.EQ.-2) cycle
         
         !Add voxel to the CoM calculation of neighbor grain 
         cMass(neighID,:)=cMass(neighID,:)+(/l,m,n/)*voxel_size
         nCounts(neighID)=nCounts(neighID)+1
       end do; end do; end do
       
       !find phyisical COMs
       NeighCOMTot=0._k_real
       do nn=2,myGBID
         NeighCOMtot=NeighCOMtot+cMass(nn,:)
         allCOM(nn,:)=cMass(nn,:)/int2real(nCounts(nn))
       end do
       NeighCOMTot=NeighCOMtot/int2real(sum(nCounts(2:myGBID)))
       allCOM(1,:)=cMass(1,:)/int2real(nCounts(1))
       
       !average normal points from MY COM to neighbor total COM
       temp=NeighCOMtot-allCOM(1,:)
       gbnormals_full(:,ix,iy,iz)=temp/sqrt(sum(temp**2))
       
       !Individual normals point from MY COM to individual neighbor COM
       do nn=2,myGBID
         temp=allCOM(nn,:)-allCOM(1,:)
         ! neinormals_full(:,nn-1,ix,iy,iz)=temp/sqrt(sum(temp**2))
       end do
       
     end do; end do; end do
     
     gbnormals = gbnormals_full(:,x_box_start_rank:x_box_end_rank,y_box_start_rank:y_box_end_rank, z_box_start_rank:z_box_end_rank)
   !   neiNormals = neinormals_full(:,:,:,:,npts3_start:npts3_end)
     deallocate(gbnormals_full)
     !deallocate(neinormals_full)
     return
   end
   
   end module
   