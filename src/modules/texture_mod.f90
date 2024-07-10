module texture_mod
  use kinds
  use print_utils_mod, only : printToScreen
  implicit none
contains

! subroutine update_rotation_matrix()
!   use global, only : npts1, npts2, npts3_rank, phase_material_array, phase_material_ptr, tdot, common_grid_data, jphase, igas
!   use tensor_math_mod, only : getSkewPart
!   use bunge_mod, only : RodriguesRotationFromDeltaRotationMatrix
!   use all_grid_data_mod, only : all_grid_data
!   implicit none
!   integer :: ix, iy, iz
!   real(k_real), dimension(3,3) :: Rdot_vel_grad, &
!                                   Rdot_inelastic, &
!                                   DeltaR, &
!                                   RodriguesRotation
!   real(k_Real), pointer, dimension(:,:,:,:,:) :: R_crystal2sample=>null(), &
!                                                  vel_grad=>null()
!   real(k_real) :: ph, th, tm
!   nullify(R_crystal2sample, vel_grad)

!   ! get pointer
!   call common_grid_data%getTensor2DataPointerByName("R_crystal2sample", R_crystal2sample)
!   call common_grid_data%getTensor2DataPointerByName("velocity_gradient", vel_grad)


!   do  iz = 1,int(npts3_rank)
!      do  iy = 1,npts2
!         do  ix = 1,npts1
!           if (.not.(igas(jphase(ix,iy,iz)))) then
!             call phase_material_array%getElementPtr(jphase(ix,iy,iz), phase_material_ptr)
!             ! get the rotation contribution from the velocity gradient
!             Rdot_vel_grad = getSkewPart(vel_grad(:,:,ix,iy,iz))

!             ! get the rotation contribution from inelastic strains
!             call phase_material_ptr%getInealsticRotationRateContributionAtMaterialPoint(Rdot_inelastic, ix,iy,iz)

!             ! compute the total Delta rotation
!             DeltaR = ( Rdot_vel_grad - Rdot_inelastic)*tdot

!             !delta R to compute a rotation matrix update using rodrigues formula
!             call RodriguesRotationFromDeltaRotationMatrix(DeltaR, RodriguesRotation)

!             ! update the actual rotation matrix
!             R_crystal2sample(:,:,ix,iy,iz) = matmul(RodriguesRotation,R_crystal2sample(:,:,ix,iy,iz))

!             call bungeAnlgesFromRotationMatrixCrystal2Sample(R_crystal2sample(:,:,ix,iy,iz), ph, th, tm )
!             call bungeRotationMatrixCrystal2Sample(ph, th, tm, R_crystal2sample(:,:,ix,iy,iz))
!           endif
!         end do
!      end do
!   end do
!   nullify(R_crystal2sample, vel_grad)
! end subroutine


! subroutine update_stiffness(update_average_stiffness)
!   use global, only : npts1, npts2, npts3_rank,  phase_material_array, grid_data, c066, c0, s0
!   use tensor_math_mod, only : rotateTensor4, getSymmetricPart
!   use print_utils_mod, only : printToScreen
!   use change_tensor_basis, only : chg_basis_tensor4_to_matrix66, chg_basis_matrix66_to_tensor4
!   use all_grid_data_mod, only : all_grid_data
!   use voigt_indicial_conversion_mod
!   use stiffness_base_mod, only : stiffness_isotropic
!   use matrix_inversion, only : matrixInverseSymmetric
!   use mpi_useful_routines_mod, only : MPISumMatrix, AmIMPIMaster
!   implicit none
!   logical update_average_stiffness
!   real(k_real) :: voxel_weight
!   real(k_real) :: stiffness66_temp(6,6), s066(6,6), c0Voigt(6,6), c066_proc(6,6), c0Voigt_proc(6,6)
!   real(k_real) :: rotated_stiffness(3,3,3,3), c0VoigtTemp(6,6)
!   real(k_real), pointer :: R_crystal2sample(:,:,:,:,:)=>null(), &
!                            stiffness_ptr(:,:,:,:,:)=>null(), &
!                            effective_porosity_grid(:,:,:)=>null(), &
!                            phase_stiffness_crystal_axes(:,:,:,:,:) => null()
!   integer, pointer :: phase(:,:,:)=>null()
!   integer :: ix, iy, iz
!   nullify(R_crystal2sample, stiffness_ptr, effective_porosity_grid, phase )

!   call grid_data%getTensor2DataPointerByName("R_crystal2sample", R_crystal2sample)
!   call grid_data%getMatrix66DataPointerByName("stiffness", stiffness_ptr)
!   call grid_data%getScalarDataPointerByName("effective_porosity", effective_porosity_grid)
!   call grid_data%getScalarIntegerDataPointerByName("phase_id", phase)
!   voxel_weight = grid_data%voxel_weight

!   call phase_material_array%updatePhaseParameters()
!   call phase_material_array%getStiffnessArray(phase_stiffness_crystal_axes)
!   if (update_average_stiffness) c066_proc = 0._k_real
!   if (update_average_stiffness) c0Voigt_proc = 0._k_real
!   do  iz = 1,int(npts3_rank)
!     do  iy = 1,npts2
!       do  ix = 1,npts1
!         associate(p=>effective_porosity_grid(ix,iy,iz), &
!                   ph=>phase(ix,iy,iz), &
!                   R=>R_crystal2sample(:,:,ix,iy,iz), &
!                   stiffness=>stiffness_ptr(:,:,ix,iy,iz) )
!           call rotateTensor4(phase_stiffness_crystal_axes(:,:,:,:,ph),R, rotated_stiffness)
!           c0VoigtTemp = Tensor4ToMatrixVoigt(rotated_stiffness)
!           call chg_basis_tensor4_to_matrix66(rotated_stiffness, stiffness66_temp)
!           stiffness = getSymmetricPart(stiffness66_temp*(1._k_real-p))
!           if (update_average_stiffness) c066_proc = c066_proc + stiffness*voxel_weight
!           if (update_average_stiffness) c0Voigt_proc = c0Voigt_proc + c0VoigtTemp*voxel_weight
!       end associate
!       end do
!     end do
!   end do

!   deallocate(phase_stiffness_crystal_axes)
!   nullify(phase_stiffness_crystal_axes)


!   call MPISumMatrix(c066_proc, C066)
!   call MPISumMatrix(c0Voigt_proc, c0Voigt)

!   if (AmIMPIMaster()) then
!   if (update_average_stiffness) call printToScreen(c066, "average stiffness C066 NEW")
!   if (update_average_stiffness) call printToScreen(c0Voigt, "average stiffness c0Voigt NEW")
!   endif
!   nullify(R_crystal2sample, stiffness_ptr, effective_porosity_grid, phase)

!   call matrixInverseSymmetric(c066,s066)
!   call chg_basis_matrix66_to_tensor4(c066, c0)
!   call chg_basis_matrix66_to_tensor4(s066, s0)

! end subroutine

end module
