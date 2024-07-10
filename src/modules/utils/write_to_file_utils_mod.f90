module write_to_file_utils_mod
use kinds
implicit none

contains
  subroutine writeTensorIndecesToFile(file_id, basename, new_line, my_stat_in, my_err_msg_in)
    use string_module
    implicit none
    integer, intent(in) :: file_id
    character(len=*), intent(in) :: basename
    logical, intent(in), optional :: new_line
    logical :: new_line_
    integer, intent(out), optional :: my_stat_in
    character(len=256), intent(out), optional :: my_err_msg_in
    character(len=256) :: my_err_msg
    integer :: my_stat 

    my_stat = 0
    my_err_msg = ""

    
    if (present(new_line)) then
      new_line_ = new_line
    else
      new_line_=.FALSE.
    end if
    IF (.NOT.(new_line_))  write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) ","
    write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) "eq_VM_", basename,","
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename,"_H,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_xx,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_xy,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_xz,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_yx,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_yy,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_yz,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_zx,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_zy,"
    if (my_stat.eq.0) write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) basename, "_zz"

    if (present(my_stat_in)) my_stat_in = my_stat
    if (present(my_err_msg_in)) my_err_msg_in = my_err_msg

  end subroutine

  subroutine write2DArrayToFileInline(file_id, array2D, stress_strain, new_line, my_stat_in, my_err_msg_in)
    use tensor_math_mod, only : computeVMEquivalentStrain, computeVMEquivalentStress, tensor2Norm, computePressure
    implicit none
    integer, intent(in) :: file_id
    real(k_real), dimension(:,:), intent(in) :: array2D
    logical, intent(in), optional :: new_line
    logical :: new_line_
    real(k_real) :: eq_val, H_val
    character(len=*) :: stress_strain
    integer :: i,j
    integer, intent(out), optional :: my_stat_in
    character(len=256), intent(out), optional :: my_err_msg_in
    character(len=256) :: my_err_msg
    integer :: my_stat 

    my_stat = 0
    my_err_msg = ""
    if (present(new_line)) then
      new_line_ = new_line
    else
      new_line_=.FALSE.
    end if
    IF (.NOT.(new_line_))  write(file_id, fmt='(*(A))', advance='no', iostat=my_stat, iomsg = my_err_msg) ","

    select case(stress_strain)
      case("stress")
        eq_val = computeVMEquivalentStress(array2D)
      case("strain")
        eq_val = computeVMEquivalentStrain(array2D)
      case("L2")
        eq_val = tensor2Norm(array2D)
      case default
        error stop "you can only enter stress or strain"
    end select

    H_val = computePressure(array2D)
    write(file_id, fmt='(*(G0,1A))', advance='no', iostat=my_stat, iomsg = my_err_msg) eq_val,","
    if (my_stat.eq.0)  write(file_id, fmt='(*(G0,1A))', advance='no', iostat=my_stat, iomsg = my_err_msg) H_val,","

    if (my_stat.eq.0) then
    do i = 1,3
      do j = 1,3
        if ((i*j)<9)  then
          if (my_stat.eq.0) write(file_id, fmt='(*(G0,1A))', advance='no', iostat=my_stat, iomsg = my_err_msg) array2D(i,j),","
        endif
    enddo
    enddo
    endif

    if (my_stat.eq.0) write(file_id, fmt='(*(G0))', advance='no', iostat=my_stat, iomsg = my_err_msg) array2D(3,3)


    if (present(my_stat_in)) my_stat_in = my_stat
    if (present(my_err_msg_in)) my_err_msg_in = my_err_msg

  end subroutine

  subroutine writeVectorIndecesToFile(file_id, basename, n_comp, add_trailing_coma)
    use string_module
    implicit none
    integer, intent(in) :: file_id, n_comp
    character(len=*), intent(in) :: basename
    logical, intent(in), optional :: add_trailing_coma
    logical :: add_trailing_coma_
    integer :: i

    if (present(add_trailing_coma)) then
      add_trailing_coma_ = add_trailing_coma
    else
      add_trailing_coma_=.TRUE.
    end if

    do i =1,n_comp-1
    write(file_id, fmt='(*(A))', advance='no') basename, i
    enddo 

    if (add_trailing_coma_) then
      write(file_id, fmt='(*(A))', advance='no') basename, i, ","
    else
      write(file_id, fmt='(*(A))', advance='no') basename, i
    endif

  end subroutine

  subroutine writeVectorToFileInline(file_id, vector, add_trailing_coma)
    use tensor_math_mod, only : computeVMEquivalentStrain, computeVMEquivalentStress, tensor2Norm, computePressure
    implicit none
    integer, intent(in) :: file_id
    real(k_real), dimension(:), intent(in) :: vector
    logical, intent(in), optional :: add_trailing_coma
    logical :: add_trailing_coma_
    integer :: i, n_comp

    n_comp = size(vector)
    if (present(add_trailing_coma)) then
      add_trailing_coma_ = add_trailing_coma
    else
      add_trailing_coma_=.TRUE.
    end if

    do i = 1,n_comp
        if ((i)<n_comp)   write(file_id, fmt='(*(G0,1A))', advance='no') vector(i),","
    enddo

    if (add_trailing_coma_) then
      write(file_id, fmt='(*(G0,1A))', advance='no') vector(n_comp),","
    else 
      write(file_id, fmt='(*(G0))') vector(n_comp)
    endif

  end subroutine

  subroutine writeScalarHeaderToFile(file_id, scalar_name, add_trailing_coma)
    use string_module
    implicit none
    integer, intent(in) :: file_id
    character(len=*), intent(in) :: scalar_name
    logical, intent(in), optional :: add_trailing_coma
    logical :: add_trailing_coma_

    if (present(add_trailing_coma)) then
      add_trailing_coma_ = add_trailing_coma
    else
      add_trailing_coma_=.TRUE.
    end if

    if (add_trailing_coma_) then
      write(file_id, fmt='(*(A))', advance='no') scalar_name,","
    else
      write(file_id, fmt='(*(A))') scalar_name
    endif

  end subroutine


  subroutine writeScalarToFileInline(file_id, scalar, add_trailing_coma)
    use tensor_math_mod, only : computeVMEquivalentStrain, computeVMEquivalentStress, tensor2Norm, computePressure
    implicit none
    integer, intent(in) :: file_id
    real(k_real), intent(in) :: scalar
    logical, intent(in), optional :: add_trailing_coma
    logical :: add_trailing_coma_

    if (present(add_trailing_coma)) then
      add_trailing_coma_ = add_trailing_coma
    else
      add_trailing_coma_=.TRUE.
    end if

    if (add_trailing_coma_) then
      write(file_id, fmt='(*(G0,1A))', advance='no') scalar,","
    else 
      write(file_id, fmt='(*(G0))', advance='no') scalar
    endif

  end subroutine

end module
