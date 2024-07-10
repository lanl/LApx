module grid_data_var_type_mod
use string_module
implicit none

integer, parameter :: nunmber_of_var_types = 14
! the grid_data_variable_types
! ALERT!!! if you add a griddata variable type you need also to add it to the enum definition,
! to the all_grid_data_var_type_enum array and to all_grid_data_var_type_enum
public :: grid_data_var_type_enum, &
          scalar, scalar_integer, real_vector, vector5, &
          vector6, generic_vector, ss_scalar, &
          tensor2, matrix66, ss_generic_vector, &
          ss_vector5, ss_vector6, ss_generic_matrix, generic_matrix


! enum of grid_data variable types
enum, bind(c)
enumerator :: grid_data_var_type_enum = 0 !-> used as kind for utilizing enumeration
enumerator :: scalar=1, scalar_integer=2, real_vector=3, vector5=4, &
              vector6=5, generic_vector=6, ss_scalar=7, &
              tensor2=8, matrix66=9, ss_generic_vector=10, &
              ss_vector5 = 11, ss_vector6 = 12, ss_generic_matrix=13, &
              generic_matrix = 14
end enum

! a vector parameter including all variable types
integer(kind(grid_data_var_type_enum)), dimension(nunmber_of_var_types), parameter :: all_grid_data_var_type_enum= &
              (/scalar, scalar_integer, real_vector, vector5, &
              vector6, generic_vector, ss_scalar, &
              tensor2, matrix66, ss_generic_vector, &
              ss_vector5, ss_vector6, ss_generic_matrix, &
              generic_matrix /)

! this array is just for error printying error puposes, do not use it anywhere else
character(len=18), dimension(nunmber_of_var_types), parameter, private ::&
 all_grid_data_var_type_enum_names=(/&
  "scalar            ",&
  "scalar_integer    ",&
  "real_vector       ",&
  "vector5           ",&
  "vector6           ",&
  "generic_vector    ",&
  "ss_scalar         ",&
  "tensor2           ",&
  "matrix66          ",&
  "ss_generic_vector ",&
  "ss_vector5        ",&
  "ss_vector6        ",&
  "ss_generic_matrix ",&
  "generic_matrix    " &
  /)
! END OF ENUM DECLARATION

! the type managing information of a single variable.
! This type does notn contain data, bubt only variable metadata that will be used
! by a grid_data object to properly discern between variable types
type grid_data_var_type
  type(string_type), private :: var_name
  integer(kind(grid_data_var_type_enum)), private :: var_type
  integer, dimension(:), allocatable, private :: additional_var_dimensions
  integer, private :: stateful_level
contains
  procedure :: getType => getVarType
  procedure :: getName => getVarName
  procedure :: checkName => checkVarName
  procedure :: getSize => getVarSize
  procedure :: setSize => setVarSize
  procedure :: getStatefulLevel => getVarStatefulLevel
  procedure :: setStatefulLevel => setVarStatefulLevel
  procedure :: updateVarDefinition
  procedure :: init => initGridDataVarType
end type

! a type managing a list of variables. This type is only used to ad
type grid_data_var_type_array
  type(grid_data_var_type), pointer, dimension(:) :: var_type_array =>null()
  integer, private :: num_vars= 0
contains
  ! Procedure to add a new string to the list
  procedure :: addVar => varTypeArrayAddVar
  ! procedure to ask the container how many strings it holds
  procedure :: getNumVars => varTypeArrayGetNumVars
  ! procedure to retrieve a string given its position in the array
  procedure :: getVarByIndex => varTypeArrayGetVarByIndex
  ! count the number of variables by type
  procedure :: countVarByType
  ! return a variable pointer given the index
  procedure :: varTypeArrayGetVarByIndex
  ! check if a variable with the given name and type exists
  procedure :: checkVarExists
  ! check if a variable with the same name but different type exists.
  procedure :: checkVarSameNameDifferentTypeExists
  ! concatenates two lists
  procedure, private :: addVarTypeArrays
  ! add operator
  generic, public :: operator(+)   => addVarTypeArrays

end type

contains

subroutine initGridDataVarType(this, var_name, var_type, additional_var_dimensions, stateful_level)
  class(grid_data_var_type), intent(inout) :: this
  character(len=*), intent(in) :: var_name
  integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
  integer, dimension(:), optional, intent(in) :: additional_var_dimensions
  integer, intent(in), optional :: stateful_level

  if ((var_type==generic_vector).or.(var_type==ss_generic_vector).or.(var_type==generic_matrix)) then
    if (.not.(present(additional_var_dimensions))) then
      write(*,*) "var name: ", var_name
      write(*,*) "You cannot add a ", var_type, " varaible wihtout providing its size "
      write(*,*) "for a generic vector or slip-system generic vector or generic matrix the size is the number of components "
      stop
    else
      if ((var_type==generic_vector).or.(var_type==ss_generic_vector)) then
        if (size(additional_var_dimensions).ne.1) then
          write(*,*) "var name: ", var_name
          write(*,*) "additional_var_dimensions input for generic vector or slip-system  generic vector should have size==1, you provided ", additional_var_dimensions
          stop
        end if
        if (additional_var_dimensions(1)<1) then
          write(*,*) "var name: ", var_name
          write(*,*) "number of components ", additional_var_dimensions
          write(*,*) "the number of components of a for a generic vector or slip-system  generic vector vector must be > 1 "
          stop
        endif
      endif
      if ((var_type==generic_matrix)) then
        if (size(additional_var_dimensions).ne.2) then
          write(*,*) "var name: ", var_name
          write(*,*) "additional_var_dimensions input for generic matrix should have size==2, you provided ", additional_var_dimensions
          stop
        end if
        if (any(additional_var_dimensions<1)) then
          write(*,*) "var name: ", var_name
          write(*,*) "number of components ", additional_var_dimensions
          write(*,*) "the number of components of a for a generic matrix must be > 1 "
          stop
        endif
      endif
    end if
  else
    if (present(additional_var_dimensions)) then
      if (sum(additional_var_dimensions)/= 0) then
        write(*,*) "only a generic vector or slip-system  generic vector need a additional_var_dimensions value"
        stop
      end if
    end if
  end if

  ! set the variable name, type and size. The latter only if provided
  call this%var_name%setString(var_name)
  this%var_type = var_type
  if (present(additional_var_dimensions)) then
    if (allocated(this%additional_var_dimensions)) deallocate(this%additional_var_dimensions)

    allocate(this%additional_var_dimensions(size(additional_var_dimensions)))
    this%additional_var_dimensions = additional_var_dimensions
  end if

  ! set the stateful level of this variable
  if (present(stateful_level)) then
    call this%setStatefulLevel(stateful_level)
  else
    ! if no stateful level is provided we set to 2
    this%stateful_level = 2
  end if
end subroutine

! get the variable name
function getVarName(this) result(var_name)
  class(grid_data_var_type), intent(in) :: this
  character(len=:), allocatable :: var_name
  if (allocated(var_name)) deallocate(var_name)
  allocate( character(this%var_name%getStringLength()) :: var_name )
  var_name = this%var_name%getString()
end function

function checkVarName(this, str) result(is_name)
  class(grid_data_var_type), intent(in) :: this
  character(len=*), intent(in) :: str
  logical :: is_name 
  is_name = this%var_name%equal(str)
end function

! get the variable additional size information
function getVarSize(this) result(additional_var_dimensions)
  class(grid_data_var_type), intent(in) :: this
  integer, dimension(:), allocatable :: additional_var_dimensions
  if (.not.(allocated(this%additional_var_dimensions))) then
    allocate(additional_var_dimensions(1))
    additional_var_dimensions(1) = 0
  else
    if (allocated(additional_var_dimensions)) deallocate(additional_var_dimensions)
    allocate(additional_var_dimensions(size(this%additional_var_dimensions)))
    additional_var_dimensions = this%additional_var_dimensions
  end if
end function

! get the variable additional size information
subroutine setVarSize(this, additional_var_dimensions)
  class(grid_data_var_type), intent(inout) :: this
  integer, dimension(:), allocatable, intent(in) :: additional_var_dimensions
  this%additional_var_dimensions = additional_var_dimensions
end subroutine

!get the variable type, this will return the enum value
function getVarType(this) result(var_type)
  class(grid_data_var_type), intent(in) :: this
  integer(kind(grid_data_var_type_enum)) :: var_type
  var_type = this%var_type
end function

function getVarStatefulLevel(this) result(stateful_level)
  class(grid_data_var_type), intent(in) :: this
  integer :: stateful_level
  stateful_level = this%stateful_level
end function

subroutine setVarStatefulLevel(this, stateful_level)
  class(grid_data_var_type), intent(inout) :: this
  integer, intent(in) :: stateful_level
  if ((stateful_level <1).or.(stateful_level >3)) then
    write(*,*) "stateful level must be between 1 and 3, you provided ", stateful_level
    stop
  end if
  this%stateful_level = stateful_level
end subroutine

subroutine updateVarDefinition(this, additional_var_dimensions, stateful_level)
  class(grid_data_var_type), intent(inout) :: this
  integer, dimension(:), optional, intent(in) :: additional_var_dimensions
  integer, dimension(:), allocatable :: current_var_size
  integer, intent(in), optional :: stateful_level
  integer :: i
  if (present(stateful_level)) then
    if (this%getStatefulLevel() < stateful_level)  call this%setStatefulLevel(stateful_level)
  end if
  ! additional_var_dimensions is a tricky argument
  if (present(additional_var_dimensions)) then
    current_var_size = this%getSize()
    if (size(current_var_size)==size(additional_var_dimensions)) then
      do i=1,size(additional_var_dimensions)
        if (current_var_size(i) < additional_var_dimensions(i)) current_var_size(i) = additional_var_dimensions(i)
      end do
      call this%setSize(current_var_size)
    else
      write(*,*) "You can't change the dimensionality of an existing variable: orignal shape ",current_var_size, "requested shape", additional_var_dimensions
    end if
  end if
end subroutine

!********************************************************************
! grid_data_var_type_array functions
!********************************************************************

function varTypeArrayGetNumVars(this) result(num_vars)
  class(grid_data_var_type_array), intent(in) :: this
  integer :: num_vars
  num_vars = this%num_vars
end function


  ! add string to the list
  subroutine varTypeArrayAddVar(this, var_name, var_type, additional_var_dimensions, stateful_level)
    use mpi_variables_mod, only : i_am_mpi_master
    use log_file_mod, only : write_detailed_log_to_screen
    class(grid_data_var_type_array), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
    integer, dimension(:), optional, intent(in) :: additional_var_dimensions
    integer, intent(in), optional :: stateful_level
    type(grid_data_var_type), allocatable, dimension(:) :: var_type_array_temp
    integer, parameter :: allocation_batch=10
    integer :: i
    type(grid_data_var_type), pointer :: existing_var  => null()

    ! check if allocated
    if (associated(this%var_type_array)) then

      ! check a variable with the same name but different type exists
      call this%checkVarSameNameDifferentTypeExists(var_name, var_type)
      ! check if a variable with the same name and type already exists
      call this%checkVarExists(var_name, var_type, existing_var)
      if (associated(existing_var)) then
        if (i_am_mpi_master.and.write_detailed_log_to_screen) &
         write(*,*) "!//!//!// WARNING: a variable named ", var_name, " of type", &
        stringFromEnumVal(var_type), &
        " already exists. Stateful level and dimensionality will be updated !//!//!//"
        ! if the variable exists we modify it using the most conservative parameters
        if (present(stateful_level)) call existing_var%updateVarDefinition(stateful_level=stateful_level)
        if (present(additional_var_dimensions)) call existing_var%updateVarDefinition(additional_var_dimensions=additional_var_dimensions)
        return
      else
      ! if we don't have enough space we need to play the tree card game:
      ! 1 copy the strings in an temporary container:
      ! 2 deallocate the current container and allocate it with more space
      ! 3 copy back the string to the extended container
        if (size(this%var_type_array) == (this%num_vars)) then
          allocate(var_type_array_temp(this%num_vars))
          do i=1,this%num_vars
            var_type_array_temp(i)=this%var_type_array(i)
          end do
          deallocate(this%var_type_array)
          allocate(this%var_type_array(this%num_vars+allocation_batch))
          do i=1,this%num_vars
            this%var_type_array(i)=var_type_array_temp(i)
          end do
          deallocate(var_type_array_temp)
        end if
      end if
    ! if not allocated, allocate
    else
      allocate(this%var_type_array(allocation_batch))
    end if

    ! if the variable does not exists, then add it
    if (.not.associated(existing_var)) then
      ! increase the counter and store the string
      this%num_vars = this%num_vars+1
      ! the call to set the variable is depends on the parameter provided to this subroutine
      if (present(additional_var_dimensions).and.present(stateful_level)) then !-> provided all parameters
        call this%var_type_array(this%num_vars)%init( var_name, var_type, additional_var_dimensions, stateful_level)
      else if  ((.not.(present(additional_var_dimensions))).and.present(stateful_level)) then ! -> didn't provide additional_var_dimensions
        call this%var_type_array(this%num_vars)%init( var_name, var_type, stateful_level=stateful_level)
      else if  ((.not.(present(stateful_level))).and.present(additional_var_dimensions)) then ! -> didn't provide stateful_level
        call this%var_type_array(this%num_vars)%init( var_name, var_type, additional_var_dimensions=additional_var_dimensions)
      else ! -> didn't provide additional parameters
        call this%var_type_array(this%num_vars)%init( var_name, var_type)
      end if
    end if
  end subroutine

  ! check if a variable with a certain name and type exists. If it find the variable
  ! it also return a pointer to it
  subroutine checkVarExists(this, var_name, var_type, existing_var)
    class(grid_data_var_type_array), intent(in) :: this
    character(len=*), intent(in) :: var_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
    type(grid_data_var_type), pointer, intent(inout) :: existing_var
    integer :: i
    existing_var => null()

    do i=1,this%getNumVars()
      if ( this%var_type_array(i)%checkName(var_name)) then
        if (var_type == this%var_type_array(i)%getType()) then
          call this%varTypeArrayGetVarByIndex(i, existing_var)
          exit
        end if
      end if
    end do
  end subroutine

  ! check if a variable with a ceratin name already exists and have a different type
  ! if this is true we stop as we don't allow the same variable name to be used for different types
  subroutine checkVarSameNameDifferentTypeExists(this, var_name, var_type)
    class(grid_data_var_type_array), intent(in) :: this
    character(len=*), intent(in) :: var_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
    integer :: i

    do i=1,this%getNumVars()
      if (var_name == this%var_type_array(i)%getName()) then
        if (var_type /= this%var_type_array(i)%getType()) then
          write(*,*) "you requested to add ", var_name, " withy type ", stringFromEnumVal(var_type), " . "
          write(*,*) "However, a varibale named ", var_name,  " with type ", stringFromEnumVal(this%var_type_array(i)%getType()), &
           " already exists. You can't add two variable with the same name but different type!"
          stop
        end if
      end if
    end do
  end subroutine

  ! count the number of variables of a specific type within the grid_data_var_type_array
  ! optional: it also returns a list containing the variable names of the specified type
  function countVarByType(this, var_type, var_name_list) result(n_var)
    class(grid_data_var_type_array), intent(in) :: this
    integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
    type(string_array), intent(inout), optional :: var_name_list
    integer :: n_var, i
    n_var = 0

    if (present(var_name_list))  call var_name_list%reset()

    do i = 1,this%getNumVars()
      if (this%var_type_array(i)%getType()==var_type) then
        n_var=n_var+1
      if (present(var_name_list)) call var_name_list%addString(this%var_type_array(i)%getName())
    end if
    end do
  end function

  ! get a variable given its index
  subroutine varTypeArrayGetVarByIndex(this, idx, var)
    class(grid_data_var_type_array), intent(in) :: this
    type(grid_data_var_type), pointer, intent(inout) ::var
    integer, intent(in) :: idx

    if (idx > this%getNumVars()) then
      write(*,*) "index ", idx, " out of range, aborrting "
      stop
    end if
    var => this%var_type_array(idx)
  end subroutine

function stringFromEnumVal(index) result(str)
  integer, intent(in) :: index
  character(:), allocatable :: str
  str = trim(all_grid_data_var_type_enum_names(index))
end function

function addVarTypeArrays(this, other_var_type_array) result(concatenated_var_type_array)
  class(grid_data_var_type_array), intent(in) :: this, other_var_type_array
  type(grid_data_var_type_array) ::concatenated_var_type_array
  type(grid_data_var_type), pointer :: var_ptr  =>null()
  integer :: i

  if (this%getNumVars()>0) then
    do i=1,this%getNumVars()
      call this%getVarByIndex(i, var_ptr)
      call concatenated_var_type_array%addVar(var_ptr%getName(), &
                                              var_ptr%getType(), &
                                              var_ptr%getSize(), &
                                              var_ptr%getStatefulLevel() &
                                              )
    end do
  end if

  if (other_var_type_array%getNumVars()>0) then
    do i=1,other_var_type_array%getNumVars()
      call other_var_type_array%getVarByIndex(i, var_ptr)
      call concatenated_var_type_array%addVar(var_ptr%getName(), &
                                              var_ptr%getType(), &
                                              var_ptr%getSize(), &
                                              var_ptr%getStatefulLevel() &
                                              )
    end do
  end if
end function

end module
