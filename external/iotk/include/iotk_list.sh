
function LIST_TYPE () {
  local arg
  local type_name
  local type_declare=""
  for arg
  do
    case "$arg" in 
    (--type=*)         type_name="${arg#--type=}" ;;
    (--type-declare=*) type_declare="${arg#--type-declare=}" ;;
    esac
  done
  [ "$type_declare" ] || type_declare="type ($type_name)"
  echo "
type ${type_name}_list
  type (${type_name}_list), pointer :: next
  $type_declare,            pointer :: ptr
end type ${type_name}_list
"
}

function LIST_IMPLEMENTATION () {
  local arg
  local type_name
  local type_declare=""
  for arg
  do
    case "$arg" in
    (--type=*)             type_name="${arg#--type=}" ;;
    (--type-declare=*)     type_declare="${arg#--type-declare=}" ;;
    (--search=*)           search="${arg#--search=}" ;;
    esac
  done
  [ "$type_declare" ] || type_declare="type ($type_name)"
echo "
subroutine ${type_name}_list_init(list)
  type (${type_name}_list), intent(out) :: list
  nullify(list%ptr)
  nullify(list%next)
end subroutine ${type_name}_list_init

subroutine ${type_name}_list_destroy(list)
  type (${type_name}_list), intent(inout) :: list
  type (${type_name}_list), pointer       :: this,next
  if(.not.associated(list%next)) return
  this=>list%next
  do
    if(associated(this%ptr))deallocate(this%ptr)
    next=>this%next
    deallocate(this)
    if(.not.associated(next)) exit
    this=>next
  end do
end subroutine ${type_name}_list_destroy

subroutine ${type_name}_list_add(list,ptr)
  type (${type_name}_list), intent(inout) :: list
  $type_declare,            pointer       :: ptr
  type (${type_name}_list), pointer       :: this
  allocate(this)
  this%next => list%next
  list%next => this
  allocate(this%ptr)
  ptr => this%ptr
end subroutine ${type_name}_list_add

subroutine ${type_name}_list_del(list,ptr)
  type (${type_name}_list), intent(inout) :: list
  $type_declare,            pointer       :: ptr
  type (${type_name}_list), pointer       :: this,next_save
  if(.not.associated(list%next)) return
  if(associated(list%next%ptr,ptr)) then
    deallocate(list%next%ptr)
    next_save => list%next%next
    deallocate(list%next)
    list%next => next_save
    nullify(ptr)
    return
  end if
  this => list%next
  do
    if(.not.associated(this%next)) return
    if(associated(this%next%ptr,ptr)) exit
    this => this%next
  end do
  deallocate(this%next%ptr)
  next_save => this%next%next
  deallocate(this%next)
  this%next => next_save
  nullify(ptr)
end subroutine ${type_name}_list_del
"

if [ "$search" ] ; then
search="${search//::/@}"
IFS=";"
args=
for ss in $search ; do
  args="$args${ss#*@},"
done
args="${args%,}"

 echo "
 subroutine ${type_name}_list_search(list,ptr,$args)
   type (${type_name}_list), intent(in) :: list
   $type_declare,            pointer    :: ptr
 "
 for ss in $search ; do echo "  ${ss%@*}, optional,intent(in) :: ${ss#*@}" ; done
 echo "
   type (${type_name}_list), pointer    :: this
   nullify(ptr)
   this => list%next
   if(.not.associated(this)) return
   do
     if(.not.associated(this%ptr)) goto 1000
"
  for ss in $search ; do echo "
     if(present(${ss#*@})) then
       if(this%ptr%${ss#*@} /= ${ss#*@}) goto 1000
     end if
" ; done
echo "
     ptr => this%ptr
     exit
1000 continue
     if(.not.associated(this%next)) exit
     this => this%next
   end do
 end subroutine ${type_name}_list_search
 "

fi

}

