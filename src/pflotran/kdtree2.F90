!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.


module kdtree2_priority_queue_module

  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload',
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
#include "petsc/finclude/petscsys.h"  
!  use kdtree2_precision_module
  use petscsnes

  implicit none

  public :: kdtree2_result

  public :: pq, pq_create, pq_insert, pq_extract_max, pq_replace_max

!  integer, parameter :: sp = kind(0.0)
!  integer, parameter :: dp = kind(0.0d0)
  !DELETE
!  PetscInt :: kdkind = dp
!  public :: kdkind
  !
  ! a pair of distances, indexes
  type kdtree2_result
    PetscReal :: dis
    PetscInt :: idx 
  end type kdtree2_result

  ! The priority queue consists of elements
  ! priority(1:heap_size), with associated payload(:).
  !
  ! There are heap_size active elements.
  ! Assumes the allocation is always sufficient.  Will NOT increase it
  ! to match.
  type pq
    PetscInt :: heap_size = 0
     
!      integer :: heap_size = 0
      type(kdtree2_result), pointer :: elems(:)
  end type pq

contains

!************************************************************************** !

function pq_create(results_in) result(res)
  !
  ! Create a priority queue from ALREADY allocated
  ! array pointers for storage.  NOTE! It will NOT
  ! add any alements to the heap, i.e. any existing
  ! data in the input arrays will NOT be used and may
  ! be overwritten.
  !

  implicit none
  
  type(kdtree2_result), target:: results_in(:)
  type(pq) :: res
  PetscInt :: nalloc
    !
    !
!    integer :: nalloc

  nalloc = size(results_in,1)
  if (nalloc < 1) then
      print *, 'PQ_CREATE: error, input arrays must be allocated.'
  end if
  res%elems => results_in
  res%heap_size = 0
  return
end function pq_create

!************************************************************************** !

subroutine heapify(a,i_in)
  !
  ! take a heap rooted at 'i' and force it to be in the
  ! heap canonical form.   This is performance critical
  ! and has been tweaked a little to reflect this.
  !
  
  type(pq),pointer   :: a
  PetscInt :: i_in
  PetscInt :: i, l, r, largest
  PetscReal :: pri_i, pri_l, pri_r, pri_largest
  type(kdtree2_result) :: temp

  i = i_in

  bigloop:  do
              ! left(i)
              l = 2 * i
              ! right(i)
              r = l + 1
       !
       ! set 'largest' to the index of either i, l, r
       ! depending on whose priority is largest.
       !
       ! note that l or r can be larger than the heap size
       ! in which case they do not count.

              if (l > a%heap_size) then
                ! we know that i is the largest as both l and r are invalid.
                exit
              else
                pri_i = a%elems(i)%dis
                pri_l = a%elems(l)%dis
                if (pri_l > pri_i) then
                  largest = l
                  pri_largest = pri_l
                else
                  largest = i
                  pri_largest = pri_i
                endif

          !
          ! between i and l we have a winner
          ! now choose between that and r.
          !
                if (r < a%heap_size) then
                  pri_r = a%elems(r)%dis
                  if (pri_r > pri_largest) then
                    largest = r
                  endif
                endif
              endif

              if (largest /= i) then
          ! swap data in nodes largest and i, then heapify

                temp = a%elems(i)
                a%elems(i) = a%elems(largest)
                a%elems(largest) = temp
          !
          ! Canonical heapify() algorithm has tail-ecursive call:
          !
          !        call heapify(a,largest)
          ! we will simulate with cycle
          !
                i = largest
                cycle bigloop ! continue the loop
              else
                return   ! break from the loop
              end if
            enddo bigloop
  return
end subroutine heapify

!************************************************************************** !

subroutine pq_extract_max(a,e)
  !
  ! return the priority and payload of maximum priority
  ! element, and remove it from the queue.
  ! (equivalent to 'pop()' on a stack)
  !

  implicit none
  
  type(pq),pointer :: a
  type(kdtree2_result), intent(out) :: e

  if (a%heap_size >= 1) then
    !
    ! return max as first element
    !
    e = a%elems(1)

    !
    ! move last element to first
    !
    a%elems(1) = a%elems(a%heap_size)
    a%heap_size = a%heap_size-1
    call heapify(a,1)
    return
  else
    print *, 'PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ'
    stop
  end if

end subroutine pq_extract_max

!************************************************************************** !

function pq_insert(a,dis,idx)
  !
  ! Insert a new element and return the new maximum priority,
  ! which may or may not be the same as the old maximum priority.
  !
  implicit none
  
  type(pq),pointer  :: a
  PetscReal :: dis
  PetscInt :: idx
  PetscInt :: i, isparent
  PetscReal:: parentdis
  PetscReal :: pq_insert
    !

    !    if (a%heap_size .ge. a%max_elems) then
    !       write (*,*) 'PQ_INSERT: error, attempt made to insert element on full PQ'
    !       stop
    !    else
  a%heap_size = a%heap_size + 1
  i = a%heap_size

  do while (i > 1)
    isparent = int(i/2)
    parentdis = a%elems(isparent)%dis
    if (dis > parentdis) then
      ! move what was in i's parent into i.
      a%elems(i)%dis = parentdis
      a%elems(i)%idx = a%elems(isparent)%idx
      i = isparent
    else
      exit
    endif
  end do

  ! insert the element at the determined position
  a%elems(i)%dis = dis
  a%elems(i)%idx = idx

  pq_insert = a%elems(1)%dis
  return

end function pq_insert

!************************************************************************** !

function pq_replace_max(a,dis,idx)
  !
  ! Replace the extant maximum priority element
  ! in the PQ with (dis,idx).  Return
  ! the new maximum priority, which may be larger
  ! or smaller than the old one.
  !

  implicit none
  
  type(pq), pointer :: a
  PetscReal :: dis
  PetscInt :: idx
  PetscInt :: parent, child, N
  PetscReal :: prichild, prichildp1

  type(kdtree2_result) :: etmp
  PetscReal :: pq_replace_max

  if (PETSC_TRUE) then
    N = a%heap_size
    if (N >= 1) then
      parent = 1
      child = 2

      loop: do while (child <= N)
              prichild = a%elems(child)%dis

              ! posibly child+1 has higher priority, and if
              ! so, get it, and increment child.

              if (child < N) then
                prichildp1 = a%elems(child + 1)%dis
                if (prichild < prichildp1) then
                  child = child + 1
                  prichild = prichildp1
                endif
              endif

              if (dis >= prichild) then
                exit loop
                ! we have a proper place for our new element,
                ! bigger than either children's priority.
              else
                ! move child into parent.
                a%elems(parent) = a%elems(child)
                parent = child
                child = 2 * parent
              end if
            end do loop
      a%elems(parent)%dis = dis
      a%elems(parent)%idx = idx
      pq_replace_max = a%elems(1)%dis
    else
      a%elems(1)%dis = dis
      a%elems(1)%idx = idx
      pq_replace_max = dis
    endif
  else
       !
       ! slower version using elementary pop and push operations.
       !
    call pq_extract_max(a,etmp)
    etmp%dis = dis
    etmp%idx = idx
    pq_replace_max = pq_insert(a,dis,idx)
  endif
  return
end function pq_replace_max

end module kdtree2_priority_queue_module


module kdtree2_module

#include "petsc/finclude/petscsys.h"
  
!  use kdtree2_precision_module
  use kdtree2_priority_queue_module
 
  use petscsnes

  implicit none
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported.
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module.
  !

  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy
  public :: kdtree2_n_nearest
  public :: kdtree2_sort_results

  ! The maximum number of points to keep in a terminal node.
  PetscInt :: bucket_size = 12   

  type interval
    PetscReal :: lower, upper 
  end type interval

  ! an internal tree node
  type :: tree_node
    ! the dimension to cut
    PetscInt :: cut_dim
    ! where to cut the dimension
    PetscReal :: cut_val
    ! improved cutoffs knowing the spread in child boxes.
    PetscReal :: cut_val_left, cut_val_right
    PetscInt :: l, u
    type(tree_node), pointer :: left, right
    type(interval), pointer :: box(:) => null()
  end type tree_node

  ! Global information about the tree, one per tree
  type :: kdtree2
    ! dimensionality and total # of points
    PetscInt :: dimen = 0, n = 0
    ! pointer to the actual data array
    PetscReal, pointer :: the_data(:,:) => null()
    ! permuted index into the data, so that indexes[l..u] of some
    ! bucket represent the indexes of the actual points in that
    ! bucket.
    PetscInt, pointer :: ind(:) => null()
    PetscBool :: sort = PETSC_FALSE
    PetscBool :: rearrange = PETSC_FALSE
    ! if (rearrange .eqv. .true.) then rearranged_data has been
    ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
    ! permitting search to use more cache-friendly rearranged_data, at
    ! some initial computation and storage cost.
    PetscReal, pointer :: rearranged_data(:,:) => null()
    ! root pointer of the tree
    type (tree_node), pointer :: root => null()
  end type kdtree2

  type :: tree_search_record    
    ! One of these is created for each search.
    ! Many fields are copied from the tree structure, in order to
    ! speed up the search.
    PetscInt :: dimen, nn, nfound
    PetscReal :: ballsize
    PetscInt :: centeridx = 999, correltime = 9999
    ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
    PetscInt :: nalloc
    PetscBool :: rearrange, overflow
    PetscReal, pointer :: qv(:)
    type(kdtree2_result), pointer :: results(:)
    type(pq) :: pq
     PetscReal, pointer :: data(:,:)
    PetscInt, pointer :: ind(:)
  end type tree_search_record

  ! A GLOBAL VARIABLE for search
  type(tree_search_record), target :: sr   

contains

! ************************************************************************** !
  
function kdtree2_create(input_data,dim,sort,rearrange) result(mr)

!  use petscvec
  
  implicit none  

  type (kdtree2), pointer :: mr
  PetscInt, optional :: dim
  PetscBool, optional :: sort, rearrange
  PetscReal,target :: input_data(:,:)
  PetscInt :: i

  allocate (mr)

  mr%the_data => input_data
    
  if (present(dim)) then
    mr%dimen = dim
  else
    mr%dimen = size(input_data,1)
  end if
  mr%n = size(input_data,2)

  if (mr%dimen > mr%n) then
    !  unlikely to be correct
    print *, 'something is wrong'
    stop
  end if


  call build_tree(mr)

  if (present(sort)) then
    mr%sort = sort
  else
    mr%sort = PETSC_FALSE
  endif

  if (present(rearrange)) then
    mr%rearrange = rearrange
  else
    mr%rearrange = PETSC_TRUE
  endif

  if (mr%rearrange) then
    allocate(mr%rearranged_data(mr%dimen,mr%n))
    do i=1,mr%n
      mr%rearranged_data(:,i) = mr%the_data(:, &
      mr%ind(i))
    enddo
  else
    nullify(mr%rearranged_data)
  endif

end function kdtree2_create

! ************************************************************************** !

subroutine build_tree(tp)

  implicit none
  
  type (kdtree2), pointer :: tp
  PetscInt :: j
  type(tree_node), pointer :: dummy => null()

  allocate (tp%ind(tp%n))
  forall (j=1:tp%n)
    tp%ind(j) = j
  end forall

  tp%root => build_tree_for_range(tp,1,tp%n, dummy)
end subroutine build_tree

! ************************************************************************** !

recursive function build_tree_for_range(tp,l,u,parent) result (res)

  implicit none
  
  ! Function Return Cut_value 
  type (tree_node), pointer :: res

  ! Structure Arguments 
  type (kdtree2), pointer :: tp
  type (tree_node),pointer :: parent

  ! Scalar Arguments

  PetscInt :: l, u

  ! Local Scalars
  PetscInt :: i, c, m, dimen
  PetscBool :: recompute

  PetscReal :: average

  ! first compute min and max
  dimen = tp%dimen
  allocate (res)
  allocate(res%box(dimen))

  ! First, compute an APPROXIMATE bounding box of all points &
  ! associated with this node.
  if (u < l) then
    ! no points in this box
    nullify(res)

    return
  end if

  if ((u - l) <= bucket_size) then
         
    ! always compute true bounding box for terminal nodes.        
    do i=1,dimen
      call spread_in_coordinate(tp,i,l,u,res%box(i))
    end do

    res%cut_dim = 0
    res%cut_val = 0.0
    res%l = l
    res%u = u
    res%left => null()
    res%right => null()
  else
         
    ! modify approximate bounding box.  This will be an
    ! overestimate of the true bounding box, as we are only recomputing
    ! the bounding box for the dimension that the parent split on.
    !
    ! Going to a true bounding box computation would significantly
    ! increase the time necessary to build the tree, and usually
    ! has only a very small difference.  This box is not used
    ! for searching but only for deciding which coordinate to split on.
         
    do i=1,dimen
      recompute= PETSC_TRUE !.true.
      if (associated(parent)) then
        if (i .ne. parent%cut_dim) then
          recompute= PETSC_FALSE !.false.
        end if
      endif
      if (recompute) then
        call spread_in_coordinate(tp,i,l,u,res%box(i))
      else
        res%box(i) = parent%box(i)
      endif
    end do
  
    c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         
    ! c is the identity of which coordinate has the greatest spread.
         
!!!!??????????????????
    if (PETSC_FALSE) then
      ! select exact median to have fully balanced tree.      
      m = (l+u)/2
      call select_on_coordinate(tp%the_data,tp%ind,c,m,l,u)
    else
            
      ! select point halfway between min and max, as per A. Moore,
      ! who says this helps in some degenerate cases, or
      ! actual arithmetic average.

       !other average possible faster...!!!!!
      if (PETSC_TRUE) then
        ! actually compute average
        average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1)
      else
        average = (res%box(c)%upper + res%box(c)%lower)/2.0
      endif

      res%cut_val = average
      m = select_on_coordinate_value(tp%the_data,tp%ind,c,average,l,u)
    endif

    ! moves indexes around
    res%cut_dim = c
    res%l = l
    res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

    res%left => build_tree_for_range(tp,l,m,res)
    res%right => build_tree_for_range(tp,m+1,u,res)

    if (associated(res%right) .eqv. PETSC_FALSE) then
      res%box = res%left%box
      res%cut_val_left = res%left%box(c)%upper
      res%cut_val = res%cut_val_left
   elseif (associated(res%left) .eqv. PETSC_FALSE) then
      res%box = res%right%box
      res%cut_val_right = res%right%box(c)%lower
      res%cut_val = res%cut_val_right
    else
      res%cut_val_right = res%right%box(c)%lower
      res%cut_val_left = res%left%box(c)%upper
      res%cut_val = (res%cut_val_left + res%cut_val_right)/2


      ! now remake the true bounding box for self.
      ! Since we are taking unions (in effect) of a tree structure,
      ! this is much faster than doing an exhaustive
      ! search over all points
      res%box%upper = max(res%left%box%upper,res%right%box%upper)
      res%box%lower = min(res%left%box%lower,res%right%box%lower)
     endif
  end if
end function build_tree_for_range

! ************************************************************************** !

function select_on_coordinate_value(v,ind,c,alpha,li,ui) &
     result(res)

  ! Move elts of ind around between l and u, so that all points
  ! <= than alpha (in c cooordinate) are first, and then
  ! all points > alpha are second.

  !
  ! Algorithm (matt kennel).
  !
  ! Consider the list as having three parts: on the left,
  ! the points known to be <= alpha.  On the right, the points
  ! known to be > alpha, and in the middle, the currently unknown
  ! points.   The algorithm is to scan the unknown points, starting
  ! from the left, and swapping them so that they are added to
  ! the left stack or the right stack, as appropriate.
  !
  ! The algorithm finishes when the unknown stack is empty.

  implicit none

  PetscInt :: c, li, ui
  PetscReal :: alpha

  PetscReal :: v(1:,1:)

  PetscInt :: ind(1:), tmp, lb, rb

  PetscInt :: res
      
  ! The points known to be <= alpha are in
  ! [l,lb-1]
  !
  ! The points known to be > alpha are in
  ! [rb+1,u].
  !
  ! Therefore we add new points into lb or
  ! rb as appropriate.  When lb=rb
  ! we are done.  We return the location of the last point <= alpha.

  lb = li; rb = ui

  do while (lb < rb)
    if ( v(c,ind(lb)) <= alpha ) then
      ! it is good where it is.
      lb = lb+1
    else
      ! swap it with rb. !faster way to swap??
      tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
      rb = rb-1
    endif
  end do

  if (v(c,ind(lb)) <= alpha) then
    res = lb
  else
    res = lb-1
  endif

end function select_on_coordinate_value

! ************************************************************************** !

subroutine select_on_coordinate(v,ind,c,k,li,ui)
  ! Move elts of ind around between l and u, so that the kth
  ! element
  ! is >= those below, <= those above, in the coordinate c.

  implicit none

  PetscInt :: c, k, li, ui
  PetscInt :: i, l, m, s, t, u, ind(:)
  PetscReal ::  v(:,:)

  l = li
  u = ui
  do while (l < u)
    t = ind(l)
    m = l
    do i = l + 1, u
      if (v(c,ind(i)) < v(c,t)) then
        m = m + 1
        s = ind(m)
        ind(m) = ind(i)
        ind(i) = s
      end if
    end do
    s = ind(l)
    ind(l) = ind(m)
    ind(m) = s
    if (m <= k) l = m + 1
    if (m >= k) u = m - 1
  end do
end subroutine select_on_coordinate

! ************************************************************************** !

subroutine spread_in_coordinate(tp,c,l,u,interv)
  ! the spread in coordinate 'c', between l and u.
  !
  ! Return lower bound in 'smin', and upper in 'smax',

  implicit none
  
  type (kdtree2), pointer :: tp
  type(interval), intent(out) :: interv

  PetscInt :: c, l, u

  PetscReal :: last, lmax, lmin, t, smin, smax
  PetscInt :: i, ulocal

  PetscReal, pointer :: v(:,:)
  PetscInt, pointer :: ind(:)

  v => tp%the_data(1:,1:)
  ind => tp%ind(1:)
  smin = v(c,ind(l))
  smax = smin

  ulocal = u

  do i = l + 2, ulocal, 2
    lmin = v(c,ind(i-1))
    lmax = v(c,ind(i))
    if (lmin > lmax) then
      t = lmin
      lmin = lmax
      lmax = t
    end if
    if (smin > lmin) smin = lmin
    if (smax < lmax) smax = lmax
  end do
 
  if (i == ulocal + 1) then
    last = v(c,ind(ulocal))
    if (smin > last) smin = last
    if (smax < last) smax = last
   end if

   interv%lower = smin
   interv%upper = smax

end subroutine spread_in_coordinate

! ************************************************************************** !

subroutine kdtree2_destroy(tp)

  ! Deallocates all memory for the tree, except input data matrix
  ! .. Structure Arguments ..

  implicit none

  type(kdtree2), pointer :: tp

  call destroy_node(tp%root)

  deallocate (tp%ind)
  nullify (tp%ind)

  if (tp%rearrange) then
    deallocate(tp%rearranged_data)
    nullify(tp%rearranged_data)
  endif

  deallocate(tp)
 end subroutine kdtree2_destroy

recursive subroutine destroy_node(np)

  implicit none

  type (tree_node), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
!      intrinsic ASSOCIATED
      ! ..
  if (associated(np%left)) then
    call destroy_node(np%left)
    nullify (np%left)
  end if
  if (associated(np%right)) then
    call destroy_node(np%right)
    nullify (np%right)
  end if
  if (associated(np%box)) deallocate(np%box)
  deallocate(np)
  return

end subroutine destroy_node

! ************************************************************************** !

subroutine kdtree2_n_nearest(tp,qv,nn,results)
  
  ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
  ! returning their indexes and distances in 'indexes' and 'distances'
  ! arrays already allocated passed to this subroutine.

  implicit none
  
  type (kdtree2), pointer      :: tp
  PetscReal, target :: qv(:)
!    real(kdkind), target, intent (In)    :: qv(:)
!  integer, intent (In)         :: nn
  PetscInt :: nn
  type(kdtree2_result), target :: results(:)
   
  sr%ballsize = huge(1.0)
  sr%qv => qv
  sr%nn = nn
  sr%nfound = 0
  sr%centeridx = -1
  sr%correltime = 0
  sr%overflow = PETSC_FALSE !.false.

  sr%results => results

  sr%nalloc = nn   ! will be checked

  sr%ind => tp%ind
  sr%rearrange = tp%rearrange
  if (tp%rearrange) then
    sr%Data => tp%rearranged_data
  else
    sr%Data => tp%the_data
  endif
  sr%dimen = tp%dimen

  call validate_query_storage(nn)
  sr%pq = pq_create(results)

  call search(tp%root)   

  if (tp%sort) then
    call kdtree2_sort_results(nn, results)
  endif
!    deallocate(sr%pqp)
  return
end subroutine kdtree2_n_nearest

! ************************************************************************** !

subroutine validate_query_storage(n)
  
  !
  ! make sure we have enough storage for n
  !

  implicit none
  PetscInt :: n
!    integer, intent(in) :: n

  if (size(sr%results,1) < n) then
      print *, 'KD_TREE_TRANS:  you did not provide enough storage for results(1:n)'
       stop
       return
  endif

  return
end subroutine validate_query_storage

! ************************************************************************** !

recursive subroutine search(node)
  !
  ! This is the innermost core routine of the kd-tree search.  Along
  ! with "process_terminal_node", it is the performance bottleneck.
  !
  ! This version uses a logically complete secondary search of
  ! "box in bounds", whether the sear
  !

  implicit none
  type(Tree_node), pointer :: node
  type(tree_node),pointer :: ncloser, nfarther

  PetscInt :: cut_dim, i
  PetscReal :: qval, dis, ballsize
  PetscReal, pointer :: qv(:)
!    integer                            :: cut_dim, i
    ! ..
!    real(kdkind)                               :: qval, dis
!    real(kdkind)                               :: ballsize
!    real(kdkind), pointer           :: qv(:)
  type(interval), pointer :: box(:)

  if ((associated(node%left) .and. associated(node%right)) .eqv. PETSC_FALSE) then
    ! we are on a terminal node
    if (sr%nn == 0) then
      call process_terminal_node_fixedball(node)
    else
      call process_terminal_node(node)
    endif
  else
    ! we are not on a terminal node
    qv => sr%qv(1:)
    cut_dim = node%cut_dim
    qval = qv(cut_dim)

    if (qval < node%cut_val) then
      ncloser => node%left
      nfarther => node%right
      dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
    else
      ncloser => node%right
      nfarther => node%left
      dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
    endif

    if (associated(ncloser)) call search(ncloser)

    ! we may need to search the second node.
    if (associated(nfarther)) then
      ballsize = sr%ballsize
!          dis=extra**2
      if (dis <= ballsize) then
             
        ! we do this separately as going on the first cut dimen is often
        ! a good idea.
        ! note that if extra**2 < sr%ballsize, then the next
        ! check will also be false.             
        box => node%box(1:)
        do i=1,sr%dimen
          if (i /= cut_dim) then
            dis = dis + dis2_from_bnd(qv(i),box(i)%lower,box(i)%upper)
            if (dis > ballsize) then
              return
            endif
          endif
        end do

             !
             ! if we are still here then we need to search mroe.
             !
        call search(nfarther)
      endif
    endif
  end if
end subroutine search

! ************************************************************************** !

function dis2_from_bnd(x,amin,amax) result (res)

  implicit none
  
  PetscReal :: x, amin, amax
  PetscReal :: res
  
  if (x > amax) then
    res = (x - amax)**2;
    return
  else
    if (x < amin) then
      res = (amin - x)**2;
      return
    else
      res = 0.0
      return
    endif
  endif
  return
end function dis2_from_bnd

! ************************************************************************** !

subroutine process_terminal_node(node)
    
  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure.

  implicit none
  
  type (tree_node), pointer :: node
  !
  PetscReal, pointer :: qv(:), data(:,:)
  PetscInt, pointer :: ind(:)
!    real(kdkind), pointer          :: qv(:)
 !   integer, pointer       :: ind(:)
!    real(kdkind), pointer          :: data(:,:)
    !
!    integer                :: dimen, i, indexofi, k, centeridx, correltime
!    real(kdkind)                   :: ballsize, sd, newpri
!    logical                :: rearrange

  PetscInt :: dimen, i , indexofi, k, centeridx, correltime
  PetscReal :: ballsize, sd, newpri
  PetscBool :: rearrange
  type(pq), pointer      :: pqp
    
  ! copy values from sr to local variables
  ! Notice, making local pointers with an EXPLICIT lower bound
  ! seems to generate faster code.
  ! why?  I don't know.
  qv => sr%qv(1:)
  pqp => sr%pq
  dimen = sr%dimen
  ballsize = sr%ballsize
  rearrange = sr%rearrange
  ind => sr%ind(1:)
  data => sr%Data(1:,1:)
  centeridx = sr%centeridx
  correltime = sr%correltime

  !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window?
  !    include_point = .true.    ! by default include all points
  ! search through terminal bucket.

  mainloop: do i = node%l, node%u
              if (rearrange) then
                sd = 0.0
                do k = 1,dimen
                  sd = sd + (data(k,i) - qv(k))**2
                  if (sd>ballsize) cycle mainloop
                end do
                indexofi = ind(i)  ! only read it if we have not broken out
              else
                indexofi = ind(i)
                sd = 0.0
                do k = 1,dimen
                  sd = sd + (data(k,indexofi) - qv(k))**2
                  if (sd > ballsize) cycle mainloop
                end do
              endif

              if (centeridx > 0) then ! doing correlation interval?
                if (abs(indexofi-centeridx) < correltime) cycle mainloop
              endif


       
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

              if (sr%nfound < sr%nn) then
          
                ! add this point unconditionally to fill list.
          
                sr%nfound = sr%nfound + 1
                newpri = pq_insert(pqp,sd,indexofi)
                if (sr%nfound == sr%nn) ballsize = newpri
                ! we have just filled the working list.
                ! put the best square distance to the maximum value
                ! on the list, which is extractable from the PQ.
              else
          
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list.
          ! Hence we replace that with the current one.
          
                ballsize = pq_replace_max(pqp,sd,indexofi)
              endif
            end do mainloop
    
  ! Reset sr variables which may have changed during loop
    
  sr%ballsize = ballsize

end subroutine process_terminal_node

! ************************************************************************** !

subroutine process_terminal_node_fixedball(node)
    
  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure, i.e.
  ! save all within a fixed ball.

  implicit none
  
  type(tree_node), pointer          :: node

  PetscReal, pointer :: qv(:), data(:,:)
  PetscInt, pointer :: ind(:)
  
 !   real(kdkind), pointer          :: qv(:)
 !   integer, pointer       :: ind(:)
 !   real(kdkind), pointer          :: data(:,:)
  !

  PetscInt :: nfound, dimen, i, indexofi, k
  PetscInt :: centeridx, correltime, nn
  PetscReal :: ballsize, sd
  PetscBool :: rearrange
!    integer                :: nfound
!    integer                :: dimen, i, indexofi, k
!    integer                :: centeridx, correltime, nn
!    real(kdkind)                   :: ballsize, sd
!    logical                :: rearrange

    
  ! copy values from sr to local variables
    
  qv => sr%qv(1:)
  dimen = sr%dimen
  ballsize = sr%ballsize
  rearrange = sr%rearrange
  ind => sr%ind(1:)
  data => sr%Data(1:,1:)
  centeridx = sr%centeridx
  correltime = sr%correltime
  nn = sr%nn ! number to search for
  nfound = sr%nfound

  ! search through terminal bucket.
  mainloop: do i = node%l, node%u

       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       ! which index to the point do we use?

            if (rearrange) then
              sd = 0.0
              do k = 1,dimen
                sd = sd + (data(k,i) - qv(k))**2
                if (sd > ballsize) cycle mainloop
              end do
              indexofi = ind(i)  ! only read it if we have not broken out
            else
              indexofi = ind(i)
              sd = 0.0
              do k = 1,dimen
                sd = sd + (data(k,indexofi) - qv(k))**2
                if (sd>ballsize) cycle mainloop
              end do
            endif

            if (centeridx > 0) then ! doing correlation interval?
              if (abs(indexofi - centeridx) < correltime) cycle mainloop
            endif

            nfound = nfound+1
            if (nfound > sr%nalloc) then
            ! oh nuts, we have to add another one to the tree but
            ! there isn't enough room.
              sr%overflow = PETSC_TRUE!.true.
            else
              sr%results(nfound)%dis = sd
              sr%results(nfound)%idx = indexofi
            endif
          end do mainloop
    
  ! Reset sr variables which may have changed during loop   
  sr%nfound = nfound
end subroutine process_terminal_node_fixedball

! ************************************************************************** !

subroutine kdtree2_sort_results(nfound,results)

  !  Use after search to sort results(1:nfound) in order of increasing
  !  distance.

  implicit none

  PetscInt :: nfound
!  integer, intent(in)          :: nfound
  type(kdtree2_result), target :: results(:)

  if (nfound > 1) call heapsort_struct(results,nfound)

  return
end subroutine kdtree2_sort_results

! ************************************************************************** !

subroutine heapsort_struct(a,n)
  
  ! Sort a(1:n) in ascending order

  implicit none
               
  PetscInt :: n
  type(kdtree2_result),intent(inout) :: a(:)
  type(kdtree2_result) :: value ! temporary value

  PetscInt :: i, j, ileft, iright

  ileft = n/2 + 1
  iright = n

  if(n == 1) return

    do
      if(ileft > 1)then
        ileft = ileft - 1
        value = a(ileft)
      else
        value = a(iright)
        a(iright) = a(1)
        iright = iright - 1
        if (iright == 1) then
          a(1) = value
          return
        endif
      endif
      i = ileft
      j = 2 * ileft
      do while (j <= iright)
        if(j < iright) then
          if(a(j)%dis < a(j+1)%dis) j = j + 1
        endif
        if(value%dis < a(j)%dis) then
          a(i) = a(j);
          i = j
          j = j + j
        else
          j = iright + 1
        endif
      end do
      a(i) = value
    end do
end subroutine heapsort_struct

end module kdtree2_module
