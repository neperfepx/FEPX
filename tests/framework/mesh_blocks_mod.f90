! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Mesh distribution for test programs on any number of ranks (ROADMAP step 9.17, D108): the mesh
!> is split in equal contiguous blocks, by writing them into its partition table before
!> mesh_distribute. FEPX itself runs a mesh only on as many ranks as it has partitions; the unit
!> tests use one mesh on 1, 2, 3 and 4 ranks.
!>
!> Source : FEPX. The blocks port v2 parallel_mod.f90::mre_decomp1d (tag v2-final-comparison).
!>          Contains no numerical method.
!>
!> Public entities:
!>   even_block                 first and last global number of a rank's block (mre_decomp1d)
!>   mesh_load_blocks           read a mesh and distribute it in blocks (collective)
!>   mesh_load_table_or_blocks  the mesh's table if it has as many partitions as ranks, else blocks
module mesh_blocks_mod

  use kinds_mod, only: ik
  use parallel_mod, only: num_ranks, my_rank, root_rank
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_distribute_mod, only: mesh_distribute

  implicit none

  private

  public :: even_block, mesh_load_blocks, mesh_load_table_or_blocks

contains

  !> First and last global number of block `part` (0-based) of n items in num_parts blocks; an
  !> empty block has last = first - 1.
  !>
  !> Arguments:
  !>   n          (in)  integer(ik)  number of items
  !>   num_parts  (in)  integer(ik)  number of blocks
  !>   part       (in)  integer(ik)  block (0-based)
  !>   bounds     (out) integer(ik)(2)  first and last number (1-based)
  !> Locals:
  !>   num_local  integer(ik)  items in the block
  !>   deficit    integer(ik)  number of blocks with one more item
  pure function even_block(n, num_parts, part) result(bounds)

    integer(ik), intent(in) :: n
    integer(ik), intent(in) :: num_parts
    integer(ik), intent(in) :: part
    integer(ik) :: bounds(2)

    integer(ik) :: num_local
    integer(ik) :: deficit

    ! v2: base size n / num_parts; the first `deficit` blocks get one more item.
    num_local = n/num_parts
    deficit = mod(n, num_parts)
    bounds(1) = part*num_local + 1_ik + min(part, deficit)
    if (part < deficit) num_local = num_local + 1_ik
    bounds(2) = bounds(1) + num_local - 1_ik
    ! v2: the last block ends at n.
    if (bounds(2) > n .or. part == num_parts - 1_ik) bounds(2) = n

  end function even_block

  !> Read the mesh file on root_rank, replace its partition table by blocks, distribute.
  !>
  !> Arguments:
  !>   mesh       (out)           type(mesh_t)     this rank's mesh
  !>   mesh_name  (in)            character        mesh file name (used on root_rank)
  !>   aux_names  (in, optional)  character(*)(:)  auxiliary file names (root_rank)
  !> Locals:
  !>   file_mesh  type(mesh_file_t)  whole mesh (root_rank only)
  subroutine mesh_load_blocks(mesh, mesh_name, aux_names)

    type(mesh_t), intent(out) :: mesh
    character(len=*), intent(in) :: mesh_name
    character(len=*), intent(in), optional :: aux_names(:)

    type(mesh_file_t) :: file_mesh

    if (my_rank == root_rank) then
      call mesh_file_load(file_mesh, mesh_name, aux_names)
      call set_block_table(file_mesh)
    end if
    call mesh_distribute(file_mesh, mesh)

  end subroutine mesh_load_blocks

  !> Read the mesh file on root_rank and distribute it by its partition table if it has as many
  !> partitions as there are ranks, else in blocks.
  !>
  !> Arguments:
  !>   mesh       (out)           type(mesh_t)     this rank's mesh
  !>   mesh_name  (in)            character        mesh file name (used on root_rank)
  !>   aux_names  (in, optional)  character(*)(:)  auxiliary file names (root_rank)
  !> Locals:
  !>   file_mesh  type(mesh_file_t)  whole mesh (root_rank only)
  subroutine mesh_load_table_or_blocks(mesh, mesh_name, aux_names)

    type(mesh_t), intent(out) :: mesh
    character(len=*), intent(in) :: mesh_name
    character(len=*), intent(in), optional :: aux_names(:)

    type(mesh_file_t) :: file_mesh

    if (my_rank == root_rank) then
      call mesh_file_load(file_mesh, mesh_name, aux_names)
      if (file_mesh%num_parts /= num_ranks) call set_block_table(file_mesh)
    end if
    call mesh_distribute(file_mesh, mesh)

  end subroutine mesh_load_table_or_blocks

  !> Write the blocks into a parsed mesh's partition table.
  !>
  !> Arguments:
  !>   file_mesh  (inout)  type(mesh_file_t)  whole mesh
  !> Locals:
  !>   i_rank  integer(ik)     block index
  !>   bounds  integer(ik)(2)  block bounds
  subroutine set_block_table(file_mesh)

    type(mesh_file_t), intent(inout) :: file_mesh

    integer(ik) :: i_rank
    integer(ik) :: bounds(2)

    ! An unpartitioned mesh has no node table yet.
    if (allocated(file_mesh%node_part)) deallocate (file_mesh%node_part)
    allocate (file_mesh%node_part(file_mesh%num_nodes))
    do i_rank = 0_ik, num_ranks - 1_ik
      bounds = even_block(file_mesh%num_elts, num_ranks, i_rank)
      file_mesh%elt_part(bounds(1):bounds(2)) = i_rank + 1_ik
      bounds = even_block(file_mesh%num_nodes, num_ranks, i_rank)
      file_mesh%node_part(bounds(1):bounds(2)) = i_rank + 1_ik
    end do
    file_mesh%num_parts = num_ranks

  end subroutine set_block_table

end module mesh_blocks_mod
