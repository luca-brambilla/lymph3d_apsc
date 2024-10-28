module MOD_MPI_CUSTOM
! Do not forget to "register" the new MPI reductions at the beginning of SPEED.f90
! via (for example)
! call MPI_OP_CREATE(MPI_ZERO_OVERWRITE, .TRUE., MPI_ZERO_OVERWRITE_OP, mpi_user_reduction_error)


      integer*4 :: mpi_user_reduction_error
      integer*4 :: MPI_ZERO_OVERWRITE_OP
      integer*4 :: MPI_OVERWRITE_BY_NEW_OP


contains

      function MPI_ZERO_OVERWRITE(invec, inoutvec, len, type)
      real*8 :: invec(0:len-1), inoutvec(0:len-1)
      integer*4 :: len, type, i

      do i=0,len-1
      if (ABS(invec(i)) .lt. 1E-16) then
            !inoutvec(i)=inoutvec(i)
      elseif (ABS(inoutvec(i)) .lt. 1E-16) then
            inoutvec(i)=invec(i)
      else
            inoutvec(i)=1.d0/2.d0*(invec(i)+inoutvec(i))
      endif
      enddo
      end function MPI_ZERO_OVERWRITE

      function MPI_OVERWRITE_BY_NEW(invec, inoutvec, len ,type)

      real*8 :: invec(0:len-1), inoutvec(0:len-1)
      integer*4 :: len, type, i

      do i=0,len-1

      if (ABS(invec(i)) .lt. 1E-16) then
            ! Local vector has no entry --> Leave collective vector as it is (probably has some other value)
      else
            ! If the local vector has some entry, accept that over the existing value
            ! (they should coincide anyways)
            inoutvec(i)=invec(i)
      endif
      enddo

      end function MPI_OVERWRITE_BY_NEW

end module MOD_MPI_CUSTOM