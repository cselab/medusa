      !-------------------------------------------------------------------------
      !  Include file for checking if two subs are neighbors and adding 
      !  them to the neighbor lists if they are. This uses symmetry.
      !
      !  INPUT:   INTEGER :: ip,jp -- indices of the two subs, including
      !                               the periodic images.
      !           INTEGER :: ii,jj -- indices of the corresponding real
      !                               subs (without images).
      !  OUTPUT:  INTEGER, DIMENSION(:  ) :: nneigh
      !           INTEGER, DIMENSION(:,:) :: ineigh
      !
      !-------------------------------------------------------------------------
      !  $Log: ppm_find_neigh_subs_2d.inc,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.2  2006/03/28 21:00:24  ivos
      !  The EXIT is skipped in the __VECTOR case to allow vectorization.
      !
      !  Revision 1.1  2006/03/28 19:26:35  ivos
      !  Initial implementation. Tested for 2d/3d/freespace/periodic.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      mx1 = MIN(max_sub(1,ip),max_sub(1,jp))
      mx2 = MIN(max_sub(2,ip),max_sub(2,jp))

      mn1 = MAX(min_sub(1,ip),min_sub(1,jp))
      mn2 = MAX(min_sub(2,ip),min_sub(2,jp))

      !-------------------------------------------------------------------------
      !  Check for adjacent faces. SHARK: The funny logic is faster.
      !-------------------------------------------------------------------------
      IF (.NOT.(mx1.LT.min_sub(1,ip).OR.mx1.GT.max_sub(1,ip)).AND.    &
     &    .NOT.(mx1.LT.min_sub(1,jp).OR.mx1.GT.max_sub(1,jp)).AND.    &
     &    .NOT.(mx2.LT.min_sub(2,ip).OR.mx2.GT.max_sub(2,ip)).AND.    &
     &    .NOT.(mx2.LT.min_sub(2,jp).OR.mx2.GT.max_sub(2,jp))) THEN

          !---------------------------------------------------------------------
          !  Check that this sub is not already in the list and exclude
          !  the tested sub itself. This check only needs to be done
          !  near periodic boundaries if one of the cells is a ghost cell.
          !  Only in this configuration do we risk double entries.
          !---------------------------------------------------------------------
          isin = .FALSE.
          IF ((pbdrx.OR.pbdry).AND.(ibox.NE.jbox)) THEN
              IF (ii .EQ. jj) THEN
                  ! do not add the periodic image of myself
                  isin = .TRUE.
              ELSE
                  DO kk=1,nneigh(ii)
                      IF (ineigh(kk,ii) .EQ. jj) THEN
                          isin = .TRUE.
#ifndef __VECTOR
                          EXIT
#endif
                      ENDIF
                  ENDDO
              ENDIF
          ENDIF
          IF (.NOT. isin) THEN
              !-----------------------------------------------------------------
              !  Check if we have enough memory allocated
              !-----------------------------------------------------------------
              IF (nneigh(ii).GE.isize.OR.nneigh(jj).GE.isize) THEN
                  !-------------------------------------------------------------
                  !  if not, reallocate and update the isize variable
                  !-------------------------------------------------------------
                  iopt   = ppm_param_alloc_grow_preserve
                  ldc(1) = isize + 100
                  ldc(2) = nsubsplus
                  CALL ppm_alloc(ineigh,ldc,iopt,info)
                  IF (info.NE.0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_find_neigh',  &
     &                    'allocation of ineigh failed',__LINE__,info)
                      GOTO 9999
                  ENDIF
                  !-------------------------------------------------------------
                  !  Save the size 
                  !-------------------------------------------------------------
                  isize  = ldc(1)
              ENDIF 
      
              !-----------------------------------------------------------------
              !  Store the neighbour set
              !-----------------------------------------------------------------
              nneigh(ii)            = nneigh(ii) + 1
              ineigh(nneigh(ii),ii) = jj
              nneigh(jj)            = nneigh(jj) + 1
              ineigh(nneigh(jj),jj) = ii
          ENDIF 
      ENDIF 

