!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!


      Module io

        Integer, Save :: iu_log = 6

        Public :: output_dir, amr_log_file
        Character (Len=80) :: output_dir
        Character (Len=80) :: amr_log_file

      End Module io
