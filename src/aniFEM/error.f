C ================================================================
      Subroutine errMesFEM(iERR, routine, message)
C ================================================================
C Warnings:        1000 - the quality has been not reached
C
C Normal errors:   1000 - the quality has been not reached
C(memory errors)   1001 - not enough memory for Integer arrays
C                  1002 - not enough memory for Real*8 arrays
C                  1003 - local parameter MaxP is small
C                  1004 - local parameter MaxF is small
C                  1006 - local parameter MaxE is small
C                  1007 - local parameter MaxS is small
C                  1012 - local parameter MaxJP is small
C                  1013 - local parameter MaxIA is small
C                  1014 - local parameter MaxA  is small
C                  1101 - one of the local mesh parameters 
C                        (MaxP, MaxF or MaxE) is small
C
C User errors:     2001 : 2003 - errors due to probably incorrect
C                                user routines
C                  2011 : 2011 - incorrect or out of bounds input data
C
C Library errors:  3001 : 3002 - errors in routine DSORT
C                  3011 : 3012 - errors in routine DSYEV
C                  3021        - errors in one of the AMG routines
C
C I/O errors:      4001 : 4002 - errors in input files
C                  4101 : 4103 - errors in input data
C                  4201        - output chanel number is wrong
C
C Internal errors: 5001 : 5022 - errors in mesh checking (1st level)
C                  5101 : 5102 - errors in list updating
C                  5201 : 5201 - errors in mesh checking (2nd level)
C                  6001 : 6005 - system errors
C                  6101 : 6102 - errorneous input for LINTRP3D 
C                  6103 : 6105 - errors in LINTRP3D algorithm
C                  7001 : 7002 - errors in applications working with
C                                curvilinear boundaries
C
C Debug errors:    8001 : 8030 - debug errors in ani2_test.f
C
C ================================================================
      Integer iERR
      Character*(*) routine, message
C ================================================================

      Write(*,5000) iERR, routine, message
      Stop

 5000 Format(/,'Error:', I7,/, 'Routine: ', A,/, 'Message: ', A,/)
      End


