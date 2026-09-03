module constval_module
  ! -- modules
  use kind_module, only: I4, SP, DP

  implicit none
  integer(I4), parameter :: VARLEN = 30, CHALEN = 100, TIMELEN = 3
  integer(I4), parameter :: FACE = 6

  integer(I4), parameter :: INF_SPEC = 0, INF_CLAS = 1, INF_POIN = 2
  integer(I4), parameter :: INF_2DTX = 3, INF_2DBI = 4, INF_3DTX = 5
  integer(I4), parameter :: INF_3DBI = 6, INF_EXTR = 7

  integer(I4), parameter :: OUTF_TABL = 1, OUTF_2DBI = 2, OUTF_3DBI = 3
  integer(I4), parameter :: INOVAL = -9999

  real(SP), parameter :: MINSEC = 60_SP, HOURSEC = 3600_SP
  real(SP), parameter :: DAYSEC = 86400_SP, YEARSEC = 31536000_SP

  real(SP), parameter :: SZERO = 0.00_SP, SONE = 1.00_SP, SHALF = 0.50_SP
  real(SP), parameter :: STWO = 2.00_SP, SQUA = 0.25_SP, SINFI = 1.00E+10_SP
  real(SP), parameter :: SSMAL = 1.00E-10_SP, SNOVAL = -9999_SP

  real(DP), parameter :: DZERO = 0.00_DP, DONE = 1.00_DP, DHALF = 0.50_DP
  real(DP), parameter :: DTWO = 2.00_DP, DQUA = 0.25_DP, DINFI = 1.00E+10_DP
  real(DP), parameter :: DSMAL = 1.00E-10_DP, DNOVAL = -9999_DP

  real(DP), parameter :: NEWPER_BASE = epsilon(1.00_SP)
  real(DP), parameter :: UROUND = epsilon(1.00_DP)

  character(7), parameter :: OUTFORM = "(E13.6)"

end module constval_module
