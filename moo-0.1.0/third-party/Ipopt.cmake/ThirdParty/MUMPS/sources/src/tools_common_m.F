C
C  This file is part of MUMPS 5.8.1, released
C  on Wed Jul 30 16:49:18 UTC 2025
C
C
C  Copyright 1991-2025 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      MODULE MUMPS_TOOLS_COMMON_M
C
C     This module with arithmetic-independant utility routines
C     such as integer sorting algorithms
C
      PRIVATE
      PUBLIC :: MUMPS_MERGESORT, MUMPS_MERGESWAP, MUMPS_MERGESWAP1
      CONTAINS
      SUBROUTINE MUMPS_MERGESORT(N, K, L)
      IMPLICIT NONE
C     Plain implementation of the merge-sort algorithm
C     as described in:
C
C     D. E. Knuth "The Art of Computer Programming,"
C     vol.3: Sorting and Searching, Addison-Wesley, 1973
      INTEGER    :: N
      INTEGER    :: K(:), L(0:)
      INTEGER    :: P, Q, S, T
      L(0) = 1
      T = N + 1
      DO  P = 1,N - 1
         IF (K(P) <= K(P+1)) THEN
            L(P) = P + 1
         ELSE
            L(T) = - (P+1)
            T = P
       END IF
      END DO
      L(T) = 0
      L(N) = 0
      IF (L(N+1) == 0) THEN
         RETURN 
      ELSE
         L(N+1) = iabs(L(N+1))
      END IF
 200  CONTINUE
      S = 0
      T = N+1
      P = L(S)
      Q = L(T)
      IF(Q .EQ. 0) RETURN
 300  CONTINUE
      IF(K(P) .GT. K(Q)) GOTO 600 
      L(S) = sign(P,L(S))
      S = P
      P = L(P)
      IF (P .GT. 0) GOTO 300
      L(S) = Q
      S = T
      DO
         T = Q
         Q = L(Q)
         IF (Q .LE. 0) EXIT
      END DO
      GOTO 800
 600  CONTINUE
      L(S) = sign(Q, L(S))
      S = Q
      Q = L(Q)
      IF (Q .GT. 0) GOTO 300
      L(S) = P
      S = T
      DO
         T = P
         P = L(P)
         IF (P .LE. 0) EXIT
      END DO
 800  CONTINUE
      P = -P
      Q = -Q
      IF(Q.EQ.0) THEN
         L(S) = sign(P, L(S))
         L(T) = 0
         GOTO 200
      END IF
      GOTO 300
      END SUBROUTINE MUMPS_MERGESORT
      SUBROUTINE MUMPS_MERGESWAP1(N, L, A)
      IMPLICIT NONE
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A(LP)
         A(LP)   = A(I)
         A(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP1
      SUBROUTINE MUMPS_MERGESWAP(N, L, A1, A2)
      IMPLICIT NONE
      INTEGER   :: I, LP, ISWAP, N
      INTEGER   :: L(0:), A1(:), A2(:)
      LP = L(0)
      I  = 1
      DO 
         IF ((LP==0).OR.(I>N)) EXIT
         DO 
            IF (LP >= I) EXIT
            LP = L(LP)
         END DO
         ISWAP    = A1(LP)
         A1(LP)   = A1(I)
         A1(I)    = ISWAP
         ISWAP    = A2(LP)
         A2(LP)   = A2(I)
         A2(I)    = ISWAP
         ISWAP    = L(LP)
         L(LP) = L(I)
         L(I)  = LP
         LP = ISWAP 
         I  = I + 1
      ENDDO
      END SUBROUTINE MUMPS_MERGESWAP
      END MODULE MUMPS_TOOLS_COMMON_M
