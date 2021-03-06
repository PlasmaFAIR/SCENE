Caveat receptor.  (Jack) dongarra@anl-mcs, (Eric Grosse) research!ehg
Compliments of netlib   Sun Jul  6 09:34:18 CDT 1986
C
C***********************************************************************
C
C     File of the  DOUBLE PRECISION  Level 2 BLAS routines:
C
C      DGEMV, DGBMV, DSYMV, DSBMV, DSPMV, DTRMV, DTBMV, DTPMV,
C      DGER , DSYR , DSPR ,
C      DSYR2, DSPR2,
C      DTRSV, DTBSV, DTPSV.
C
C     See:
C
C        Dongarra J. J., Du Croz J. J., Hammarling S. and Hanson R. J..
C        A proposal for an extended set of Fortran Basic Linear Algebra
C        Subprograms. Technical Memorandum No.41 (revision 1),
C        Mathematics and Computer Science Division, Argone National
C        Laboratory, 9700 South Cass Avenue, Argonne, Illinois 60439,
C        USA or NAG Technical Report TR4/85, Nuemrical Algorithms Group
C        Inc., 1101 31st Street, Suite 100, Downers Grove, Illinois
C        60515-1263, USA.
C
C***********************************************************************
C
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
      CHARACTER*1        TRANS
      INTEGER            M, N, LDA, INCX, INCY
      DOUBLE PRECISION   ALPHA, A( LDA, * ), X( * ), BETA, Y( * )
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y
*.
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           m.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Note that TRANS, M, N and LDA must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*     OK = ( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ).OR.
*    $       ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
*    $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )
*    $     .AND.
*    $     ( M.GE.0 )
*    $     .AND.
*    $     ( N.GE.0 )
*    $     .AND.
*    $     ( LDA.GE.M )
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
*
      INTEGER            I     , IX    , IY    , J     , JX    , JY
      INTEGER            KX    , KY    , LENX  , LENY
      DOUBLE PRECISION   ONE   ,         ZERO
      PARAMETER        ( ONE   = 1.0D+0, ZERO  = 0.0D+0 )
      DOUBLE PRECISION   TEMP
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.
     $    ( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
         IF( BETA.NE.ONE )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         END IF
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( LENX - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( LENY - 1 )*INCY
         END IF
         IF( BETA.NE.ONE )THEN
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP
  100       CONTINUE
         ELSE
            JY = KY
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP  = TEMP + A( I, J )*X( IX )
                  IX    = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGEMV .
*
      END
