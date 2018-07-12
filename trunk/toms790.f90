SUBROUTINE CSHEP2 (N,X,Y,F,NC,NW,NR, LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A,IER)

!*****************************************************************************80
!
!! CSHEP2 defines a C2 function F(X,Y) interpolating scattered data.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/13/97
!
!   This subroutine computes a set of parameters defining a
! C2 (twice continuously differentiable) bivariate function
! C(X,Y) which interpolates data values F at a set of N
! arbitrarily distributed points (X,Y) in the plane (nodes).
! The interpolant C may be evaluated at an arbitrary point
! by function CS2VAL, and its first partial derivatives are
! computed by Subroutine CS2GRD.
!
!   The interpolation scheme is a modified Cubic Shepard
! method:
!
! C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]
!
! for bivariate functions W(k) and C(k).  The nodal func-
! tions are given by
!
!  C(k)(x,y) = A(1,k)*(x-X(k))**3 +
!              A(2,k)*(x-X(k))**2*(y-Y(k)) +
!              A(3,k)*(x-X(k))*(y-Y(k))**2 +
!              A(4,k)*(y-Y(k))**3 + A(5,k)*(x-X(k))**2 +
!              A(6,k)*(x-X(k))*(y-Y(k)) + A(7,k)*(y-Y(k))**2
!              + A(8,k)*(x-X(k)) + A(9,k)*(y-Y(k)) + F(k) .
!
! Thus, C(k) is a cubic function which interpolates the data
! value at node k.  Its coefficients A(,k) are obtained by a
! weighted least squares fit to the closest NC data points
! with weights similar to W(k).  Note that the radius of
! influence for the least squares fit is fixed for each k,
! but varies with k.
!
! The weights are taken to be
!
!   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,
!
! where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is
! the Euclidean distance between (x,y) and (X(k),Y(k)).  The
! radius of influence R(k) varies with k and is chosen so
! that NW nodes are within the radius.  Note that W(k) is
! not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
! as (x,y) approaches (X(k),Y(k)).
!
! On input:
!
!       N = Number of nodes and data values.  N .GE. 10.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       F = Array of length N containing the data values
!           in one-to-one correspondence with the nodes.
!
!       NC = Number of data points to be used in the least
!            squares fit for coefficients defining the nodal
!            functions C(k).  Values found to be optimal for
!            test data sets ranged from 11 to 25.  A recom-
!            mended value for general data sets is NC = 17.
!            For nodes lying on (or close to) a rectangular
!            grid, the recommended value is NC = 11.  In any
!            case, NC must be in the range 9 to Min(40,N-1).
!
!       NW = Number of nodes within (and defining) the radii
!            of influence R(k) which enter into the weights
!            W(k).  For N sufficiently large, a recommended
!            value is NW = 30.  In general, NW should be
!            about 1.5*NC.  1 .LE. NW .LE. Min(40,N-1).
!
!       NR = Number of rows and columns in the cell grid de-
!            fined in Subroutine STORE2.  A rectangle con-
!            taining the nodes is partitioned into cells in
!            order to increase search efficiency.  NR =
!            Sqrt(N/3) is recommended.  NR .GE. 1.
!
! The above parameters are not altered by this routine.
!
!       LCELL = Array of length .GE. NR**2.
!
!       LNEXT = Array of length .GE. N.
!
!       RW = Array of length .GE. N.
!
!       A = Array of length .GE. 9N.
!
! On output:
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k).
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Note that the output parameters described above are not
! defined unless IER = 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NC, NW, or NR is outside its
!                     valid range.
!             IER = 2 if duplicate nodes were encountered.
!             IER = 3 if all nodes are collinear.
!
! Modules required by CSHEP2:  GETNP2, GIVENS, ROTATE,
!                                SETUP2, STORE2
!
! Intrinsic functions called by CSHEP2:  ABS, DBLE, MAX,
!                                          MIN, SQRT
!
  implicit none

  INTEGER N, NC, NW, NR, LCELL(NR,NR), LNEXT(N), IER
  DOUBLE PRECISION  X(N), Y(N), F(N), XMIN, YMIN, DX, &
                    DY, RMAX, RW(N), A(9,N)
  INTEGER LMX
  PARAMETER (LMX=40)
  INTEGER I, IERR, IP1, IRM1, IROW, J, JP1, K, LMAX, &
          LNP, NEQ, NN, NNC, NNR, NNW, NP, NPTS(LMX), &
          NCWMAX
  DOUBLE PRECISION B(10,10), C, DDX, DDY, DMIN, DTOL, &
                   FK, RC, RS, RSMX, RSOLD, RTOL, RWS, &
                   S, SF, SFC, SFS, STF, SUM, T, XK, &
                   XMN, YK, YMN
!
  DATA    RTOL/1.D-5/, DTOL/.01/
!
! Local parameters:
!
! B =          Transpose of the augmented regression matrix
! C =          First component of the plane rotation used to
!                zero the lower triangle of B**T -- computed
!                by Subroutine GIVENS
! DDX,DDY =    Local variables for DX and DY
! DMIN =       Minimum of the magnitudes of the diagonal
!                elements of the regression matrix after
!                zeros are introduced below the diagonal
! DTOL =       Tolerance for detecting an ill-conditioned
!                system.  The system is accepted when
!                DMIN*RC .GE. DTOL.
! FK =         Data value at mode K -- F(K)
! I =          Index for A, B, and NPTS
! IERR =       Error flag for the call to Subroutine STORE2
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       Row index for B
! J =          Index for A and B
! JP1 =        J+1
! K =          Nodal function index and column index for A
! LMAX =       Maximum number of NPTS elements
! LMX =        Maximum value of LMAX
! LNP =        Current length of NPTS
! NEQ =        Number of equations in the least squares fit
! NN,NNC,NNR = Local copies of N, NC, and NR
! NNW =        Local copy of NW
! NP =         NPTS element
! NPTS =       Array containing the indexes of a sequence of
!                nodes to be used in the least squares fit
!                or to compute RW.  The nodes are ordered
!                by distance from K, and the last element
!                (usually indexed by LNP) is used only to
!                determine RC, or RW(K) if NW > NC.
! NCWMAX =     Max(NC,NW)
! RC =         Radius of influence which enters into the
!                weights for C(K) (see Subroutine SETUP2)
! RS =         Squared distance between K and NPTS(LNP) --
!                used to compute RC and RW(K)
! RSMX =       Maximum squared RW element encountered
! RSOLD =      Squared distance between K and NPTS(LNP-1) --
!                used to compute a relative change in RS
!                between succeeding NPTS elements
! RTOL =       Tolerance for detecting a sufficiently large
!                relative change in RS.  If the change is
!                not greater than RTOL, the nodes are
!                treated as being the same distance from K
! RWS =        Current squared value of RW(K)
! S =          Second component of the plane rotation deter-
!                mined by subroutine GIVENS
! SF =        Scale factor for the linear terms (columns 8
!               and 9) in the least squares fit -- inverse
!               of the root-mean-square distance between K
!               and the nodes (other than K) in the least
!               squares fit
! SFS =       Scale factor for the quadratic terms (columns
!               5, 6, and 7) in the least squares fit --
!               SF*SF
! SFC =       Scale factor for the cubic terms (first 4
!               columns) in the least squares fit -- SF**3
! STF =        Marquardt stabilization factor used to damp
!                out the first 4 solution components (third
!                partials of the cubic) when the system is
!                ill-conditioned.  As STF increases, the
!                fitting function approaches a quadratic
!                polynomial.
! SUM =        Sum of squared Euclidean distances between
!                node K and the nodes used in the least
!                squares fit (unless additional nodes are
!                added for stability)
! T =          Temporary variable for accumulating a scalar
!                product in the back solve
! XK,YK =      Coordinates of node K -- X(K), Y(K)
! XMN,YMN =    Local variables for XMIN and YMIN
!
  NN = N
  NNC = NC
  NNW = NW
  NNR = NR
  NCWMAX = MAX(NNC,NNW)
  LMAX = MIN(LMX,NN-1)
  IF (NNC .LT. 9  .OR.  NNW .LT. 1  .OR.  NCWMAX .GT. &
      LMAX  .OR.  NNR .LT. 1) GO TO 21
!
! Create the cell data structure, and initialize RSMX.
!
  CALL STORE2 (NN,X,Y,NNR, LCELL,LNEXT,XMN,YMN,DDX,DDY,IERR)
  IF (IERR .NE. 0) GO TO 23
  RSMX = 0.
!
! Outer loop on node K:
!
  DO K = 1,NN

    XK = X(K)
    YK = Y(K)
    FK = F(K)
!
! Mark node K to exclude it from the search for nearest
!   neighbors.
!
    LNEXT(K) = -LNEXT(K)
!
! Initialize for loop on NPTS.
!
    RS = 0.
    SUM = 0.
    RWS = 0.
    RC = 0.
    LNP = 0
!
! Compute NPTS, LNP, RWS, NEQ, RC, and SFS.
!
1   SUM = SUM + RS
      IF (LNP .EQ. LMAX) GO TO 2
      LNP = LNP + 1
      RSOLD = RS
      CALL GETNP2 (XK,YK,X,Y,NNR,LCELL,LNEXT,XMN,YMN, &
                   DDX,DDY, NP,RS)
      IF (RS .EQ. 0.) GO TO 22
      NPTS(LNP) = NP
      IF ( (RS-RSOLD)/RS .LT. RTOL ) GO TO 1
      IF (RWS .EQ. 0.  .AND.  LNP .GT. NNW) RWS = RS
      IF (RC .EQ. 0.  .AND.  LNP .GT. NNC) THEN
!
!   RC = 0 (not yet computed) and LNP > NC.  RC = Sqrt(RS)
!     is sufficiently large to (strictly) include NC nodes.
!     The least squares fit will include NEQ = LNP - 1
!     equations for 9 .LE. NC .LE. NEQ .LT. LMAX .LE. N-1.
!
        NEQ = LNP - 1
        RC = SQRT(RS)
        SFS = DBLE(NEQ)/SUM
      end if
!
!   Bottom of loop -- test for termination.
!
      IF (LNP .GT. NCWMAX) GO TO 3
      GO TO 1
!
! All LMAX nodes are included in NPTS.  RWS and/or RC**2 is
!   (arbitrarily) taken to be 10 percent larger than the
!   distance RS to the last node included.
!
2   IF (RWS .EQ. 0.) RWS = 1.1*RS
    IF (RC .EQ. 0.) THEN
      NEQ = LMAX
      RC = SQRT(1.1*RS)
      SFS = DBLE(NEQ)/SUM
    end if
!
! Store RW(K), update RSMX if necessary, and compute SF
!   and SFC.
!
3   RW(K) = SQRT(RWS)
    IF (RWS .GT. RSMX) RSMX = RWS
    SF = SQRT(SFS)
    SFC = SF*SFS
!
! A Q-R decomposition is used to solve the least squares
!   system.  The transpose of the augmented regression
!   matrix is stored in B with columns (rows of B) defined
!   as follows:  1-4 are the cubic terms, 5-7 are the quad-
!   ratic terms, 8 and 9 are the linear terms, and the last
!   column is the right hand side.
!
! Set up the equations and zero out the lower triangle with
!   Givens rotations.
!
    I = 0
4   continue

      I = I + 1
      NP = NPTS(I)
      IROW = MIN(I,10)
      CALL SETUP2 (XK,YK,FK,X(NP),Y(NP),F(NP),SF,SFS, &
                   SFC,RC, B(1,IROW))
      IF (I .EQ. 1) GO TO 4
      IRM1 = IROW-1
      DO J = 1,IRM1
        JP1 = J + 1
        CALL GIVENS (B(J,J),B(J,IROW),C,S)
        CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,IROW))
      end do
      IF (I .LT. NEQ) GO TO 4
!
! Test the system for ill-conditioning.
!
    DMIN = MIN( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)), &
                ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)), &
                ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)) )
    IF (DMIN*RC .GE. DTOL) GO TO 11
    IF (NEQ .EQ. LMAX) GO TO 7
!
! Increase RC and add another equation to the system to
!   improve the conditioning.  The number of NPTS elements
!   is also increased if necessary.
!
6   RSOLD = RS
    NEQ = NEQ + 1
    IF (NEQ .EQ. LMAX) THEN
      RC = SQRT(1.1*RS)
      GO TO 4
    end if
    IF (NEQ .LT. LNP) THEN
!
!   NEQ < LNP.
!
      NP = NPTS(NEQ+1)
      RS = (X(NP)-XK)**2 + (Y(NP)-YK)**2
      IF ( (RS-RSOLD)/RS .LT. RTOL ) GO TO 6
      RC = SQRT(RS)
      GO TO 4
    end if
!
!   NEQ = LNP.  Add an element to NPTS.
!
    LNP = LNP + 1
    CALL GETNP2 (XK,YK,X,Y,NNR,LCELL,LNEXT,XMN,YMN, &
                 DDX,DDY, NP,RS)
    IF (NP .EQ. 0) GO TO 22
    NPTS(LNP) = NP
    IF ( (RS-RSOLD)/RS .LT. RTOL ) GO TO 6
    RC = SQRT(RS)
    GO TO 4
!
! Stabilize the system by damping third partials -- add
!   multiples of the first four unit vectors to the first
!   four equations.
!
7   STF = 1.0/RC

    DO I = 1,4

      B(I,10) = STF
      IP1 = I + 1
      DO J = IP1,10
        B(J,10) = 0.
      end do
      DO J = I,9
        JP1 = J + 1
        CALL GIVENS (B(J,J),B(J,10),C,S)
        CALL ROTATE (10-J,C,S,B(JP1,J),B(JP1,10))
      end do

    end do
!
! Test the damped system for ill-conditioning.
!
    DMIN = MIN( ABS(B(5,5)),ABS(B(6,6)),ABS(B(7,7)), &
                ABS(B(8,8)),ABS(B(9,9)) )
    IF (DMIN*RC .LT. DTOL) GO TO 23
!
! Solve the 9 by 9 triangular system for the coefficients.
!
   11   continue

    DO I = 9,1,-1
      T = 0.
      IF (I .NE. 9) THEN
        IP1 = I + 1
        DO J = IP1,9
          T = T + B(J,I)*A(J,K)
        end do
      end if
      A(I,K) = (B(10,I)-T)/B(I,I)
    end do
!
! Scale the coefficients to adjust for the column scaling.
!
    DO I = 1,4
      A(I,K) = A(I,K)*SFC
    end do
    A(5,K) = A(5,K)*SFS
    A(6,K) = A(6,K)*SFS
    A(7,K) = A(7,K)*SFS
    A(8,K) = A(8,K)*SF
    A(9,K) = A(9,K)*SF
!
! Unmark K and the elements of NPTS.
!
    LNEXT(K) = -LNEXT(K)
    DO I = 1,LNP
      NP = NPTS(I)
      LNEXT(NP) = -LNEXT(NP)
    end do

  end do
!
! No errors encountered.
!
  XMIN = XMN
  YMIN = YMN
  DX = DDX
  DY = DDY
  RMAX = SQRT(RSMX)
  IER = 0
  RETURN
!
! N, NC, NW, or NR is outside its valid range.
!
   21 IER = 1
  RETURN
!
! Duplicate nodes were encountered by GETNP2.
!
   22 IER = 2
  RETURN
!
! No unique solution due to collinear nodes.
!
   23 XMIN = XMN
  YMIN = YMN
  DX = DDX
  DY = DDY
  IER = 3
  RETURN
END
FUNCTION CS2VAL (PX,PY,N,X,Y,F,NR, LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A)

!*****************************************************************************80
!
!! CS2VAL evaluates an interpolant defined by CSHEP2.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
!
!   This function returns the value C(PX,PY), where C is the
! weighted sum of cubic nodal functions defined in Subrou-
! tine CSHEP2.  CS2GRD may be called to compute a gradient
! of C along with the value, and/or to test for errors.
! CS2HES may be called to compute a value, first partial
! derivatives, and second partial derivatives at a point.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C is to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N .GE. 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this function.  The
! parameters other than PX and PY should be input unaltered
! from their values on output from CSHEP2.  This function
! should not be called if a nonzero error flag was returned
! by CSHEP2.
!
! On output:
!
!       CS2VAL = Function value C(PX,PY) unless N, NR, DX,
!                DY, or RMAX is invalid, in which case no
!                value is returned.
!
! Modules required by CS2VAL:  NONE
!
! Intrinsic functions called by CS2VAL:  INT, SQRT
!
  implicit none

  DOUBLE PRECISION CS2VAL
  INTEGER N, NR, LCELL(NR,NR), LNEXT(N)
  DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN, &
                   DX, DY, RMAX, RW(N), A(9,N)
  INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
  DOUBLE PRECISION D, DELX, DELY, R, SW, SWC, W, XP, YP
!
! Local parameters:
!
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! XP,YP =     Local copies of PX and PY -- coordinates of P
!
  XP = PX
  YP = PY
  IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR. &
      DY .LE. 0.  .OR.  RMAX .LT. 0.) RETURN
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.
!
  IMIN = INT((XP-XMIN-RMAX)/DX) + 1
  IMAX = INT((XP-XMIN+RMAX)/DX) + 1
  IF (IMIN .LT. 1) IMIN = 1
  IF (IMAX .GT. NR) IMAX = NR
  JMIN = INT((YP-YMIN-RMAX)/DY) + 1
  JMAX = INT((YP-YMIN+RMAX)/DY) + 1
  IF (JMIN .LT. 1) JMIN = 1
  IF (JMAX .GT. NR) JMAX = NR
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
  IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 6
!
! Accumulate weight values in SW and weighted nodal function
!   values in SWC.  The weights are W(K) = ((R-D)+/(R*D))**3
!   for R = RW(K) and D = distance between P and node K.
!
  SW = 0.
  SWC = 0.
!
! Outer loop on cells (I,J).
!
  DO 4 J = JMIN,JMAX
    DO 3 I = IMIN,IMAX
      K = LCELL(I,J)
      IF (K .EQ. 0) GO TO 3
!
! Inner loop on nodes K.
!
1     DELX = XP - X(K)
      DELY = YP - Y(K)
      D = SQRT(DELX*DELX + DELY*DELY)
      R = RW(K)
      IF (D .GE. R) GO TO 2
      IF (D .EQ. 0.) GO TO 5
      W = (1.0/D - 1.0/R)**3
      SW = SW + W
      SWC = SWC + W*( ( (A(1,K)*DELX+A(2,K)*DELY+ &
                         A(5,K))*DELX + (A(3,K)*DELY+ &
                         A(6,K))*DELY + A(8,K) )*DELX + &
                      ( (A(4,K)*DELY+A(7,K))*DELY + &
                        A(9,K) )*DELY + F(K) )
!
! Bottom of loop on nodes in cell (I,J).
!
2     KP = K
      K = LNEXT(KP)
      IF (K .NE. KP) GO TO 1
3     CONTINUE
4   CONTINUE
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
  IF (SW .EQ. 0.) GO TO 6
  CS2VAL = SWC/SW
  RETURN
!
! (PX,PY) = (X(K),Y(K)).
!
5 CS2VAL = F(K)
  RETURN
!
! All weights are 0 at P.
!
6 CS2VAL = 0.
  RETURN
END
SUBROUTINE CS2GRD (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A, &
  C,CX,CY,IER)

!*****************************************************************************80
!
!! CS2GRD returns the value and gradients of a function defined by CSHEP2.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
!
!   This subroutine computes the value and gradient at P =
! (PX,PY) of the interpolatory function C defined in Sub-
! routine CSHEP2.  C is a weighted sum of cubic nodal
! functions.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C and its partial derivatives are
!               to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N .GE. 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.
!
! On output:
!
!       C = Value of C at (PX,PY) unless IER .EQ. 1, in
!           which case no values are returned.
!
!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER .EQ. 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C=CX=CY=0).
!
! Modules required by CS2GRD:  None
!
! Intrinsic functions called by CS2GRD:  INT, SQRT
!
  implicit none

  INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
  DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN, &
                   DX, DY, RMAX, RW(N), A(9,N), C, CX, &
                   CY
  INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
  DOUBLE PRECISION CK, CKX, CKY, D, DELX, DELY, R, SW, &
                   SWC, SWCX, SWCY, SWS, SWX, SWY, T, W, &
                   WX, WY, XP, YP
!
! Local parameters:
!
! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! T =         Temporary variable
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! XP,YP =     Local copies of PX and PY -- coordinates of P
!
  XP = PX
  YP = PY
  IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR. &
      DY .LE. 0.  .OR.  RMAX .LT. 0.) GO TO 6
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.
!
  IMIN = INT((XP-XMIN-RMAX)/DX) + 1
  IMAX = INT((XP-XMIN+RMAX)/DX) + 1
  IF (IMIN .LT. 1) IMIN = 1
  IF (IMAX .GT. NR) IMAX = NR
  JMIN = INT((YP-YMIN-RMAX)/DY) + 1
  JMAX = INT((YP-YMIN+RMAX)/DY) + 1
  IF (JMIN .LT. 1) JMIN = 1
  IF (JMAX .GT. NR) JMAX = NR
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
  IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 7
!
! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus
!
!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2
!
!   where SWCX and SWX are partial derivatives with respect
!   to X of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.
!
  SW = 0.
  SWX = 0.
  SWY = 0.
  SWC = 0.
  SWCX = 0.
  SWCY = 0.
!
! Outer loop on cells (I,J).
!
  DO 4 J = JMIN,JMAX
    DO 3 I = IMIN,IMAX
      K = LCELL(I,J)
      IF (K .EQ. 0) GO TO 3
!
! Inner loop on nodes K.
!
1     DELX = XP - X(K)
      DELY = YP - Y(K)
      D = SQRT(DELX*DELX + DELY*DELY)
      R = RW(K)
      IF (D .GE. R) GO TO 2
      IF (D .EQ. 0.) GO TO 5
      T = (1.0/D - 1.0/R)
      W = T**3
      T = -3.0*T*T/(D**3)
      WX = DELX*T
      WY = DELY*T
      T = A(2,K)*DELX + A(3,K)*DELY + A(6,K)
      CKY = ( 3.0*A(4,K)*DELY + A(3,K)*DELX + &
              2.0*A(7,K) )*DELY + T*DELX + A(9,K)
      T = T*DELY + A(8,K)
      CKX = ( 3.0*A(1,K)*DELX + A(2,K)*DELY + &
              2.0*A(5,K) )*DELX + T
      CK = ( (A(1,K)*DELX+A(5,K))*DELX + T )*DELX + &
           ( (A(4,K)*DELY+A(7,K))*DELY + A(9,K) )*DELY + &
           F(K)
      SW = SW + W
      SWX = SWX + WX
      SWY = SWY + WY
      SWC = SWC + W*CK
      SWCX = SWCX + WX*CK + W*CKX
      SWCY = SWCY + WY*CK + W*CKY
!
! Bottom of loop on nodes in cell (I,J).
!
2     KP = K
      K = LNEXT(KP)
      IF (K .NE. KP) GO TO 1
3     CONTINUE
4   CONTINUE
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
  IF (SW .EQ. 0.) GO TO 7
  C = SWC/SW
  SWS = SW*SW
  CX = (SWCX*SW - SWC*SWX)/SWS
  CY = (SWCY*SW - SWC*SWY)/SWS
  IER = 0
  RETURN
!
! (PX,PY) = (X(K),Y(K)).
!
5 C = F(K)
  CX = A(8,K)
  CY = A(9,K)
  IER = 0
  RETURN
!
! Invalid input parameter.
!
6 IER = 1
  RETURN
!
! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D .GE. RW(K) for all K.
!
7 C = 0.
  CX = 0.
  CY = 0.
  IER = 2
  RETURN
END
SUBROUTINE CS2HES (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY,RMAX,RW,A, &
  C,CX,CY,CXX,CXY,CYY,IER)

!*****************************************************************************80
!
!! CS2HES returns value, gradients and Hessian of a function defined by CSHEP2.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
!
!   This subroutine computes the value, gradient, and
! Hessian at P = (PX,PY) of the interpolatory function C
! defined in Subroutine CSHEP2.  C is a weighted sum of
! cubic nodal functions.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P at
!               which C and its partial derivatives are
!               to be evaluated.
!
!       N = Number of nodes and data values defining C.
!           N .GE. 10.
!
!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!       RMAX = Largest element in RW -- maximum radius R(k).
!
!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.
!
!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.
!
!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input
! unaltered from their values on output from CSHEP2.  This
! subroutine should not be called if a nonzero error flag
! was returned by CSHEP2.
!
! On output:
!
!       C = Value of C at (PX,PY) unless IER .EQ. 1, in
!           which case no values are returned.
!
!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER .EQ. 1.
!
!       CXX,CXY,CYY = Second partial derivatives of C at
!                     (PX,PY) unless IER .EQ. 1.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C = 0).
!
! Modules required by CS2HES:  None
!
! Intrinsic functions called by CS2HES:  INT, SQRT
!
  implicit none

  INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
  DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN, &
                   DX, DY, RMAX, RW(N), A(9,N), C, CX, &
                   CY, CXX, CXY, CYY
  INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
  DOUBLE PRECISION CK, CKX, CKXX, CKXY, CKY, CKYY, D, &
                   DELX, DELY, DXSQ, DYSQ, R, SW, SWC, &
                   SWCX, SWCXX, SWCXY, SWCY, SWCYY, SWS, &
                   SWX, SWXX, SWXY, SWY, SWYY, T1, T2, &
                   T3, T4, W, WX, WXX, WXY, WY, WYY, XP, &
                   YP
!
! Local parameters:
!
! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X
!               and Y, respectively
! CKXX,CKXY,CKYY = Second partial derivatives of CK
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! DXSQ,DYSQ = DELX**2, DELY**2
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells
!               intersected by a disk of radius RMAX
!               centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes
!               in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X
!               and Y, respectively
! SWCXX,SWCXY,SWCYY = Second partial derivatives of SWC
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X
!               and Y, respectively
! SWXX,SWXY,SWYY = Second partial derivatives of SW
! T1,T2,T3,T4 = Temporary variables
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X
!               and Y, respectively
! WXX,WXY,WYY = Second partial derivatives of W
! XP,YP =     Local copies of PX and PY -- coordinates of P
!
  XP = PX
  YP = PY
  IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR. &
      DY .LE. 0.  .OR.  RMAX .LT. 0.) GO TO 6
!
! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
!   the range of the search for nodes whose radii include
!   P.  The cells which must be searched are those inter-
!   sected by (or contained in) a circle of radius RMAX
!   centered at P.
!
  IMIN = INT((XP-XMIN-RMAX)/DX) + 1
  IMAX = INT((XP-XMIN+RMAX)/DX) + 1
  IF (IMIN .LT. 1) IMIN = 1
  IF (IMAX .GT. NR) IMAX = NR
  JMIN = INT((YP-YMIN-RMAX)/DY) + 1
  JMAX = INT((YP-YMIN+RMAX)/DY) + 1
  IF (JMIN .LT. 1) JMIN = 1
  IF (JMAX .GT. NR) JMAX = NR
!
! The following is a test for no cells within the circle
!   of radius RMAX.
!
  IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 7
!
! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
!   from K = 1 to N, C(K) is the cubic nodal function value,
!   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
!   ance D(K).  Thus
!
!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2
!
!   where SWCX and SWX are partial derivatives with respect
!   to x of SWC and SW, respectively.  SWCY and SWY are de-
!   fined similarly.  The second partials are
!
!        CXX = ( SW*(SWCXX -    2*SWX*CX) - SWC*SWXX )/SW**2
!        CXY = ( SW*(SWCXY-SWX*CY-SWY*CX) - SWC*SWXY )/SW**2
!        CYY = ( SW*(SWCYY -    2*SWY*CY) - SWC*SWYY )/SW**2
!
!   where SWCXX and SWXX are second partials with respect
!   to x, SWCXY and SWXY are mixed partials, and SWCYY and
!   SWYY are second partials with respect to y.
!
  SW = 0.
  SWX = 0.
  SWY = 0.
  SWXX = 0.
  SWXY = 0.
  SWYY = 0.
  SWC = 0.
  SWCX = 0.
  SWCY = 0.
  SWCXX = 0.
  SWCXY = 0.
  SWCYY = 0.
!
! Outer loop on cells (I,J).
!
  DO 4 J = JMIN,JMAX
    DO 3 I = IMIN,IMAX
      K = LCELL(I,J)
      IF (K .EQ. 0) GO TO 3
!
! Inner loop on nodes K.
!
1     DELX = XP - X(K)
      DELY = YP - Y(K)
      DXSQ = DELX*DELX
      DYSQ = DELY*DELY
      D = SQRT(DXSQ + DYSQ)
      R = RW(K)
      IF (D .GE. R) GO TO 2
      IF (D .EQ. 0.) GO TO 5
      T1 = (1.0/D - 1.0/R)
      W = T1**3
      T2 = -3.0*T1*T1/(D**3)
      WX = DELX*T2
      WY = DELY*T2
      T1 = 3.0*T1*(2.0+3.0*D*T1)/(D**6)
      WXX = T1*DXSQ + T2
      WXY = T1*DELX*DELY
      WYY = T1*DYSQ + T2
      T1 = A(1,K)*DELX + A(2,K)*DELY + A(5,K)
      T2 = T1 + T1 + A(1,K)*DELX
      T3 = A(4,K)*DELY + A(3,K)*DELX + A(7,K)
      T4 = T3 + T3 + A(4,K)*DELY
      CK = (T1*DELX + A(6,K)*DELY + A(8,K))*DELX + &
           (T3*DELY + A(9,K))*DELY + F(K)
      CKX = T2*DELX + (A(3,K)*DELY+A(6,K))*DELY + A(8,K)
      CKY = T4*DELY + (A(2,K)*DELX+A(6,K))*DELX + A(9,K)
      CKXX = T2 + 3.0*A(1,K)*DELX
      CKXY = 2.0*(A(2,K)*DELX + A(3,K)*DELY) + A(6,K)
      CKYY = T4 + 3.0*A(4,K)*DELY
      SW = SW + W
      SWX = SWX + WX
      SWY = SWY + WY
      SWXX = SWXX + WXX
      SWXY = SWXY + WXY
      SWYY = SWYY + WYY
      SWC = SWC + W*CK
      SWCX = SWCX + WX*CK + W*CKX
      SWCY = SWCY + WY*CK + W*CKY
      SWCXX = SWCXX + W*CKXX + 2.0*WX*CKX + CK*WXX
      SWCXY = SWCXY + W*CKXY + WX*CKY + WY*CKX + CK*WXY
      SWCYY = SWCYY + W*CKYY + 2.0*WY*CKY + CK*WYY
!
! Bottom of loop on nodes in cell (I,J).
!
2     KP = K
      K = LNEXT(KP)
      IF (K .NE. KP) GO TO 1
3     CONTINUE
4   CONTINUE
!
! SW = 0 iff P is not within the radius R(K) for any node K.
!
  IF (SW .EQ. 0.) GO TO 7
  C = SWC/SW
  SWS = SW*SW
  CX = (SWCX*SW - SWC*SWX)/SWS
  CY = (SWCY*SW - SWC*SWY)/SWS
  CXX = (SW*(SWCXX-2.0*SWX*CX) - SWC*SWXX)/SWS
  CXY = (SW*(SWCXY-SWY*CX-SWX*CY) - SWC*SWXY)/SWS
  CYY = (SW*(SWCYY-2.0*SWY*CY) - SWC*SWYY)/SWS
  IER = 0
  RETURN
!
! (PX,PY) = (X(K),Y(K)).
!
5 C = F(K)
  CX = A(8,K)
  CY = A(9,K)
  CXX = 2.0*A(5,K)
  CXY = A(6,K)
  CYY = 2.0*A(7,K)
  IER = 0
  RETURN
!
! Invalid input parameter.
!
6 IER = 1
  RETURN
!
! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D .GE. RW(K) for all K.
!
7 C = 0.
  CX = 0.
  CY = 0.
  CXX = 0.
  CXY = 0.
  CYY = 0.
  IER = 2
  RETURN
END
SUBROUTINE GETNP2 (PX,PY,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,DX,DY, NP,DSQ)

!*****************************************************************************80
!
!! GETNP2 gets the closest unmarked node to a specified point.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
!
!   Given a set of N nodes and the data structure defined in
! Subroutine STORE2, this subroutine uses the cell method to
! find the closest unmarked node NP to a specified point P.
! NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
! node is marked if and only if the corresponding LNEXT ele-
! ment is negative.  The absolute values of LNEXT elements,
! however, must be preserved.)  Thus, the closest M nodes to
! P may be determined by a sequence of M calls to this rou-
! tine.  Note that if the nearest neighbor to node K is to
! be determined (PX = X(K) and PY = Y(K)), then K should be
! marked before the call to this routine.
!
!   The search is begun in the cell containing (or closest
! to) P and proceeds outward in rectangular layers until all
! cells which contain points within distance R of P have
! been searched, where R is the distance from P to the first
! unmarked node encountered (infinite if no unmarked nodes
! are present).
!
!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.
!
! On input:
!
!       PX,PY = Cartesian coordinates of the point P whose
!               nearest unmarked neighbor is to be found.
!
!       X,Y = Arrays of length N, for N .GE. 2, containing
!             the Cartesian coordinates of the nodes.
!
!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR .GE. 1.
!
!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.
!
!       LNEXT = Array of length N containing next-node
!               indexes (or their negatives).  Refer to
!               Subroutine STORE2.
!
!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.
!
!   Input parameters other than LNEXT are not altered by
! this routine.  With the exception of (PX,PY) and the signs
! of LNEXT elements, these parameters should be unaltered
! from their values on output from Subroutine STORE2.
!
! On output:
!
!       NP = Index (for X and Y) of the nearest unmarked
!            node to P, or 0 if all nodes are marked or NR
!            .LT. 1 or DX .LE. 0 or DY .LE. 0.  LNEXT(NP)
!            .LT. 0 IF NP .NE. 0.
!
!       DSQ = Squared Euclidean distance between P and node
!             NP, or 0 if NP = 0.
!
! Modules required by GETNP2:  None
!
! Intrinsic functions called by GETNP2:  ABS, INT, SQRT
!
  implicit none

  INTEGER NR, LCELL(NR,NR), LNEXT(*), NP
  DOUBLE PRECISION PX, PY, X(*), Y(*), XMIN, YMIN, DX, &
                   DY, DSQ
  INTEGER I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2, &
          JMAX, JMIN, L, LMIN, LN
  LOGICAL FIRST
  DOUBLE PRECISION DELX, DELY, R, RSMIN, RSQ, XP, YP
!
! Local parameters:
!
! DELX,DELY =   PX-XMIN, PY-YMIN
! FIRST =       Logical variable with value TRUE iff the
!                 first unmarked node has yet to be
!                 encountered
! I,J =         Cell indexes in the range [I1,I2] X [J1,J2]
! I0,J0 =       Indexes of the cell containing or closest
!                 to P
! I1,I2,J1,J2 = Range of cell indexes defining the layer
!                 whose intersection with the range
!                 [IMIN,IMAX] X [JMIN,JMAX] is currently
!                 being searched
! IMIN,IMAX =   Cell row indexes defining the range of the
!                 search
! JMIN,JMAX =   Cell column indexes defining the range of
!                 the search
! L,LN =        Indexes of nodes in cell (I,J)
! LMIN =        Current candidate for NP
! R =           Distance from P to node LMIN
! RSMIN =       Squared distance from P to node LMIN
! RSQ =         Squared distance from P to node L
! XP,YP =       Local copy of PX,PY -- coordinates of P
!
  XP = PX
  YP = PY
!
! Test for invalid input parameters.
!
  IF (NR .LT. 1  .OR.  DX .LE. 0.  .OR.  DY .LE. 0.) &
    GO TO 9
!
! Initialize parameters.
!
  FIRST = .TRUE.
  IMIN = 1
  IMAX = NR
  JMIN = 1
  JMAX = NR
  DELX = XP - XMIN
  DELY = YP - YMIN
  I0 = INT(DELX/DX) + 1
  IF (I0 .LT. 1) I0 = 1
  IF (I0 .GT. NR) I0 = NR
  J0 = INT(DELY/DY) + 1
  IF (J0 .LT. 1) J0 = 1
  IF (J0 .GT. NR) J0 = NR
  I1 = I0
  I2 = I0
  J1 = J0
  J2 = J0
!
! Outer loop on layers, inner loop on layer cells, excluding
!   those outside the range [IMIN,IMAX] X [JMIN,JMAX].
!
1 DO 6 J = J1,J2
    IF (J .GT. JMAX) GO TO 7
    IF (J .LT. JMIN) GO TO 6
    DO 5 I = I1,I2
      IF (I .GT. IMAX) GO TO 6
      IF (I .LT. IMIN) GO TO 5
      IF (J .NE. J1  .AND.  J .NE. J2  .AND.  I .NE. I1 &
          .AND.  I .NE. I2) GO TO 5
!
! Search cell (I,J) for unmarked nodes L.
!
      L = LCELL(I,J)
      IF (L .EQ. 0) GO TO 5
!
!   Loop on nodes in cell (I,J).
!
2     LN = LNEXT(L)
      IF (LN .LT. 0) GO TO 4
!
!   Node L is not marked.
!
      RSQ = (X(L)-XP)**2 + (Y(L)-YP)**2
      IF (.NOT. FIRST) GO TO 3
!
!   Node L is the first unmarked neighbor of P encountered.
!     Initialize LMIN to the current candidate for NP, and
!     RSMIN to the squared distance from P to LMIN.  IMIN,
!     IMAX, JMIN, and JMAX are updated to define the smal-
!     lest rectangle containing a circle of radius R =
!     Sqrt(RSMIN) centered at P, and contained in [1,NR] X
!     [1,NR] (except that, if P is outside the rectangle
!     defined by the nodes, it is possible that IMIN > NR,
!     IMAX < 1, JMIN > NR, or JMAX < 1).  FIRST is reset to
!     FALSE.
!
      LMIN = L
      RSMIN = RSQ
      R = SQRT(RSMIN)
      IMIN = INT((DELX-R)/DX) + 1
      IF (IMIN .LT. 1) IMIN = 1
      IMAX = INT((DELX+R)/DX) + 1
      IF (IMAX .GT. NR) IMAX = NR
      JMIN = INT((DELY-R)/DY) + 1
      IF (JMIN .LT. 1) JMIN = 1
      JMAX = INT((DELY+R)/DY) + 1
      IF (JMAX .GT. NR) JMAX = NR
      FIRST = .FALSE.
      GO TO 4
!
!   Test for node L closer than LMIN to P.
!
3     IF (RSQ .GE. RSMIN) GO TO 4
!
!   Update LMIN and RSMIN.
!
      LMIN = L
      RSMIN = RSQ
!
!   Test for termination of loop on nodes in cell (I,J).
!
4     IF (ABS(LN) .EQ. L) GO TO 5
      L = ABS(LN)
      GO TO 2
5     CONTINUE
6   CONTINUE
!
! Test for termination of loop on cell layers.
!
7 IF (I1 .LE. IMIN  .AND.  I2 .GE. IMAX  .AND. &
      J1 .LE. JMIN  .AND.  J2 .GE. JMAX) GO TO 8
  I1 = I1 - 1
  I2 = I2 + 1
  J1 = J1 - 1
  J2 = J2 + 1
  GO TO 1
!
! Unless no unmarked nodes were encountered, LMIN is the
!   closest unmarked node to P.
!
8 IF (FIRST) GO TO 9
  NP = LMIN
  DSQ = RSMIN
  LNEXT(LMIN) = -LNEXT(LMIN)
  RETURN
!
! Error:  NR, DX, or DY is invalid or all nodes are marked.
!
9 NP = 0
  DSQ = 0.
  RETURN
END
SUBROUTINE GIVENS ( A,B, C,S)

!*****************************************************************************80
!
!! GIVENS constructs a Givens plane rotation.
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!   This subroutine constructs the Givens plane rotation,
!
!           ( C  S)
!       G = (     ) , where C*C + S*S = 1,
!           (-S  C)
!
! which zeros the second component of the vector (A,B)**T
! (transposed).  Subroutine ROTATE may be called to apply
! the transformation to a 2 by N matrix.
!
!   This routine is identical to subroutine SROTG from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).
!
! On input:
!
!       A,B = Components of the vector defining the rota-
!             tion.  These are overwritten by values R
!             and Z (described below) which define C and S.
!
! On output:
!
!       A = Signed Euclidean norm R of the input vector:
!           R = +/-SQRT(A*A + B*B)
!
!       B = Value Z such that:
!             C = SQRT(1-Z*Z) and S=Z if ABS(Z) .LE. 1, and
!             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.
!
!       C = +/-(A/R) or 1 if R = 0.
!
!       S = +/-(B/R) or 0 if R = 0.
!
! Modules required by GIVENS:  None
!
! Intrinsic functions called by GIVENS:  ABS, SQRT
!
  implicit none

  DOUBLE PRECISION A, B, C, S
  DOUBLE PRECISION AA, BB, R, U, V
!
! Local parameters:
!
! AA,BB = Local copies of A and B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   Variables used to scale A and B for computing R
!
  AA = A
  BB = B
  IF (ABS(AA) .LE. ABS(BB)) GO TO 1
!
! ABS(A) > ABS(B).
!
  U = AA + AA
  V = BB/U
  R = SQRT(.25 + V*V) * U
  C = AA/R
  S = V * (C + C)
!
! Note that R has the sign of A, C > 0, and S has
!   SIGN(A)*SIGN(B).
!
  B = S
  A = R
  RETURN
!
! ABS(A) .LE. ABS(B).
!
1 IF (BB .EQ. 0.) GO TO 2
  U = BB + BB
  V = AA/U
!
! Store R in A.
!
  A = SQRT(.25 + V*V) * U
  S = BB/A
  C = V * (S + S)
!
! Note that R has the sign of B, S > 0, and C has
!   SIGN(A)*SIGN(B).
!
  B = 1.
  IF (C .NE. 0.) B = 1./C
  RETURN
!
! A = B = 0.
!
2 C = 1.
  S = 0.
  RETURN
END
SUBROUTINE ROTATE (N,C,S, X,Y )

!*****************************************************************************80
!
!! ROTATE applies a Givens rotation.
!
!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88
!
!                                                ( C  S)
!   This subroutine applies the Givens rotation  (     )  to
!                                                (-S  C)
!                    (X(1) ... X(N))
! the 2 by N matrix  (             ) .
!                    (Y(1) ... Y(N))
!
!   This routine is identical to subroutine SROT from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).
!
! On input:
!
!       N = Number of columns to be rotated.
!
!       C,S = Elements of the Givens rotation.  Refer to
!             subroutine GIVENS.
!
! The above parameters are not altered by this routine.
!
!       X,Y = Arrays of length .GE. N containing the compo-
!             nents of the vectors to be rotated.
!
! On output:
!
!       X,Y = Arrays containing the rotated vectors (not
!             altered if N < 1).
!
! Modules required by ROTATE:  None
!
  implicit none

  INTEGER N
  DOUBLE PRECISION C, S, X(N), Y(N)
  INTEGER I
  DOUBLE PRECISION XI, YI

  DO I = 1,N
    XI = X(I)
    YI = Y(I)
    X(I) = C*XI + S*YI
    Y(I) = -S*XI + C*YI
  end do

  RETURN
END
SUBROUTINE SETUP2 (XK,YK,ZK,XI,YI,ZI,S1,S2,S3,R, ROW)

!*****************************************************************************80
!
!! SETUP2 sets a row of a regression matrix to fit a cubic function to data.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97
!
!   This subroutine sets up the I-th row of an augmented re-
! gression matrix for a weighted least squares fit of a
! cubic function f(x,y) to a set of data values z, where
! f(XK,YK) = ZK.  The first four columns (cubic terms) are
! scaled by S3, the next three columns (quadratic terms)
! are scaled by S2, and the eighth and ninth columns (lin-
! ear terms) are scaled by S1.
!
! On input:
!
!       XK,YK = Coordinates of node K.
!
!       ZK = Data value at node K to be interpolated by f.
!
!       XI,YI,ZI = Coordinates and data value at node I.
!
!       S1,S2,S3 = Scale factors.
!
!       R = Radius of influence about node K defining the
!           weight.
!
! The above parameters are not altered by this routine.
!
!       ROW = Array of length 10.
!
! On output:
!
!       ROW = Array containing a row of the augmented re-
!             gression matrix.
!
! Modules required by SETUP2:  None
!
! Intrinsic function called by SETUP2:  SQRT
!
  implicit none

  DOUBLE PRECISION XK, YK, ZK, XI, YI, ZI, S1, S2, S3, &
                   R, ROW(10)
  INTEGER I
  DOUBLE PRECISION D, DX, DXSQ, DY, DYSQ, W, W1, W2, W3
!
! Local parameters:
!
! D =    Distance between nodes K and I
! DX =   XI - XK
! DXSQ = DX*DX
! DY =   YI - YK
! DYSQ = DY*DY
! I =    DO-loop index
! W =    Weight associated with the row:  (R-D)/(R*D)
!          (0 if D = 0 or D > R)
! W1 =   S1*W
! W2 =   S2*W
! W3 =   W3*W
!
  DX = XI - XK
  DY = YI - YK
  DXSQ = DX*DX
  DYSQ = DY*DY
  D = SQRT(DXSQ + DYSQ)
!
!  Nodes K and I coincide or node I is outside of the radius
!  of influence.  Set ROW to the zero vector.
!
  IF (D .LE. 0.  .OR.  D .GE. R) then
    DO I = 1,10
      ROW(I) = 0.
    end do
    return
  end if

  W = (R-D)/R/D
  W1 = S1*W
  W2 = S2*W
  W3 = S3*W
  ROW(1) = DXSQ*DX*W3
  ROW(2) = DXSQ*DY*W3
  ROW(3) = DX*DYSQ*W3
  ROW(4) = DYSQ*DY*W3
  ROW(5) = DXSQ*W2
  ROW(6) = DX*DY*W2
  ROW(7) = DYSQ*W2
  ROW(8) = DX*W1
  ROW(9) = DY*W1
  ROW(10) = (ZI - ZK)*W

  RETURN
END
SUBROUTINE STORE2 (N,X,Y,NR, LCELL,LNEXT,XMIN,YMIN,DX,DY,IER)

!*****************************************************************************80
!
!! STORE2 creates a closest-point data structure given points in the plane.
!
!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/28/97
!
!   Given a set of N arbitrarily distributed nodes in the
! plane, this subroutine creates a data structure for a
! cell-based method of solving closest-point problems.  The
! smallest rectangle containing the nodes is partitioned
! into an NR by NR uniform grid of cells, and nodes are as-
! sociated with cells.  In particular, the data structure
! stores the indexes of the nodes contained in each cell.
! For a uniform random distribution of nodes, the nearest
! node to an arbitrary point can be determined in constant
! expected time.
!
!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.
!
! On input:
!
!       N = Number of nodes.  N .GE. 2.
!
!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.
!
!       NR = Number of rows and columns in the grid.  The
!            cell density (average number of nodes per cell)
!            is D = N/(NR**2).  A recommended value, based
!            on empirical evidence, is D = 3 -- NR =
!            Sqrt(N/3).  NR .GE. 1.
!
! The above parameters are not altered by this routine.
!
!       LCELL = Array of length .GE. NR**2.
!
!       LNEXT = Array of length .GE. N.
!
! On output:
!
!       LCELL = NR by NR cell array such that LCELL(I,J)
!               contains the index (for X and Y) of the
!               first node (node with smallest index) in
!               cell (I,J), or LCELL(I,J) = 0 if no nodes
!               are contained in the cell.  The upper right
!               corner of cell (I,J) has coordinates (XMIN+
!               I*DX,YMIN+J*DY).  LCELL is not defined if
!               IER .NE. 0.
!
!       LNEXT = Array of next-node indexes such that
!               LNEXT(K) contains the index of the next node
!               in the cell which contains node K, or
!               LNEXT(K) = K if K is the last node in the
!               cell for K = 1,...,N.  (The nodes contained
!               in a cell are ordered by their indexes.)
!               If, for example, cell (I,J) contains nodes
!               2, 3, and 5 (and no others), then LCELL(I,J)
!               = 2, LNEXT(2) = 3, LNEXT(3) = 5, and
!               LNEXT(5) = 5.  LNEXT is not defined if
!               IER .NE. 0.
!
!       XMIN,YMIN = Cartesian coordinates of the lower left
!                   corner of the rectangle defined by the
!                   nodes (smallest nodal coordinates) un-
!                   less IER = 1.  The upper right corner is
!                   (XMAX,YMAX) for XMAX = XMIN + NR*DX and
!                   YMAX = YMIN + NR*DY.
!
!       DX,DY = Dimensions of the cells unless IER = 1.  DX
!               = (XMAX-XMIN)/NR and DY = (YMAX-YMIN)/NR,
!               where XMIN, XMAX, YMIN, and YMAX are the
!               extrema of X and Y.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N < 2 or NR < 1.
!             IER = 2 if DX = 0 or DY = 0.
!
! Modules required by STORE2:  None
!
! Intrinsic functions called by STORE2:  DBLE, INT
!
  implicit none

  INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
  DOUBLE PRECISION X(N), Y(N), XMIN, YMIN, DX, DY
  INTEGER I, J, K, L, NN, NNR
  DOUBLE PRECISION DELX, DELY, XMN, XMX, YMN, YMX
!
! Local parameters:
!
! DELX,DELY = Components of the cell dimensions -- local
!               copies of DX,DY
! I,J =       Cell indexes
! K =         Nodal index
! L =         Index of a node in cell (I,J)
! NN =        Local copy of N
! NNR =       Local copy of NR
! XMN,XMX =   Range of nodal X coordinates
! YMN,YMX =   Range of nodal Y coordinates
!
  NN = N
  NNR = NR

  IF (NN .LT. 2  .OR.  NNR .LT. 1) then
    ier = 1
    return
  end if
!
! Compute the dimensions of the rectangle containing the
!   nodes.
!
  XMN = X(1)
  XMX = XMN
  YMN = Y(1)
  YMX = YMN
  DO K = 2,NN
    IF (X(K) .LT. XMN) XMN = X(K)
    IF (X(K) .GT. XMX) XMX = X(K)
    IF (Y(K) .LT. YMN) YMN = Y(K)
    IF (Y(K) .GT. YMX) YMX = Y(K)
  end do
  XMIN = XMN
  YMIN = YMN
!
! Compute cell dimensions and test for zero area.
!
  DELX = (XMX-XMN)/DBLE(NNR)
  DELY = (YMX-YMN)/DBLE(NNR)
  DX = DELX
  DY = DELY

  IF (DELX .EQ. 0.  .OR.  DELY .EQ. 0.) then
    ier = 2
    return
  end if
!
! Initialize LCELL.
!
  DO J = 1,NNR
    DO I = 1,NNR
      LCELL(I,J) = 0
    end do
  end do
!
! Loop on nodes, storing indexes in LCELL and LNEXT.
!
  DO K = NN,1,-1
    I = INT((X(K)-XMN)/DELX) + 1
    IF (I .GT. NNR) I = NNR
    J = INT((Y(K)-YMN)/DELY) + 1
    IF (J .GT. NNR) J = NNR
    L = LCELL(I,J)
    LNEXT(K) = L
    IF (L .EQ. 0) LNEXT(K) = K
    LCELL(I,J) = K
  end do

  IER = 0

  RETURN
END
SUBROUTINE TESTDT (K, N,X,Y)

!*****************************************************************************80
!
!! TESTDT returns one of five sets of test nodes.
!
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/28/97
!
!   This subroutine returns one of five sets of nodes used
! for testing scattered data fitting methods.  All five sets
! approximately cover the unit square [0,1] X [0,1]:  the
! convex hulls of sets 1 and 3 extend slightly outside the
! square but do not completely cover it, those of sets 2 and
! 5 coincide with the unit square, and the convex hull of
! set 4 is a large subset of the unit square.
!
! On input:
!
!       K = Integer in the range 1 to 5 which determines the
!           choice of data set as follows:
!
!               K = 1 - Franke's 100-node set
!               K = 2 - Franke's 33-node set
!               K = 3 - Lawson's 25-node set
!               K = 4 - Random 100-node set
!               K = 5 - Gridded 81-node set
!
!       X,Y = Arrays of length at least N(K), where
!             N(1) = 100, N(2) = 33, N(3) = 25,
!             N(4) = 100, and N(5) = 81.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       N = Number of nodes in set K, or 0 if K is outside
!           its valid range.
!
!       X,Y = Nodal coordinates of node set K.
!
! Subprograms required by TESTDT:  None
!
  implicit none

  DOUBLE PRECISION X(100), Y(100)
  INTEGER K, N
  DOUBLE PRECISION X1(100), Y1(100),  X2(33), Y2(33), &
                   X3(25), Y3(25),  X4(100), Y4(100), &
                   X5(81), Y5(81)
  INTEGER I
!
! Node set 1:  Franke's 100-node set.
!
  DATA (X1(I),Y1(I), I = 1,20)/ &
        0.0227035, -0.0310206,  0.0539888,  0.1586742, &
        0.0217008,  0.2576924,  0.0175129,  0.3414014, &
        0.0019029,  0.4943596, -0.0509685,  0.5782854, &
        0.0395408,  0.6993418, -0.0487061,  0.7470194, &
        0.0315828,  0.9107649, -0.0418785,  0.9962890, &
        0.1324189,  0.0501330,  0.1090271,  0.0918555, &
        0.1254439,  0.2592973,  0.0934540,  0.3381592, &
        0.0767578,  0.4171125,  0.1451874,  0.5615563, &
        0.0626494,  0.6552235,  0.1452734,  0.7524066, &
        0.0958668,  0.9146523,  0.0695559,  0.9632421/
  DATA (X1(I),Y1(I), I = 21,40)/ &
        0.2645602,  0.0292939,  0.2391645,  0.0602303, &
        0.2088990,  0.2668783,  0.2767329,  0.3696044, &
        0.1714726,  0.4801738,  0.2266781,  0.5940595, &
        0.1909212,  0.6878797,  0.1867647,  0.8185576, &
        0.2304634,  0.9046507,  0.2426219,  0.9805412, &
        0.3663168,  0.0396955,  0.3857662,  0.0684484, &
        0.3832392,  0.2389548,  0.3179087,  0.3124129, &
        0.3466321,  0.4902989,  0.3776591,  0.5199303, &
        0.3873159,  0.6445227,  0.3812917,  0.8203789, &
        0.3795364,  0.8938079,  0.2803515,  0.9711719/
  DATA (X1(I),Y1(I), I = 41,60)/ &
        0.4149771, -0.0284618,  0.4277679,  0.1560965, &
        0.4200010,  0.2262471,  0.4663631,  0.3175094, &
        0.4855658,  0.3891417,  0.4092026,  0.5084949, &
        0.4792578,  0.6324247,  0.4812279,  0.7511007, &
        0.3977761,  0.8489712,  0.4027321,  0.9978728, &
        0.5848691, -0.0271948,  0.5730076,  0.1272430, &
        0.6063893,  0.2709269,  0.5013894,  0.3477728, &
        0.5741311,  0.4259422,  0.6106955,  0.6084711, &
        0.5990105,  0.6733781,  0.5380621,  0.7235242, &
        0.6096967,  0.9242411,  0.5026188,  1.0308762/
  DATA (X1(I),Y1(I), I = 61,80)/ &
        0.6616928,  0.0255959,  0.6427836,  0.0707835, &
        0.6396475,  0.2008336,  0.6703963,  0.3259843, &
        0.7001181,  0.4890704,  0.6333590,  0.5096324, &
        0.6908947,  0.6697880,  0.6895638,  0.7759569, &
        0.6718889,  0.9366096,  0.6837675,  1.0064516, &
        0.7736939,  0.0285374,  0.7635332,  0.1021403, &
        0.7410424,  0.1936581,  0.8258981,  0.3235775, &
        0.7306034,  0.4714228,  0.8086609,  0.6091595, &
        0.8214531,  0.6685053,  0.7290640,  0.8022808, &
        0.8076643,  0.8476790,  0.8170951,  1.0512371/
  DATA (X1(I),Y1(I), I = 81,100)/ &
        0.8424572,  0.0380499,  0.8684053,  0.0902048, &
        0.8366923,  0.2083092,  0.9418461,  0.3318491, &
        0.8478122,  0.4335632,  0.8599583,  0.5910139, &
        0.9175700,  0.6307383,  0.8596328,  0.8144841, &
        0.9279871,  0.9042310,  0.8512805,  0.9696030, &
        1.0449820, -0.0120900,  0.9670631,  0.1334114, &
        0.9857884,  0.2695844,  0.9676313,  0.3795281, &
        1.0129299,  0.4396054,  0.9657040,  0.5044425, &
        1.0019855,  0.6941519,  1.0359297,  0.7459923, &
        1.0414677,  0.8682081,  0.9471506,  0.9801409/
!
! Node set 2:  Franke's 33-node set.
!
  DATA (X2(I),Y2(I), I = 1,33)/ &
        0.05,  0.45,  0.00,  0.50, &
        0.00,  1.00,  0.00,  0.00, &
        0.10,  0.15,  0.10,  0.75, &
        0.15,  0.30,  0.20,  0.10, &
        0.25,  0.20,  0.30,  0.35, &
        0.35,  0.85,  0.50,  0.00, &
        0.50,  1.00,  0.55,  0.95, &
        0.60,  0.25,  0.60,  0.65, &
        0.60,  0.85,  0.65,  0.70, &
        0.70,  0.20,  0.70,  0.65, &
        0.70,  0.90,  0.75,  0.10, &
        0.75,  0.35,  0.75,  0.85, &
        0.80,  0.40,  0.80,  0.65, &
        0.85,  0.25,  0.90,  0.35, &
        0.90,  0.80,  0.95,  0.90, &
        1.00,  0.00,  1.00,  0.50, &
        1.00,  1.00/
!
! Node set 3:  Lawson's 25-node set.
!
  DATA (X3(I),Y3(I), I = 1,25)/ &
        0.13750,  0.97500,   0.91250,  0.98750, &
        0.71250,  0.76250,   0.22500,  0.83750, &
       -0.05000,  0.41250,   0.47500,  0.63750, &
        0.05000, -0.05000,   0.45000,  1.03750, &
        1.08750,  0.55000,   0.53750,  0.80000, &
       -0.03750,  0.75000,   0.18750,  0.57500, &
        0.71250,  0.55000,   0.85000,  0.43750, &
        0.70000,  0.31250,   0.27500,  0.42500, &
        0.45000,  0.28750,   0.81250,  0.18750, &
        0.45000, -0.03750,   1.00000,  0.26250, &
        0.50000,  0.46250,   0.18750,  0.26250, &
        0.58750,  0.12500,   1.05000, -0.06125, &
        0.10000,  0.11250/
!
! Node set 4:  Random 100-node set.
!
  DATA (X4(I),Y4(I), I = 1,20)/ &
        0.0096326,  0.3083158,  0.0216348,  0.2450434, &
        0.0298360,  0.8613847,  0.0417447,  0.0977864, &
        0.0470462,  0.3648355,  0.0562965,  0.7156339, &
        0.0646857,  0.5311312,  0.0740377,  0.9755672, &
        0.0873907,  0.1781117,  0.0934832,  0.5452797, &
        0.1032216,  0.1603881,  0.1110176,  0.7837139, &
        0.1181193,  0.9982015,  0.1251704,  0.6910589, &
        0.1327330,  0.1049580,  0.1439536,  0.8184662, &
        0.1564861,  0.7086405,  0.1651043,  0.4456593, &
        0.1786039,  0.1178342,  0.1886405,  0.3189021/
  DATA (X4(I),Y4(I), I = 21,40)/ &
        0.2016706,  0.9668446,  0.2099886,  0.7571834, &
        0.2147003,  0.2016598,  0.2204141,  0.3232444, &
        0.2343715,  0.4368583,  0.2409660,  0.8907869, &
        0.2527740,  0.0647260,  0.2570839,  0.5692618, &
        0.2733365,  0.2947027,  0.2853833,  0.4332426, &
        0.2901755,  0.3347464,  0.2964854,  0.7436284, &
        0.3019725,  0.1066265,  0.3125695,  0.8845357, &
        0.3307163,  0.5158730,  0.3378504,  0.9425637, &
        0.3439061,  0.4799701,  0.3529922,  0.1783069, &
        0.3635507,  0.1146760,  0.3766172,  0.8225797/
  DATA (X4(I),Y4(I), I = 41,60)/ &
        0.3822429,  0.2270688,  0.3869838,  0.4073598, &
        0.3973137,  0.8875080,  0.4170708,  0.7631616, &
        0.4255588,  0.9972804,  0.4299218,  0.4959884, &
        0.4372839,  0.3410421,  0.4705033,  0.2498120, &
        0.4736655,  0.6409007,  0.4879299,  0.1058690, &
        0.4940260,  0.5411969,  0.5055324,  0.0089792, &
        0.5162593,  0.8784268,  0.5219219,  0.5515874, &
        0.5348529,  0.4038952,  0.5483213,  0.1654023, &
        0.5569571,  0.2965158,  0.5638611,  0.3660356, &
        0.5784908,  0.0366554,  0.5863950,  0.9502420/
  DATA (X4(I),Y4(I), I = 61,80)/ &
        0.5929148,  0.2638101,  0.5987839,  0.9277386, &
        0.6117561,  0.5377694,  0.6252296,  0.7374676, &
        0.6331381,  0.4674627,  0.6399048,  0.9186109, &
        0.6488972,  0.0416884,  0.6558537,  0.1291029, &
        0.6677405,  0.6763676,  0.6814074,  0.8444238, &
        0.6887812,  0.3273328,  0.6940896,  0.1893879, &
        0.7061687,  0.0645923,  0.7160957,  0.0180147, &
        0.7317445,  0.8904992,  0.7370798,  0.4160648, &
        0.7462030,  0.4688995,  0.7566957,  0.2174508, &
        0.7699998,  0.5734231,  0.7879347,  0.8853319/
  DATA (X4(I),Y4(I), I = 81,100)/ &
        0.7944014,  0.8018436,  0.8164468,  0.6388941, &
        0.8192794,  0.8931002,  0.8368405,  0.1000558, &
        0.8500993,  0.2789506,  0.8588255,  0.9082948, &
        0.8646496,  0.3259159,  0.8792329,  0.8318747, &
        0.8837536,  0.0508513,  0.8900077,  0.9708450, &
        0.8969894,  0.5120548,  0.9044917,  0.2859716, &
        0.9083947,  0.9581641,  0.9203972,  0.6183429, &
        0.9347906,  0.3779934,  0.9434519,  0.4010423, &
        0.9490328,  0.9478657,  0.9569571,  0.7425486, &
        0.9772067,  0.8883287,  0.9983493,  0.5496750/
!
! Node set 5:  9 by 9 uniform grid.
!
  DATA (X5(I),Y5(I), I = 1,20)/ &
        0.125,  0.000,  0.000,  0.125, &
        0.000,  0.250,  0.000,  0.375, &
        0.000,  0.500,  0.000,  0.625, &
        0.000,  0.750,  0.000,  0.875, &
        0.000,  1.000,  0.000,  0.000, &
        0.125,  0.125,  0.125,  0.250, &
        0.125,  0.375,  0.125,  0.500, &
        0.125,  0.625,  0.125,  0.750, &
        0.125,  0.875,  0.125,  1.000, &
        0.250,  0.000,  0.250,  0.125/
  DATA (X5(I),Y5(I), I = 21,40)/ &
        0.250,  0.250,  0.250,  0.375, &
        0.250,  0.500,  0.250,  0.625, &
        0.250,  0.750,  0.250,  0.875, &
        0.250,  1.000,  0.375,  0.000, &
        0.375,  0.125,  0.375,  0.250, &
        0.375,  0.375,  0.375,  0.500, &
        0.375,  0.625,  0.375,  0.750, &
        0.375,  0.875,  0.375,  1.000, &
        0.500,  0.000,  0.500,  0.125, &
        0.500,  0.250,  0.500,  0.375/
  DATA (X5(I),Y5(I), I = 41,60)/ &
        0.500,  0.500,  0.500,  0.625, &
        0.500,  0.750,  0.500,  0.875, &
        0.500,  1.000,  0.625,  0.000, &
        0.625,  0.125,  0.625,  0.250, &
        0.625,  0.375,  0.625,  0.500, &
        0.625,  0.625,  0.625,  0.750, &
        0.625,  0.875,  0.625,  1.000, &
        0.750,  0.000,  0.750,  0.125, &
        0.750,  0.250,  0.750,  0.375, &
        0.750,  0.500,  0.750,  0.625/
  DATA (X5(I),Y5(I), I = 61,81)/ &
        0.750,  0.750,  0.750,  0.875, &
        0.750,  1.000,  0.875,  0.000, &
        0.875,  0.125,  0.875,  0.250, &
        0.875,  0.375,  0.875,  0.500, &
        0.875,  0.625,  0.875,  0.750, &
        0.875,  0.875,  0.875,  1.000, &
        1.000,  0.000,  1.000,  0.125, &
        1.000,  0.250,  1.000,  0.375, &
        1.000,  0.500,  1.000,  0.625, &
        1.000,  0.750,  1.000,  0.875, &
        1.000,  1.000/
!
! Store node set K in (X,Y).
!
  IF (K .EQ. 1) THEN
    DO 1 I = 1,100
      X(I) = X1(I)
      Y(I) = Y1(I)
1     CONTINUE
    N = 100
  ELSEIF (K .EQ. 2) THEN
    DO 2 I = 1,33
      X(I) = X2(I)
      Y(I) = Y2(I)
2     CONTINUE
    N = 33
  ELSEIF (K .EQ. 3) THEN
    DO 3 I = 1,25
      X(I) = X3(I)
      Y(I) = Y3(I)
3     CONTINUE
    N = 25
  ELSEIF (K .EQ. 4) THEN
    DO 4 I = 1,100
      X(I) = X4(I)
      Y(I) = Y4(I)
4     CONTINUE
    N = 100
  ELSEIF (K .EQ. 5) THEN
    DO 5 I = 1,81
      X(I) = X5(I)
      Y(I) = Y5(I)
5     CONTINUE
    N = 81
  ELSE
    N = 0
  end if
  RETURN
END
SUBROUTINE TSTFN1 (K,X,Y,IFLAG, F,FX,FY)

!*****************************************************************************80
!
!! TSTFN1 computes one of ten bivariate test functions.
!
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   10/14/98
!
!   This subroutine computes the value and, optionally, the
! first partial derivatives of one of ten bivariate test
! functions.  The first six functions were chosen by Richard
! Franke to test interpolation software (See the reference
! below).  The last four functions represent more chal-
! lenging surface fitting problems.
!
! On input:
!
!       K = Integer in the range 1 to 10 which determines
!           the choice of function as follows:
!
!               K = 1 - Exponential
!               K = 2 - Cliff
!               K = 3 - Saddle
!               K = 4 - Gentle
!               K = 5 - Steep
!               K = 6 - Sphere
!               K = 7 - Trig
!               K = 8 - Synergistic Gaussian
!               K = 9 - Cloverleaf Asymmetric Peak/Valley
!               K = 10 - Cosine Peak
!
!   Note that function 6 is only defined inside a circle of
! radius 8/9 centered at (.5,.5).  Thus, if (X-.5)**2 +
! (Y-.5)**2 .GE. 64/81, the value (and partials if IFLAG=1)
! are set to 0 for this function.  Also, the first partial
! derivatives of function 10 are not defined at (.5,.5) --
! again, zeros are returned.
!
!       X,Y = Coordinates of the point at which the selected
!             function is to be evaluated.
!
!       IFLAG = Derivative option indicator:
!               IFLAG = 0 if only a function value is
!                         required.
!               IFLAG = 1 if both the function and its first
!                         partial derivatives are to be
!                         evaluated.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       F = Value of function K at (X,Y).
!
!       FX,FY = First partial derivatives of function K at
!               (X,Y) if IFLAG = 1, unaltered otherwise.
!
! Intrinsic functions called by TSTFN1:  COS, EXP, SIN,
!                                          SQRT, TANH
!
! Reference:  R. Franke, A Critical Comparison of Some
!               Methods for Interpolation of Scattered Data,
!               Naval Postgraduate School Technical Report,
!               NPS-53-79-003, 1979.
!
  implicit none

  INTEGER K, IFLAG
  DOUBLE PRECISION X, Y, F, FX, FY
  DOUBLE PRECISION T1, T2, T3, T4
  IF (K .LT. 1  .OR.  K .GT. 10) RETURN
  GO TO (1,2,3,4,5,6,7,8,9,10), K
!
! Exponential:
!
1 F = .75*EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.) + &
      .75*EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.) + &
       .5*EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.) - &
       .2*EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
  IF (IFLAG .NE. 1) RETURN
  T1 = EXP(-((9.*X-2.)**2 + (9.*Y-2.)**2)/4.)
  T2 = EXP(-((9.*X+1.)**2)/49. - (9.*Y+1.)/10.)
  T3 = EXP(-((9.*X-7.)**2 + (9.*Y-3.)**2)/4.)
  T4 = EXP(-(9.*X-4.)**2 - (9.*Y-7.)**2)
  FX = -3.375*(9.*X-2.)*T1 - (27./98.)*(9.*X+1.)*T2 &
       -2.25*(9.*X-7.)*T3 + 3.6*(9.*X-4.)*T4
  FY = -3.375*(9.*Y-2.)*T1 - .675*T2 &
       -2.25*(9.*Y-3.)*T3 + 3.6*(9.*Y-7.)*T4
  RETURN
!
! Cliff:
!
2 F = (TANH(9.0*(Y-X)) + 1.0)/9.0
  IF (IFLAG .NE. 1) RETURN
  T1 = 18.0*(Y-X)
  FX = -4.0/(EXP(T1) + 2.0 + EXP(-T1))
  FY = -FX
  RETURN
!
! Saddle:
!
3 F = (1.25 + COS(5.4*Y))/(6.0 + 6.0*(3.0*X-1.0)**2)
  IF (IFLAG .NE. 1) RETURN
  T1 = 5.4*Y
  T2 = 1.0 + (3.0*X-1.)**2
  FX = -(3.0*X-1.0)*(1.25 + COS(T1))/(T2**2)
  FY = -.9*SIN(T1)/T2
  RETURN
!
! Gentle:
!
4 F = EXP(-5.0625*((X-.5)**2 + (Y-.5)**2))/3.0
  IF (IFLAG .NE. 1) RETURN
  T1 = X - .5
  T2 = Y - .5
  T3 = -3.375*EXP(-5.0625*(T1**2 + T2**2))
  FX = T1*T3
  FY = T2*T3
  RETURN
!
! Steep:
!
5 F = EXP(-20.25*((X-.5)**2 + (Y-.5)**2))/3.0
  IF (IFLAG .NE. 1) RETURN
  T1 = X - .5
  T2 = Y - .5
  T3 = -13.5*EXP(-20.25*(T1**2 + T2**2))
  FX = T1*T3
  FY = T2*T3
  RETURN
!
! Sphere:
!
6 T4 = 64.0 - 81.0*((X-.5)**2 + (Y-.5)**2)
  F = 0.
  IF (T4 .GE. 0.) F = SQRT(T4)/9.0 - .5
  IF (IFLAG .NE. 1) RETURN
  T1 = X - .5
  T2 = Y - .5
  T3 = 0.
  IF (T4 .GT. 0.) T3 = -9.0/SQRT(T4)
  FX = T1*T3
  FY = T2*T3
  RETURN
!
! Trig:
!
7 F = 2.0*COS(10.0*X)*SIN(10.0*Y) + SIN(10.0*X*Y)
  IF (IFLAG .NE. 1) RETURN
  T1 = 10.0*X
  T2 = 10.0*Y
  T3 = 10.0*COS(10.0*X*Y)
  FX = -20.0*SIN(T1)*SIN(T2) + T3*Y
  FY = 20.0*COS(T1)*COS(T2) + T3*X
  RETURN
!
! Gaussx(1,.5,.1) + Gaussy(.75,.5,.1) + Gaussx(1,.5,.1)*
!   Gaussy(.75,.5,.1), where Gaussx(a,b,c) is the Gaussian
!   function of x with amplitude a, center (mean) b, and
!   width (standard deviation) c.
!
8 T1 = 5.0 - 10.0*X
  T2 = 5.0 - 10.0*Y
  T3 = EXP(-.5*T1*T1)
  T4 = EXP(-.5*T2*T2)
  F = T3 + .75*T4*(1.0+T3)
  IF (IFLAG .NE. 1) RETURN
  FX = T1*T3*(10.0 + 7.5*T4)
  FY = T2*T4*(7.5 + 7.5*T3)
  RETURN
!
! Cloverleaf Asymmetric Hill/Valley:
!
9 T1 = EXP((10.0 - 20.0*X)/3.0)
  T2 = EXP((10.0 - 20.0*Y)/3.0)
  T3 = 1.0/(1.0 + T1)
  T4 = 1.0/(1.0 + T2)
  F = ((20.0/3.0)**3 * T1*T2)**2 * (T3*T4)**5 * &
      (T1-2.0*T3)*(T2-2.0*T4)
  IF (IFLAG .NE. 1) RETURN
  FX = ((20.0/3.0)*T1)**2 * ((20.0/3.0)*T3)**5 * &
       (2.0*T1-3.0*T3-5.0+12.0*T3*T3)*T2*T2*T4**5 * &
       (T2-2.0*T4)
  FY = ((20.0/3.0)*T1)**2 * ((20.0/3.0)*T3)**5 * &
       (2.0*T2-3.0*T4-5.0+12.0*T4*T4)*T2*T2*T4**5 * &
       (T1-2.0*T3)
  RETURN
!
! Cosine Peak:
!
   10 T1 = SQRT( (80.0*X - 40.0)**2 + (90.0*Y - 45.0)**2 )
  T2 = EXP(-.04*T1)
  T3 = COS(.15*T1)
  F = T2*T3
  IF (IFLAG .NE. 1) RETURN
  T4 = SIN(.15*T1)
  FX = 0.
  FY = 0.
  IF (T1 .EQ. 0.) RETURN
  T4 = SIN(.15*T1)
  FX = -T2*(12.0*T4 + 3.2*T3)*(80.0*X - 40.0)/T1
  FY = -T2*(13.5*T4 + 3.6*T3)*(90.0*Y - 45.0)/T1
  RETURN
END
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
