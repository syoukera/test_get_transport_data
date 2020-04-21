*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH(I)
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  870618   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(R1MACH-S D1MACH-D),
C             MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters
C   for the local machine environment.  It is a function
C   subprogram with one (input) argument, and can be called
C   as follows, for example
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is
C   determined by the (input) value of I.  The results for
C   various values of I are discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment,
C   the desired set of DATA statements should be activated by
C   removing the C from column 1.  Also, the values of
C   D1MACH(1) - D1MACH(4) should be checked for consistency
C   with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
C     DATA SMALL(1)           / 2.23D-308              /
C     DATA LARGE(1)           / 1.79D+308              /
C     DATA RIGHT(1)           / 1.11D-16               /
C     DATA DIVER(1)           / 2.22D-16               /
C     DATA LOG10(1)           / 0.301029995663981195D0 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
      DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
      DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
      DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
      DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
      DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ('D1MACH -- I OUT OF BOUNDS', 25, 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  870618   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    6 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -125 /
C     DATA IMACH(13)/  129 /
C     DATA IMACH(14)/   53 /
C     DATA IMACH(15)/ -1021 /
C     DATA IMACH(16)/  1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     7 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -4095 /
C     DATA IMACH(13) /  4094 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -4095 /
C     DATA IMACH(16) /  4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CRAY-1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /   11 /
C     DATA IMACH( 2) /   12 /
C     DATA IMACH( 3) /    8 /
C     DATA IMACH( 4) /   10 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) /32767 /
C     DATA IMACH(10) /   16 /
C     DATA IMACH(11) /    6 /
C     DATA IMACH(12) /  -64 /
C     DATA IMACH(13) /   63 /
C     DATA IMACH(14) /   14 /
C     DATA IMACH(15) /  -64 /
C     DATA IMACH(16) /   63 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /     5/
C     DATA IMACH( 2) /     6/
C     DATA IMACH( 3) /     6/
C     DATA IMACH( 4) /     6/
C     DATA IMACH( 5) /    32/
C     DATA IMACH( 6) /     4/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    32/
C     DATA IMACH( 9) /2147483647/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -126/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    53/
C     DATA IMACH(15) / -1022/
C     DATA IMACH(16) /  1023/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /       5 /
C     DATA IMACH( 2) /       6 /
C     DATA IMACH( 3) /       0 /
C     DATA IMACH( 4) /       6 /
C     DATA IMACH( 5) /      24 /
C     DATA IMACH( 6) /       3 /
C     DATA IMACH( 7) /       2 /
C     DATA IMACH( 8) /      23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /       2 /
C     DATA IMACH(11) /      23 /
C     DATA IMACH(12) /    -127 /
C     DATA IMACH(13) /     127 /
C     DATA IMACH(14) /      38 /
C     DATA IMACH(15) /    -127 /
C     DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5/
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     39 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5 /
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     55 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH(1)  /    5 /
C     DATA IMACH(2)  /    6 /
C     DATA IMACH(3)  /    6 /
C     DATA IMACH(3)  /    7 /
C     DATA IMACH(5)  /   32 /
C     DATA IMACH(6)  /    4 /
C     DATA IMACH(7)  /    2 /
C     DATA IMACH(8)  /   32 /
C     DATA IMACH(9)  /2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -126 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   53 /
C     DATA IMACH(15) /-1015 /
C     DATA IMACH(16) / 1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     0 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -125 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1021 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
      DATA IMACH(1) /    5 /
      DATA IMACH(2) /    6 /
      DATA IMACH(3) /    6 /
      DATA IMACH(4) /    6 /
      DATA IMACH(5) /   32 /
      DATA IMACH(6) /    4 /
      DATA IMACH(7) /    2 /
      DATA IMACH(8) /   31 /
      DATA IMACH(9) /2147483647 /
      DATA IMACH(10)/    2 /
      DATA IMACH(11)/   24 /
      DATA IMACH(12)/ -125 /
      DATA IMACH(13)/  128 /
      DATA IMACH(14)/   53 /
      DATA IMACH(15)/ -1021 /
      DATA IMACH(16)/  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    1 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024 /
C     DATA IMACH(16) / 1023 /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /     1/
C     DATA IMACH( 2) /     1/
C     DATA IMACH( 3) /     0/
C     DATA IMACH( 4) /     1/
C     DATA IMACH( 5) /    16/
C     DATA IMACH( 6) /     2/
C     DATA IMACH( 7) /     2/
C     DATA IMACH( 8) /    15/
C     DATA IMACH( 9) / 32767/
C     DATA IMACH(10) /     2/
C     DATA IMACH(11) /    24/
C     DATA IMACH(12) /  -127/
C     DATA IMACH(13) /   127/
C     DATA IMACH(14) /    56/
C     DATA IMACH(15) /  -127/
C     DATA IMACH(16) /   127/
C
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
C
      STOP
      END
*DECK DAXPY
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SAXPY-S DAXPY-D CAXPY-C),
C             LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DAXPY
C
      DOUBLE PRECISION DX(1),DY(1),DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
*DECK DASUM
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C***BEGIN PROLOGUE  DASUM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A3A
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SASUM-S DASUM-D SCASUM-C),ADD,
C             LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Sum of magnitudes of d.p. vector components
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DASUM  double precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of double precision DX.
C     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DASUM
C
      DOUBLE PRECISION DX(1)
C***FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
*DECK DCOPY
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DCOPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SCOPY-S DCOPY-D CCOPY-C),COPY,
C             LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector copy y = x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DCOPY
C
      DOUBLE PRECISION DX(1),DY(1)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END
*DECK DDIFF
      DOUBLE PRECISION FUNCTION DDIFF(X,Y)
C***BEGIN PROLOGUE  DDIFF
C***REFER TO  DHFTI,DLSEI
C***ROUTINES CALLED  (NONE)
C***DESCRIPTION
C
C   ****Double Precision version of DIFF ****
C     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 June 7
C     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
C***END PROLOGUE  DDIFF
      DOUBLE PRECISION X, Y
C***FIRST EXECUTABLE STATEMENT  DDIFF
      DDIFF = X - Y
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SDOT-S DDOT-D CDOTU-C),
C             INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. inner product of d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DDOT
C
      DOUBLE PRECISION DX(1),DY(1)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
*DECK DGBCO
      SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  DGBCO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGBCO-S DGBCO-D CGBCO-C),BANDED,
C             CONDITION NUMBER,LINEAR ALGEBRA,MATRIX,
C             MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision BAND matrix by Gaussian
C            elimination and estimates the condition of the matrix.
C***DESCRIPTION
C
C     DGBCO factors a double precision band matrix by Gaussian
C     elimination and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGBFA is slightly faster.
C     To solve  A*X = B , follow DGBCO by DGBSL.
C     To compute  INVERSE(A)*C , follow DGBCO by DGBSL.
C     To compute  DETERMINANT(A) , follow DGBCO by DGBDI.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     Example:  If the original matrix is
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain
C
C            *  *  *  +  +  +  , * = not used
C            *  * 13 24 35 46  , + = used for pivoting
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and functions used:
C
C     LINPACK DGBFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,MAX0,MIN0,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGBFA,DSCAL
C***END PROLOGUE  DGBCO
      INTEGER LDA,N,ML,MU,IPVT(1)
      DOUBLE PRECISION ABD(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGBCO
      ANORM = 0.0D0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C
C     FACTOR
C
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30
            S = DABS(ABD(M,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40
            WK = WK/ABD(M,K)
            WKM = WKM/ABD(M,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML,N-K)
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = W
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150
            S = DABS(ABD(M,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0
         LM = MIN0(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
*DECK DGBFA
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C***BEGIN PROLOGUE  DGBFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGBFA-S DGBFA-D CGBFA-C),BANDED,
C             LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision BAND matrix by elimination.
C***DESCRIPTION
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C     Fortran MAX0,MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      DOUBLE PRECISION ABD(LDA,1)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGBSL
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
C***BEGIN PROLOGUE  DGBSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGBSL-S DGBSL-D CGBSL-C),BANDED,
C             LINEAR ALGEBRA,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision BAND system  A*X=B or
C            TRANS(A)*X=B using the factors computed by DGBCO or DGBFA.
C***DESCRIPTION
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C     Fortran MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      DOUBLE PRECISION ABD(LDA,1),B(1)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DGECO
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  DGECO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGECO-S DGECO-D CGECO-C),
C             CONDITION NUMBER,GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,
C             MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination
C            and estimates the condition of the matrix.
C***DESCRIPTION
C
C     DGECO factors a double precision matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGEFA is slightly faster.
C     To solve  A*X = B , follow DGECO by DGESL.
C     To compute  INVERSE(A)*C , follow DGECO by DGESL.
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
C     To compute  INVERSE(A) , follow DGECO by DGEDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an INTEGER vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
C***END PROLOGUE  DGECO
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
*DECK DGEDI
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
C***BEGIN PROLOGUE  DGEDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D3A1,D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGEDI-S DGEDI-D CGEDI-C),
C             DETERMINANT,INVERSE,LINEAR ALGEBRA,MATRIX,
C             MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant and inverse of a matrix using
C            factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGEDI computes the determinant and inverse of a matrix
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
C        INFO .EQ. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,DSWAP
C     Fortran DABS,MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,DSWAP
C***END PROLOGUE  DGEDI
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),DET(2),WORK(1)
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  DGEDI
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
*DECK DGEFA
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
C***BEGIN PROLOGUE  DGEFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGEFA-S DGEFA-D CGEFA-C),
C             GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination.
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
C***BEGIN PROLOGUE  DGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),
C             TYPE=DOUBLE PRECISION(SGESL-S DGESL-D CGESL-C),
C             LINEAR ALGEBRA,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
C            using the factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DH12
      SUBROUTINE DH12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C***BEGIN PROLOGUE  DH12
C***REFER TO  DHFTI,DLSEI,DWNNLS
C***ROUTINES CALLED  DAXPY,DDOT,DSWAP
C***DESCRIPTION
C
C      *** DOUBLE PRECISION VERSION OF H12 ******
C
C     SUBROUTINE DH12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C
C     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
C     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
C
C     Modified at SANDIA LABS, May 1977, to --
C
C     1)  Remove double precision accumulation, and
C     2)  Include usage of the Basic Linear Algebra Package for
C         vectors longer than a particular threshold.
C
C     Construction and/or application of a single
C     Householder transformation..     Q = I + U*(U**T)/B
C
C     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
C     LPIVOT is the index of the pivot element.
C     L1,M   If L1 .LE. M   the transformation will be constructed to
C            zero elements indexed from L1 through M.   If L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
C                   IUE is the storage increment between elements.
C                                       On exit from H1 U() and UP
C                   contain quantities defining the vector U of the
C                   Householder transformation.   On entry to H2 U()
C                   and UP should contain quantities previously computed
C                   by H1.  These will not be modified by H2.
C     C()    On entry to H1 or H2 C() contains a matrix which will be
C            regarded as a set of vectors to which the Householder
C            transformation is to be applied.  On exit C() contains the
C            set of transformed vectors.
C     ICE    Storage increment between elements of vectors in C().
C     ICV    Storage increment between vectors in C().
C     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
C            no operations will be done on C().
C***END PROLOGUE  DH12
      INTEGER I, I2, I3, I4, ICE, ICV, INCR, IUE, J, KL1, KL2, KLP,
     *     L1, L1M1, LPIVOT, M, MML1P2, MODE, NCV
      DOUBLE PRECISION B, C, CL, CLINV, ONE, UL1M1, SM, U, UP, DDOT
      DIMENSION U(IUE,M), C(*)
C     BEGIN BLOCK PERMITTING ...EXITS TO 140
C***FIRST EXECUTABLE STATEMENT  DH12
         ONE = 1.0D0
C
C     ...EXIT
         IF (0 .GE. LPIVOT .OR. LPIVOT .GE. L1 .OR. L1 .GT. M) GO TO 140
         CL = DABS(U(1,LPIVOT))
         IF (MODE .EQ. 2) GO TO 40
C           ****** CONSTRUCT THE TRANSFORMATION. ******
            DO 10 J = L1, M
               CL = DMAX1(DABS(U(1,J)),CL)
   10       CONTINUE
            IF (CL .GT. 0.0D0) GO TO 20
C     .........EXIT
               GO TO 140
   20       CONTINUE
            CLINV = ONE/CL
            SM = (U(1,LPIVOT)*CLINV)**2
            DO 30 J = L1, M
               SM = SM + (U(1,J)*CLINV)**2
   30       CONTINUE
            CL = CL*DSQRT(SM)
            IF (U(1,LPIVOT) .GT. 0.0D0) CL = -CL
            UP = U(1,LPIVOT) - CL
            U(1,LPIVOT) = CL
         GO TO 50
   40    CONTINUE
C        ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
         IF (CL .GT. 0.0D0) GO TO 50
C     ......EXIT
            GO TO 140
   50    CONTINUE
C     ...EXIT
         IF (NCV .LE. 0) GO TO 140
         B = UP*U(1,LPIVOT)
C        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
         IF (B .LT. 0.0D0) GO TO 60
C     ......EXIT
            GO TO 140
   60    CONTINUE
         B = ONE/B
         MML1P2 = M - L1 + 2
         IF (MML1P2 .LE. 20) GO TO 80
            L1M1 = L1 - 1
            KL1 = 1 + (L1M1 - 1)*ICE
            KL2 = KL1
            KLP = 1 + (LPIVOT - 1)*ICE
            UL1M1 = U(1,L1M1)
            U(1,L1M1) = UP
            IF (LPIVOT .NE. L1M1) CALL DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
            DO 70 J = 1, NCV
               SM = DDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
               SM = SM*B
               CALL DAXPY(MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
               KL1 = KL1 + ICV
   70       CONTINUE
            U(1,L1M1) = UL1M1
C     ......EXIT
            IF (LPIVOT .EQ. L1M1) GO TO 140
            KL1 = KL2
            CALL DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
         GO TO 130
   80    CONTINUE
            I2 = 1 - ICV + ICE*(LPIVOT - 1)
            INCR = ICE*(L1 - LPIVOT)
            DO 120 J = 1, NCV
               I2 = I2 + ICV
               I3 = I2 + INCR
               I4 = I3
               SM = C(I2)*UP
               DO 90 I = L1, M
                  SM = SM + C(I3)*U(1,I)
                  I3 = I3 + ICE
   90          CONTINUE
               IF (SM .EQ. 0.0D0) GO TO 110
                  SM = SM*B
                  C(I2) = C(I2) + SM*UP
                  DO 100 I = L1, M
                     C(I4) = C(I4) + SM*U(1,I)
                     I4 = I4 + ICE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
*DECK DHFTI
      SUBROUTINE DHFTI(A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)
C***BEGIN PROLOGUE  DHFTI
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D9
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(HFTI-S DHFTI-D),
C             CURVE FITTING,LEAST SQUARES
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C***PURPOSE  Solve the least squares problem Ax = b for
C            banded matrices A using sequential accumulation of rows of
C            the data matrix. Exactly one right-handed side vector is
C            permitted.
C***DESCRIPTION
C
C       ***** Double Precision version of HFTI *****
C
C     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
C
C     This subroutine solves a linear least squares problem or a set of
C     linear least squares problems having the same matrix but different
C     right-side vectors.  The problem data consists of an M by N matrix
C     A, an M by NB matrix B, and an absolute tolerance parameter TAU
C     whose usage is described below.  The NB column vectors of B
C     represent right-side vectors for NB distinct linear least squares
C     problems.
C
C     This set of problems can also be written as the matrix least
C     squares problem
C
C                       AX = B,
C
C     where X is the N by NB solution matrix.
C
C     Note that if B is the M by M identity matrix, then X will be the
C     pseudo-inverse of A.
C
C     This subroutine first transforms the augmented matrix (A B) to a
C     matrix (R C) using premultiplying Householder transformations with
C     column interchanges.  All subdiagonal elements in the matrix R are
C     zero and its diagonal elements satisfy
C
C                       DABS(R(I,I)).GE.DABS(R(I+1,I+1)),
C
C                       I = 1,...,L-1, where
C
C                       L = MIN(M,N).
C
C     The subroutine will compute an integer, KRANK, equal to the number
C     of diagonal terms of R that exceed TAU in magnitude. Then a
C     solution of minimum Euclidean length is computed using the first
C     KRANK rows of (R C).
C
C     To be specific we suggest that the user consider an easily
C     computable matrix norm, such as, the maximum of all column sums of
C     magnitudes.
C
C     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
C     norm of B), it is suggested that TAU be set approximately equal to
C     EPS*(norm of A).
C
C     The user must dimension all arrays appearing in the call list..
C     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
C     permits the solution of a range of problems in the same array
C     space.
C
C     The entire set of parameters for DHFTI are
C
C     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
C                       matrix A of the least squares problem AX = B.
C                       The first dimensioning parameter of the array
C                       A(*,*) is MDA, which must satisfy MDA.GE.M
C                       Either M.GE.N or M.LT.N is permitted.  There
C                       is no restriction on the rank of A.  The
C                       condition MDA.LT.M is considered an error.
C
C     B(*),MDB,NB       If NB = 0 the subroutine will perform the
C                       orthogonal decomposition but will make no
C                       references to the array B(*).  If NB.GT.0
C                       the array B(*) must initially contain the M by
C                       NB matrix B of the least squares problem AX =
C                       B.  If NB.GE.2 the array B(*) must be doubly
C                       subscripted with first dimensioning parameter
C                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
C                       be either doubly or singly subscripted.  In
C                       the latter case the value of MDB is arbitrary
C                       but it should be set to some valid integer
C                       value such as MDB = M.
C
C                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
C                       is considered an error.
C
C     TAU               Absolute tolerance parameter provided by user
C                       for pseudorank determination.
C
C     H(*),G(*),IP(*)   Arrays of working space used by DHFTI.
C
C     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     A(*,*)            The contents of the array A(*,*) will be
C                       modified by the subroutine. These contents
C                       are not generally required by the user.
C
C     B(*)              On return the array B(*) will contain the N by
C                       NB solution matrix X.
C
C     KRANK             Set by the subroutine to indicate the
C                       pseudorank of A.
C
C     RNORM(*)          On return, RNORM(J) will contain the Euclidean
C                       norm of the residual vector for the problem
C                       defined by the J-th column vector of the array
C                       B(*,*) for J = 1,...,NB.
C
C     H(*),G(*)         On return these arrays respectively contain
C                       elements of the pre- and post-multiplying
C                       Householder transformations used to compute
C                       the minimum Euclidean length solution.
C
C     IP(*)             Array in which the subroutine records indices
C                       describing the permutation of column vectors.
C                       The contents of arrays H(*),G(*) and IP(*)
C                       are not generally required by the user.
C***REFERENCES  C. L. LAWSON AND R. J. HANSON,
C                 SOLVING LEAST SQUARES PROBLEMS, PRENTICE-HALL,INC
C                 (1974), CHAPTER 14.
C***ROUTINES CALLED  DDIFF,DH12,XERROR
C***END PROLOGUE  DHFTI
      INTEGER I, II, IOPT, IP(N), IP1, J, JB, JJ, K, KP1, KRANK, L,
     *     LDIAG, LMAX, M, MAX0, MDA, MDB, MIN0, N, NB, NERR
      DOUBLE PRECISION A, B, DABS, DBLE, DDIFF, DSQRT, DZERO, FACTOR,
     *     G, H, HMAX, RNORM, SM, SM1, SZERO, TAU, TMP
      DIMENSION A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB)
C     BEGIN BLOCK PERMITTING ...EXITS TO 360
C***FIRST EXECUTABLE STATEMENT  DHFTI
         SZERO = 0.0D0
         DZERO = 0.0D0
         FACTOR = 0.001D0
C
         K = 0
         LDIAG = MIN0(M,N)
         IF (LDIAG .LE. 0) GO TO 350
C           BEGIN BLOCK PERMITTING ...EXITS TO 130
C              BEGIN BLOCK PERMITTING ...EXITS TO 120
                  IF (MDA .GE. M) GO TO 10
                     NERR = 1
                     IOPT = 2
                     CALL XERROR('DHFTI MDA.LT.M.. PROBABLE ERROR.',33,
     *                           NERR,IOPT)
C     ...............EXIT
                     GO TO 360
   10             CONTINUE
C
                  IF (NB .LE. 1 .OR. MAX0(M,N) .LE. MDB) GO TO 20
                     NERR = 2
                     IOPT = 2
                     CALL XERROR(
     *               ' DHFTI MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERRO
     *R.',51,NERR,IOPT)
C     ...............EXIT
                     GO TO 360
   20             CONTINUE
C
                  DO 100 J = 1, LDIAG
C                    BEGIN BLOCK PERMITTING ...EXITS TO 70
                        IF (J .EQ. 1) GO TO 40
C
C                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C                          ..
                           LMAX = J
                           DO 30 L = J, N
                              H(L) = H(L) - A(J-1,L)**2
                              IF (H(L) .GT. H(LMAX)) LMAX = L
   30                      CONTINUE
C                    ......EXIT
                           IF (DDIFF(HMAX+FACTOR*H(LMAX),HMAX)
     *                         .GT. 0.0D0) GO TO 70
   40                   CONTINUE
C
C                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C                       ..
                        LMAX = J
                        DO 60 L = J, N
                           H(L) = 0.0D0
                           DO 50 I = J, M
                              H(L) = H(L) + A(I,L)**2
   50                      CONTINUE
                           IF (H(L) .GT. H(LMAX)) LMAX = L
   60                   CONTINUE
                        HMAX = H(LMAX)
   70                CONTINUE
C                    ..
C                     LMAX HAS BEEN DETERMINED
C
C                     DO COLUMN INTERCHANGES IF NEEDED.
C                    ..
                     IP(J) = LMAX
                     IF (IP(J) .EQ. J) GO TO 90
                        DO 80 I = 1, M
                           TMP = A(I,J)
                           A(I,J) = A(I,LMAX)
                           A(I,LMAX) = TMP
   80                   CONTINUE
                        H(LMAX) = H(J)
   90                CONTINUE
C
C                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
C                     AND B.
C                    ..
                     CALL DH12(1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,
     *                         N-J)
                     CALL DH12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
  100             CONTINUE
C
C                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
C                  TAU.
C                 ..
                  DO 110 J = 1, LDIAG
C              ......EXIT
                     IF (DABS(A(J,J)) .LE. TAU) GO TO 120
  110             CONTINUE
                  K = LDIAG
C           ......EXIT
                  GO TO 130
  120          CONTINUE
               K = J - 1
  130       CONTINUE
            KP1 = K + 1
C
C           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
            IF (NB .LT. 1) GO TO 170
            DO 160 JB = 1, NB
               TMP = SZERO
               IF (M .LT. KP1) GO TO 150
               DO 140 I = KP1, M
                  TMP = TMP + B(I,JB)**2
  140          CONTINUE
  150          CONTINUE
               RNORM(JB) = DSQRT(TMP)
  160       CONTINUE
  170       CONTINUE
C           SPECIAL FOR PSEUDORANK = 0
            IF (K .GT. 0) GO TO 210
               IF (NB .LT. 1) GO TO 200
               DO 190 JB = 1, NB
                  DO 180 I = 1, N
                     B(I,JB) = SZERO
  180             CONTINUE
  190          CONTINUE
  200          CONTINUE
            GO TO 340
  210       CONTINUE
C
C               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C               DECOMPOSITION OF FIRST K ROWS.
C              ..
               IF (K .EQ. N) GO TO 230
                  DO 220 II = 1, K
                     I = KP1 - II
                     CALL DH12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  220             CONTINUE
  230          CONTINUE
C
C
               IF (NB .LT. 1) GO TO 330
               DO 320 JB = 1, NB
C
C                  SOLVE THE K BY K TRIANGULAR SYSTEM.
C                 ..
                  DO 260 L = 1, K
                     SM = DZERO
                     I = KP1 - L
                     IP1 = I + 1
                     IF (K .LT. IP1) GO TO 250
                     DO 240 J = IP1, K
                        SM = SM + A(I,J)*B(J,JB)
  240                CONTINUE
  250                CONTINUE
                     SM1 = SM
                     B(I,JB) = (B(I,JB) - SM1)/A(I,I)
  260             CONTINUE
C
C                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
C                 ..
                  IF (K .EQ. N) GO TO 290
                     DO 270 J = KP1, N
                        B(J,JB) = SZERO
  270                CONTINUE
                     DO 280 I = 1, K
                        CALL DH12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,
     *                            MDB,1)
  280                CONTINUE
  290             CONTINUE
C
C                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C                   COLUMN INTERCHANGES.
C                 ..
                  DO 310 JJ = 1, LDIAG
                     J = LDIAG + 1 - JJ
                     IF (IP(J) .EQ. J) GO TO 300
                        L = IP(J)
                        TMP = B(L,JB)
                        B(L,JB) = B(J,JB)
                        B(J,JB) = TMP
  300                CONTINUE
  310             CONTINUE
  320          CONTINUE
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
C        ..
C         THE SOLUTION VECTORS, X, ARE NOW
C         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
C
         KRANK = K
  360 CONTINUE
      RETURN
      END
*DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2(N,DX,INCX)
C***BEGIN PROLOGUE  DNRM2
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A3B
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SNRM2-S DNRM2-D SCNRM2-C),
C             EUCLIDEAN LENGTH,EUCLIDEAN NORM,L2,LINEAR ALGEBRA,UNITARY,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Euclidean length (L2 norm) of d.p. vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX() with storage
C     increment INCX .
C     If    N .LE. 0 return with result = 0.
C     If N .GE. 1 then INCX must be .GE. 1
C
C           C.L. Lawson, 1978 Jan 08
C
C     Four phase method     using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  DSQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  DSQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm..
C
C     Phase 1    scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI..
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as followS..
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
      SAVE CUTLO, CUTHI, ZERO, ONE
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DNRM2
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C***FIRST EXECUTABLE STATEMENT  DNRM2
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
*DECK DP1VLU
      SUBROUTINE DP1VLU(L,NDER,X,YFIT,YP,A)
C***BEGIN PROLOGUE  DP1VLU
C***DATE WRITTEN   740601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  K6
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(PVALUE-S DP1VLU-D),
C             CURVE FITTING,LEAST SQUARES,POLYNOMIAL APPROXIMATION
C***AUTHOR  SHAMPINE, L. F., (SNLA)
C           DAVENPORT, S. M., (SNLA)
C***PURPOSE  Use the coefficients generated by DPOLFT to evaluate the
C            polynomial fit of degree L, along with the first NDER of
C            its derivatives, at a specified point.
C***DESCRIPTION
C
C       *** Double Precision Version of PVALUE ***
C
C     Written by L. F. Shampine and S. M. Davenport.
C
C     Abstract
C
C     The subroutine  DP1VLU  uses the coefficients generated by  DPOLFT
C     to evaluate the polynomial fit of degree  L , along with the first
C     NDER  of its derivatives, at a specified point.  Computationally
C     stable recurrence relations are used to perform this task.
C
C     The parameters for  DP1VLU  are
C
C     Input -- ALL TYPE REAL variables are DOUBLE PRECISION
C         L -      the degree of polynomial to be evaluated.  L  may be
C                  any non-negative integer which is less than or equal
C                  to  NDEG , the highest degree polynomial provided
C                  by  DPOLFT .
C         NDER -   the number of derivatives to be evaluated.  NDER
C                  may be 0 or any positive value.  If NDER is less
C                  than 0, it will be treated as 0.
C         X -      the argument at which the polynomial and its
C                  derivatives are to be evaluated.
C         A -      work and output array containing values from last
C                  call to  DPOLFT .
C
C     Output -- ALL TYPE REAL variables are DOUBLE PRECISION
C         YFIT -   value of the fitting polynomial of degree  L  at  X
C         YP -     array containing the first through  NDER  derivatives
C                  of the polynomial of degree  L .  YP  must be
C                  dimensioned at least  NDER  in the calling program.
C***REFERENCES  SHAMPINE L.F., DAVENPORT S.M., HUDDLESTON R.E., *CURVE
C                 FITTING BY POLYNOMIALS IN ONE VARIABLE*, SLA-74-0270,
C                 SANDIA LABORATORIES, JUNE 1974.
C***ROUTINES CALLED  XERROR,XERRWV
C***END PROLOGUE  DP1VLU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I,IC,ILO,IN,INP1,IUP,K1,K1I,K2,K3,K3P1,K3PN,K4,K4P1,K4PN,
     * KC,L,LM1,LP1,MAXORD,N,NDER,NDO,NDP1,NORD
      DOUBLE PRECISION A(*),CC,DIF,VAL,X,YFIT,YP(*)
C***FIRST EXECUTABLE STATEMENT  DP1VLU
      IF (L .LT. 0) GO TO 12
      NDO = MAX0(NDER,0)
      NDO = MIN0(NDO,L)
      MAXORD = A(1) + 0.5D0
      K1 = MAXORD + 1
      K2 = K1 + MAXORD
      K3 = K2 + MAXORD + 2
      NORD = A(K3) + 0.5D0
      IF (L .GT. NORD) GO TO 11
      K4 = K3 + L + 1
      IF (NDER .LT. 1) GO TO 2
      DO 1 I = 1,NDER
 1      YP(I) = 0.0D0
 2    IF (L .GE. 2) GO TO 4
      IF (L .EQ. 1) GO TO 3
C
C L IS 0
C
      VAL = A(K2+1)
      GO TO 10
C
C L IS 1
C
 3    CC = A(K2+2)
      VAL = A(K2+1) + (X-A(2))*CC
      IF (NDER .GE. 1) YP(1) = CC
      GO TO 10
C
C L IS GREATER THAN 1
C
 4    NDP1 = NDO + 1
      K3P1 = K3 + 1
      K4P1 = K4 + 1
      LP1 = L + 1
      LM1 = L - 1
      ILO = K3 + 3
      IUP = K4 + NDP1
      DO 5 I = ILO,IUP
 5      A(I) = 0.0D0
      DIF = X - A(LP1)
      KC = K2 + LP1
      A(K4P1) = A(KC)
      A(K3P1) = A(KC-1) + DIF*A(K4P1)
      A(K3+2) = A(K4P1)
C
C EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
C
      DO 9 I = 1,LM1
        IN = L - I
        INP1 = IN + 1
        K1I = K1 + INP1
        IC = K2 + IN
        DIF = X - A(INP1)
        VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
        IF (NDO .LE. 0) GO TO 8
        DO 6 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
 6        YP(N) = DIF*A(K3PN) + DBLE(N)*A(K3PN-1) - A(K1I)*A(K4PN)
C
C SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
C
        DO 7 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
          A(K4PN) = A(K3PN)
 7        A(K3PN) = YP(N)
 8      A(K4P1) = A(K3P1)
 9      A(K3P1) = VAL
C
C NORMAL RETURN OR ABORT DUE TO ERROR
C
 10   YFIT = VAL
      RETURN
   11 CALL XERRWV ( 'DP1VLU-THE ORDER OF POLYNOMIAL EVALUATION, L(I1), R
     1EQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD(I2), COMPUTED BY POLF
     2IT -- EXECUTION TERMINATED.',144,8,2,2,L,NORD,0,0.0,0.0)
      RETURN
   12 CALL XERROR ( 'DP1VLU-INVALID INPUT PARAMETER.  ORDER OF POLYNOMIA
     1L EVALUATION REQUESTED IS NEGATIVE -- EXECUTION TERMINATED.',110,
     2 2,2)
      RETURN
      END
*DECK DPCOEF
      SUBROUTINE DPCOEF(L,C,TC,A)
C***BEGIN PROLOGUE  DPCOEF
C***DATE WRITTEN   740601   (YYMMDD)
C***REVISION DATE  870701   (YYMMDD)
C***CATEGORY NO.  K1A1A2
C***KEYWORDS  CURVE FITTING,DATA FITTING,LEAST SQUARES,POLYNOMIAL FIT,
C             LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(PCOEF-S DPCOEF-D)
C***AUTHOR  SHAMPINE, L. F., (SNLA)
C           DAVENPORT, S. M., (SNLA)
C***PURPOSE  Convert the DPOLFT coefficients to Taylor series form.
C***DESCRIPTION
C
C  **** Double Precision version of PCOEF ****
C     Written BY L. F. Shampine and S. M. Davenport.
C
C     Abstract
C
C     DPOLFT  computes the least squares polynomial fit of degree  L  as
C     a sum of orthogonal polynomials.  DPCOEF  changes this fit to its
C     Taylor expansion about any point  C , i.e. writes the polynomial
C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
C     in powers of X, but a suitable non-zero  C  often leads to
C     polynomials which are better scaled and more accurately evaluated.
C
C     The parameters for  DPCOEF  are
C
C     INPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         L -      Indicates the degree of polynomial to be changed to
C                  its Taylor expansion.  To obtain the Taylor
C                  coefficients in reverse order, input  L  as the
C                  negative of the degree desired.  The absolute value
C                  of L  must be less than or equal to NDEG, the highest
C                  degree polynomial fitted by  DPOLFT .
C         C -      The point about which the Taylor expansion is to be
C                  made.
C         A -      Work and output array containing values from last
C                  call to  DPOLFT .
C
C     OUTPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         TC -     Vector containing the first LL+1 Taylor coefficients
C                  where LL=IABS(L).  If  L.GT.0 , the coefficients are
C                  in the usual Taylor series order, i.e.
C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
C                  If L .LT. 0, the coefficients are in reverse order,
C                  i.e.
C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
C
C***REFERENCES  SHAMPINE L.F., DAVENPORT S.M., HUDDLESTON R.E., *CURVE
C                 FITTING BY POLYNOMIALS IN ONE VARIABLE*, SLA-74-
C                 0270, SANDIA LABORATORIES, JUNE 1974.
C***ROUTINES CALLED  DP1VLU
C***END PROLOGUE  DPCOEF
C
      INTEGER I,L,LL,LLP1,LLP2,NEW,NR
      DOUBLE PRECISION A(*),C,FAC,SAVE,TC(*)
C***FIRST EXECUTABLE STATEMENT  DPCOEF
      LL = IABS(L)
      LLP1 = LL + 1
      CALL DP1VLU (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = 1.0D0
      DO 1 I = 3,LLP1
        FAC = FAC*FLOAT(I-1)
 1      TC(I) = TC(I)/FAC
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO 3 I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
 3      TC(NEW) = SAVE
 4    RETURN
      END
*DECK DPOLFT
      SUBROUTINE DPOLFT(N,X,Y,W,MAXDEG,NDEG,EPS,R,IERR,A)
C***BEGIN PROLOGUE  DPOLFT
C***DATE WRITTEN   740601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  K1A1A2
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(POLFIT-S DPOLFT-D),
C             CURVE FITTING,DATA FITTING,LEAST SQUARES,POLYNOMIAL FIT
C***AUTHOR  SHAMPINE, L. F., (SNLA)
C           DAVENPORT, S. M., (SNLA)
C           HUDDLESTON, R. E., (SNLL)
C***PURPOSE  Fit discrete data in a least squares sense by polynomials
C            in one variable.
C***DESCRIPTION
C
C       *** Double Precision version of POLFIT ***
C
C     Written by L. F. Shampine and S. M.  Davenport.  The statistical
C     options provided were written by R. E. Huddleston.
C
C     Abstract
C
C     Given a collection of points X(I) and a set of values Y(I) which
C     correspond to some function or measurement at each of the X(I),
C     subroutine  DPOLFT  computes the weighted least-squares polynomial
C     fits of all degrees up to some degree either specified by the user
C     or determined by the routine.  The fits thus obtained are in
C     orthogonal polynomial form.  Subroutine  DP1VLU  may then be
C     called to evaluate the fitted polynomials and any of their
C     derivatives at any point.  The subroutine  DPCOEF  may be used to
C     express the polynomial fits as powers of (X-C) for any specified
C     point C.
C
C     The parameters for  DPOLFT  are
C
C     Input -- All TYPE REAL variables are DOUBLE PRECISION
C         N -      the number of data points.  The arrays X, Y, W, R
C                  must be dimensioned at least  N  (N .GE. 1).
C         X -      array of values of the independent variable.  These
C                  values may appear in any order and need not all be
C                  distinct.
C         Y -      array of corresponding function values.
C         W -      array of positive values to be used as weights.  If
C                  W(1) is negative,  DPOLFT  will set all the weights
C                  to 1.0, which means absolute error will be minimized.
C                  to minimize relative error, the user should set
C                  weights to:  W(I) = 1.0/Y(I)**2, I = 1,...,N .
C         MAXDEG - maximum degree to be allowed for polynomial fit.
C                  MAXDEG  may be any non-negative integer less than  N.
C                  Note -- MAXDEG  cannot be equal to  N-1  when a
C                  statistical test is to be used for degree selection,
C                  i.e., when input value of  EPS  is negative.
C         EPS -    specifies the criterion to be used in determining
C                  the degree of fit to be computed.
C                  (1)  If  EPS  is input negative,  DPOLFT  chooses the
C                       degree based on a statistical F test of
C                       significance.  One of three possible
C                       significance levels will be used:  .01, .05 or
C                       .10.  If  EPS=-1.0 , the routine will
C                       automatically select one of these levels based
C                       on the number of data points and the maximum
C                       degree to be considered.  If  EPS  is input as
C                       -.01, -.05, or -.10, a significance level of
C                       .01, .05, or .10, respectively, will be used.
C                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
C                       polynomials of degrees 0 through  MAXDEG .
C                  (3)  If  EPS  is input positive,  EPS  is the RMS
C                       error tolerance which must be satisfied by the
C                       fitted polynomial.  DPOLFT  will increase the
C                       degree of fit until this criterion is met or
C                       until the maximum degree is reached.
C
C     Output -- All TYPE REAL variables are DOUBLE PRECISION
C         NDEG -   degree of the highest degree fit computed.
C         EPS -    RMS error of the polynomial of degree  NDEG .
C         R -      vector containing values of the fit of degree  NDEG
C                  at each of the  X(I) .  Except when the statistical
C                  test is used, these values are more accurate than
C                  results from subroutine  DP1VLU  normally are.
C         IERR -   error flag with the following possible values.
C             1 -- indicates normal execution, i.e., either
C                  (1)  the input value of  EPS  was negative, and the
C                       computed polynomial fit of degree  NDEG
C                       satisfies the specified F test, or
C                  (2)  the input value of  EPS  was 0., and the fits of
C                       all degrees up to  MAXDEG  are complete, or
C                  (3)  the input value of  EPS  was positive, and the
C                       polynomial of degree  NDEG  satisfies the RMS
C                       error requirement.
C             2 -- invalid input parameter.  At least one of the input
C                  parameters has an illegal value and must be corrected
C                  before  DPOLFT  can proceed.  Valid input results
C                  when the following restrictions are observed
C                       N .GE. 1
C                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
C                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
C                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
C             3 -- cannot satisfy the RMS error requirement with a
C                  polynomial of degree no greater than  MAXDEG .  Best
C                  fit found is of degree  MAXDEG .
C             4 -- cannot satisfy the test for significance using
C                  current value of  MAXDEG .  Statistically, the
C                  best fit found is of order  NORD .  (In this case,
C                  NDEG will have one of the values:  MAXDEG-2,
C                  MAXDEG-1, or MAXDEG).  Using a higher value of
C                  MAXDEG  may result in passing the test.
C         A -      work and output array having at least 3N+3MAXDEG+3
C                  locations
C
C     Note - DPOLFT  calculates all fits of degrees up to and including
C            NDEG .  Any or all of these fits can be evaluated or
C            expressed as powers of (X-C) using  DP1VLU  and  DPCOEF
C            after just one call to  DPOLFT .
C***REFERENCES  SHAMPINE L.F., DAVENPORT S.M., HUDDLESTON R.E., *CURVE
C                 FITTING BY POLYNOMIALS IN ONE VARIABLE*, SLA-74-
C                 0270, SANDIA LABORATORIES, JUNE 1974.
C***ROUTINES CALLED  DP1VLU,XERROR
C***END PROLOGUE  DPOLFT
      INTEGER I,IDEGF,IERR,J,JP1,JPAS,K,K1,K1PJ,K2,K2PJ,K3,K3PI,K4,
     * K4PI,K5,K5PI,KSIG,M,MAXDEG,MOP1,NDEG,NDER,NFAIL
      DOUBLE PRECISION TEMD1,TEMD2
      DOUBLE PRECISION A(*),DEGF,DEN,EPS,ETST,F,FCRIT,R(*),SIGJ,
     * SIGJM1,SIGPAS,TEMP,X(*),XM,Y(*),W(*),W1,W11
      DOUBLE PRECISION CO(4,3)
      DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2),
     1      CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3),
     2  CO(4,3)/-13.086850,-2.4648165,-3.3846535,-1.2973162,
     3          -3.3381146,-1.7812271,-3.2578406,-1.6589279,
     4          -1.6282703,-1.3152745,-3.2640179,-1.9829776/
C***FIRST EXECUTABLE STATEMENT  DPOLFT
      M = IABS(N)
      IF (M .EQ. 0) GO TO 30
      IF (MAXDEG .LT. 0) GO TO 30
      A(1) = MAXDEG
      MOP1 = MAXDEG + 1
      IF (M .LT. MOP1) GO TO 30
      IF (EPS .LT. 0.0D0 .AND.  M .EQ. MOP1) GO TO 30
      XM = M
      ETST = EPS*EPS*XM
      IF (W(1) .LT. 0.0D0) GO TO 2
      DO 1 I = 1,M
        IF (W(I) .LE. 0.0D0) GO TO 30
 1      CONTINUE
      GO TO 4
 2    DO 3 I = 1,M
 3      W(I) = 1.0D0
 4    IF (EPS .GE. 0.0D0) GO TO 8
C
C DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
C CHOOSING DEGREE OF POLYNOMIAL FIT
C
      IF (EPS .GT. (-.55D0)) GO TO 5
      IDEGF = M - MAXDEG - 1
      KSIG = 1
      IF (IDEGF .LT. 10) KSIG = 2
      IF (IDEGF .LT. 5) KSIG = 3
      GO TO 8
 5    KSIG = 1
      IF (EPS .LT. (-.03D0)) KSIG = 2
      IF (EPS .LT. (-.07D0)) KSIG = 3
C
C INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
C
 8    K1 = MAXDEG + 1
      K2 = K1 + MAXDEG
      K3 = K2 + MAXDEG + 2
      K4 = K3 + M
      K5 = K4 + M
      DO 9 I = 2,K4
 9      A(I) = 0.0D0
      W11 = 0.0D0
      IF (N .LT. 0) GO TO 11
C
C UNCONSTRAINED CASE
C
      DO 10 I = 1,M
        K4PI = K4 + I
        A(K4PI) = 1.0D0
 10     W11 = W11 + W(I)
      GO TO 13
C
C CONSTRAINED CASE
C
 11   DO 12 I = 1,M
        K4PI = K4 + I
 12     W11 = W11 + W(I)*A(K4PI)**2
C
C COMPUTE FIT OF DEGREE ZERO
C
 13   TEMD1 = 0.0D0
      DO 14 I = 1,M
        K4PI = K4 + I
        TEMD1 = TEMD1 + W(I)*Y(I)*A(K4PI)
 14     CONTINUE
      TEMD1 = TEMD1/W11
      A(K2+1) = TEMD1
      SIGJ = 0.0D0
      DO 15 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
 15     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
      J = 0
C
C SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
C
      IF (EPS) 24,26,27
C
C INCREMENT DEGREE
C
 16   J = J + 1
      JP1 = J + 1
      K1PJ = K1 + J
      K2PJ = K2 + J
      SIGJM1 = SIGJ
C
C COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
C
      IF (J .GT. 1) A(K1PJ) = W11/W1
C
C COMPUTE NEW A COEFFICIENT
C
      TEMD1 = 0.0D0
      DO 18 I = 1,M
        K4PI = K4 + I
        TEMD2 = A(K4PI)
        TEMD1 = TEMD1 + X(I)*W(I)*TEMD2*TEMD2
 18     CONTINUE
      A(JP1) = TEMD1/W11
C
C EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
C
      W1 = W11
      W11 = 0.0D0
      DO 19 I = 1,M
        K3PI = K3 + I
        K4PI = K4 + I
        TEMP = A(K3PI)
        A(K3PI) = A(K4PI)
        A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
 19     W11 = W11 + W(I)*A(K4PI)**2
C
C GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
C PRECISION
C
      TEMD1 = 0.0D0
      DO 20 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = W(I)*((Y(I)-R(I))-A(K5PI))*A(K4PI)
 20     TEMD1 = TEMD1 + TEMD2
      TEMD1 = TEMD1/W11
      A(K2PJ+1) = TEMD1
C
C UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
C ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
C COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
C THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
C SIGNIFICANT BITS ARE IN  A(K5PI) .
C
      SIGJ = 0.0D0
      DO 21 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = R(I) + A(K5PI) + TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
 21     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
C
C SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
C MAXDEG  HAS BEEN REACHED
C
      IF (EPS) 23,26,27
C
C COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
C
 23   IF (SIGJ .EQ. 0.0D0) GO TO 29
      DEGF = M - J - 1
      DEN = (CO(4,KSIG)*DEGF + 1.0D0)*DEGF
      FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
      FCRIT = FCRIT*FCRIT
      F = (SIGJM1 - SIGJ)*DEGF/SIGJ
      IF (F .LT. FCRIT) GO TO 25
C
C POLYNOMIAL OF DEGREE J SATISFIES F TEST
C
 24   SIGPAS = SIGJ
      JPAS = J
      NFAIL = 0
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
C SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
C
 25   NFAIL = NFAIL + 1
      IF (NFAIL .GE. 3) GO TO 29
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
C EPS = 0.)
C
 26   IF (MAXDEG .EQ. J) GO TO 28
      GO TO 16
C
C SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
C
 27   IF (SIGJ .LE. ETST) GO TO 28
      IF (MAXDEG .EQ. J) GO TO 31
      GO TO 16
C
C RETURNS
C
 28   IERR = 1
      NDEG = J
      SIG = SIGJ
      GO TO 33
 29   IERR = 1
      NDEG = JPAS
      SIG = SIGPAS
      GO TO 33
 30   IERR = 2
      CALL XERROR ( 'DPOLFT-INVALID INPUT PARAMETER.',31,2,1)
      GO TO 37
 31   IERR = 3
      NDEG = MAXDEG
      SIG = SIGJ
      GO TO 33
 32   IERR = 4
      NDEG = JPAS
      SIG = SIGPAS
C
 33   A(K3) = NDEG
C
C WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
C ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
C
      IF(EPS .GE. 0.0  .OR.  NDEG .EQ. MAXDEG) GO TO 36
      NDER = 0
      DO 35 I = 1,M
        CALL DP1VLU (NDEG,NDER,X(I),R(I),YP,A)
 35     CONTINUE
 36   EPS = DSQRT(SIG/XM)
 37   RETURN
      END
*DECK DROTM
      SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
C***BEGIN PROLOGUE  DROTM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A8
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SROTM-S DROTM-D),LINEAR ALGEBRA,
C             MODIFIED GIVENS ROTATION,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Apply d.p. modified Givens transformation
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below.
C            Elements 2-5 form the transformation matrix H.
C
C     --Output--
C       DX  rotated vector (unchanged if N .LE. 0)
C       DY  rotated vector (unchanged if N .LE. 0)
C
C     Apply the modified Givens transformation, H, to the 2 by N matrix
C
C     (DX**T), where **T indicates transpose.  The elements of DX are in
C     (DY**T)
C
C     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = (-INCX)*N, and similarly for SY using LY and INCY.
C     With DPARAM(1)=DFLAG, H has one of the following forms.
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C
C     See DROTMG for a description of data storage in DPARAM.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DROTM
C
      DOUBLE PRECISION DFLAG,DH12,DH22,DX,TWO,Z,DH11,DH21,
     1 DPARAM,DY,W,ZERO
      DIMENSION DX(*),DY(*),DPARAM(5)
      SAVE ZERO, TWO
      DATA ZERO,TWO/0.D0,2.D0/
C***FIRST EXECUTABLE STATEMENT  DROTM
      DFLAG=DPARAM(1)
      IF(N .LE. 0 .OR.(DFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF(DFLAG) 50,10,30
   10          CONTINUE
               DH12=DPARAM(4)
               DH21=DPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W+Z*DH12
                    DY(I)=W*DH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               DH11=DPARAM(2)
               DH22=DPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z
                    DY(I)=-W+DH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               DH11=DPARAM(2)
               DH12=DPARAM(4)
               DH21=DPARAM(3)
               DH22=DPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z*DH12
                    DY(I)=W*DH21+Z*DH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
C
          IF(DFLAG)120,80,100
   80     CONTINUE
          DH12=DPARAM(4)
          DH21=DPARAM(3)
               DO 90 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W+Z*DH12
               DY(KY)=W*DH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          DH11=DPARAM(2)
          DH22=DPARAM(5)
               DO 110 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z
               DY(KY)=-W+DH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          DH11=DPARAM(2)
          DH12=DPARAM(4)
          DH21=DPARAM(3)
          DH22=DPARAM(5)
               DO 130 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z*DH12
               DY(KY)=W*DH21+Z*DH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
      END
*DECK DROTMG
      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
C***BEGIN PROLOGUE  DROTMG
C***DATE WRITTEN   780301   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1B10
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SROTMG-S DROTMG-D),LINEAR ALGEBRA,
C             MODIFIED GIVENS ROTATION,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Construct d.p. modified Givens transformation
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C      DD1  double precision scalar
C      DD2  double precision scalar
C      DX1  double precision scalar
C      DX2  double precision scalar
C   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below.
C             Elements 2-5  define the transformation matrix H.
C
C     --Output--
C      DD1  changed to represent the effect of the transformation
C      DD2  changed to reflect the transformation
C      DX1  changed to reflect the transformation
C      DX2  unchanged
C
C     Construct the modified Givens transformation matrix H which zeros
C     the second component of the 2-vector  (DSQRT(DD1)*DX1,DSQRT(DD2)*
C     DY2)**T.
C     With DPARAM(1)=DFLAG, H has one of the following forms..
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C
C     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22
C     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the
C     value of DPARAM(1) are not stored in DPARAM.)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DROTMG
C
      DOUBLE PRECISION GAM,ONE,RGAMSQ,DD2,DH11,DH21,DPARAM,DP2,
     1 DQ2,DU,DY1,ZERO,GAMSQ,DD1,DFLAG,DH12,DH22,DP1,DQ1,
     2 DTEMP,DX1,TWO
      DIMENSION DPARAM(5)
      SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
      DATA ZERO,ONE,TWO /0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
C***FIRST EXECUTABLE STATEMENT  DROTMG
      IF(.NOT. DD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-DX1..
          GO TO 60
   10 CONTINUE
C     CASE-DD1-NONNEGATIVE
      DP2=DD2*DY1
      IF(.NOT. DP2 .EQ. ZERO) GO TO 20
          DFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      DP1=DD1*DX1
      DQ2=DP2*DY1
      DQ1=DP1*DX1
C
      IF(.NOT. DABS(DQ1) .GT. DABS(DQ2)) GO TO 40
          DH21=-DY1/DX1
          DH12=DP2/DP1
C
          DU=ONE-DH12*DH21
C
          IF(.NOT. DU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   30     CONTINUE
               DFLAG=ZERO
               DD1=DD1/DU
               DD2=DD2/DU
               DX1=DX1*DU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. DQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   50     CONTINUE
               DFLAG=ONE
               DH11=DP1/DP2
               DH22=DX1/DY1
               DU=ONE+DH11*DH22
               DTEMP=DD2/DU
               DD2=DD1/DU
               DD1=DTEMP
               DX1=DY1*DU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
          DFLAG=-ONE
          DH11=ZERO
          DH12=ZERO
          DH21=ZERO
          DH22=ZERO
C
          DD1=ZERO
          DD2=ZERO
          DX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. DFLAG .GE. ZERO) GO TO 90
C
          IF(.NOT. DFLAG .EQ. ZERO) GO TO 80
          DH11=ONE
          DH22=ONE
          DFLAG=-ONE
          GO TO 90
   80     CONTINUE
          DH21=-ONE
          DH12=ONE
          DFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. DD1 .LE. RGAMSQ) GO TO 130
               IF(DD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               DD1=DD1*GAM**2
               DX1=DX1/GAM
               DH11=DH11/GAM
               DH12=DH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. DD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               DD1=DD1/GAM**2
               DX1=DX1*GAM
               DH11=DH11*GAM
               DH12=DH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. DABS(DD2) .LE. RGAMSQ) GO TO 190
               IF(DD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               DD2=DD2*GAM**2
               DH21=DH21/GAM
               DH22=DH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. DABS(DD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               DD2=DD2/GAM**2
               DH21=DH21*GAM
               DH22=DH22*GAM
          GO TO 200
  220 CONTINUE
          IF(DFLAG)250,230,240
  230     CONTINUE
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               GO TO 260
  240     CONTINUE
               DPARAM(2)=DH11
               DPARAM(5)=DH22
               GO TO 260
  250     CONTINUE
               DPARAM(2)=DH11
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               DPARAM(5)=DH22
  260 CONTINUE
          DPARAM(1)=DFLAG
          RETURN
      END
*DECK DSWAP
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DSWAP
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SSWAP-S DSWAP-D CSWAP-C ISWAP-I),
C             INTERCHANGE,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Interchange d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSWAP
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP1,DTEMP2,DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
*DECK DSCAL
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SSCAL-S DSCAL-D CSCAL-C),
C             LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSCAL
C
      DOUBLE PRECISION DA,DX(1)
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
*DECK DVNRMS
      DOUBLE PRECISION FUNCTION DVNRMS(N,V,W)
C***BEGIN PROLOGUE  DVNRMS
C***REFER TO  DDEBDF
C
C   DVNRMS computes a weighted root-mean-square vector norm for the
C   integrator package DDEBDF.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DVNRMS
      INTEGER I, N
      DOUBLE PRECISION DBLE, DSQRT, SUM, V, W
      DIMENSION V(N),W(N)
C***FIRST EXECUTABLE STATEMENT  DVNRMS
      SUM = 0.0D0
      DO 10 I = 1, N
         SUM = SUM + (V(I)/W(I))**2
   10 CONTINUE
      DVNRMS = DSQRT(SUM/DBLE(FLOAT(N)))
      RETURN
C     ----------------------- END OF FUNCTION DVNRMS
C     ------------------------
      END
*DECK DWNLIT
      SUBROUTINE DWNLIT(W,MDW,M,N,L,IPIVOT,ITYPE,H,SCALE,RNORM,IDOPE,
     *   DOPE,DONE)
C***BEGIN PROLOGUE  DWNLIT
C***REFER TO  DLSEI,DWNNLS
C***ROUTINES CALLED  DCOPY,DH12,DROTM,DROTMG,DSWAP,IDAMAX
C***DESCRIPTION
C
C     THIS IS A COMPANION SUBPROGRAM TO DWNNLS( ).
C     THE DOCUMENTATION FOR DWNNLS( ) HAS MORE COMPLETE
C     USAGE INSTRUCTIONS.
C
C     NOTE  THE M BY (N+1) MATRIX W( , ) CONTAINS THE RT. HAND SIDE
C           B AS THE (N+1)ST COL.
C
C
C     TRIANGULARIZE L1 BY L1 SUBSYSTEM, WHERE L1=MIN(M,L), WITH
C     COL INTERCHANGES.
C      REVISED APRIL 19, 1981.
C***END PROLOGUE  DWNLIT
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /DOUBLE PRECISION (12 BLANKS)/DOUBLE
C     PRECISION/,/DCOPY/DCOPY/,/DROTM/DROTM/,
C     /DSWAP/DSWAP/,/DMAX1/DMAX1/,/IDAMAX/IDAMAX/,/.E-/.D-/,/E0/D0/
C
      INTEGER I, I1, IC, IDAMAX, IDOPE(8), IGO990, IGO993, IGO996,
     *     IMAX, IP1, IPIVOT(*), IR, IRP1, ITEMP, ITYPE(*), J, J1, JJ,
     *     JM1, JP, K, KK, KRANK, KRP1, L, L1, LB, LEND, LP1, M, MAX,
     *     MDW, ME, MEND, MEP1, MIN0, N, NIV, NIV1, NP1, NSOLN
      DOUBLE PRECISION ALSQ, AMAX, DMAX1, DOPE(4), EPENSQ, FAC, FACTOR,
     *     H(*), HBAR, ONE, RN, RNORM, SCALE(*), SN, SPARAM(5), T, TAU,
     *     TENM3, W(MDW,*), ZERO
      LOGICAL INDEP,DONE,RECALC
      DATA TENM3 /1.0D-3/, ZERO /0.0D0/, ONE /1.0D0/
C
C***FIRST EXECUTABLE STATEMENT  DWNLIT
C     BEGIN BLOCK PERMITTING ...EXITS TO 530
C        BEGIN BLOCK PERMITTING ...EXITS TO 510
C           BEGIN BLOCK PERMITTING ...EXITS TO 470
C              BEGIN BLOCK PERMITTING ...EXITS TO 410
C                 BEGIN BLOCK PERMITTING ...EXITS TO 380
C                    BEGIN BLOCK PERMITTING ...EXITS TO 370
C                       BEGIN BLOCK PERMITTING ...EXITS TO 250
C                          BEGIN BLOCK PERMITTING ...EXITS TO 240
C                             BEGIN BLOCK PERMITTING ...EXITS TO 140
C                                BEGIN BLOCK PERMITTING ...EXITS TO 120
                                    ME = IDOPE(1)
                                    MEP1 = IDOPE(2)
                                    KRANK = IDOPE(3)
                                    KRP1 = IDOPE(4)
                                    NSOLN = IDOPE(5)
                                    NIV = IDOPE(6)
                                    NIV1 = IDOPE(7)
                                    L1 = IDOPE(8)
C
                                    ALSQ = DOPE(1)
                                    EPENSQ = DOPE(2)
                                    FAC = DOPE(3)
                                    TAU = DOPE(4)
                                    NP1 = N + 1
                                    LB = MIN0(M-1,L)
                                    RECALC = .TRUE.
                                    RNORM = ZERO
                                    KRANK = 0
C                                   WE SET FACTOR=1.E0 SO THAT THE HEAVY
C                                   WEIGHT ALAMDA WILL BE INCLUDED IN
C                                   THE TEST FOR COL INDEPENDENCE.
                                    FACTOR = 1.0D0
                                    I = 1
                                    IP1 = 2
C                                ...EXIT
                                    GO TO 120
C
C                                   UPDATE-COL-SS-AND-FIND-PIVOT-COL
   10                               CONTINUE
                                    ASSIGN 20 TO IGO993
C        ...........................EXIT
                                    GO TO 510
C
C                                   PERFORM-COL-INTERCHANGE
C
C                                   SET IC TO POINT TO I-TH COL.
   20                               CONTINUE
                                    IC = I
                                    ASSIGN 30 TO IGO990
C           ........................EXIT
                                    GO TO 470
C
C                                   TEST-INDEP-OF-INCOMING-COL
   30                               CONTINUE
                                    IF (.NOT.INDEP) GO TO 100
C
C                                      ELIMINATE I-TH COL BELOW DIAG.
C                                      USING MOD. GIVENS TRANSFORMATIONS
C                                      APPLIED TO (A B).
                                       J = M
                                       DO 90 JJ = IP1, M
                                          JM1 = J - 1
                                          JP = JM1
C                                         WHEN OPERATING NEAR THE ME
C                                         LINE, USE THE LARGEST ELT.
C                                         ABOVE IT AS THE PIVOT.
                                          IF (J .NE. MEP1) GO TO 70
                                             IMAX = ME
                                             AMAX = SCALE(ME)*W(ME,I)**2
   40                                        IF (JP .LT. I) GO TO 60
                                                T = SCALE(JP)*W(JP,I)**2
                                                IF (T .LE. AMAX)
     *                                             GO TO 50
                                                   IMAX = JP
                                                   AMAX = T
   50                                           CONTINUE
                                                JP = JP - 1
                                             GO TO 40
   60                                        CONTINUE
                                             JP = IMAX
   70                                     CONTINUE
                                          IF (W(J,I) .EQ. ZERO) GO TO 80
                                             CALL DROTMG(SCALE(JP),
     *                                                   SCALE(J),
     *                                                   W(JP,I),W(J,I),
     *                                                   SPARAM)
                                             W(J,I) = ZERO
                                             CALL DROTM(NP1-I,W(JP,IP1),
     *                                                  MDW,W(J,IP1),
     *                                                  MDW,SPARAM)
   80                                     CONTINUE
                                          J = JM1
   90                                  CONTINUE
                                       KRANK = KRANK + 1
                                    GO TO 110
  100                               CONTINUE
C                             .........EXIT
                                       GO TO 140
  110                               CONTINUE
                                    I = IP1
                                    IP1 = IP1 + 1
  120                            CONTINUE
                                 IF (I .GT. LB) GO TO 130
C
C                                   SET IR TO POINT TO THE I-TH ROW.
                                    IR = I
                                    LEND = L
                                    MEND = M
                                    ASSIGN 10 TO IGO996
C              .....................EXIT
                                    GO TO 410
  130                            CONTINUE
                                 KRANK = L1
  140                         CONTINUE
                              KRP1 = KRANK + 1
C                       ......EXIT
                              IF (KRANK .GE. ME) GO TO 250
                              FACTOR = ALSQ
C
C                             DETERMINE THE RANK OF THE REMAINING
C                             EQUALITY CONSTRAINT EQUATIONS.  REMOVE ANY
C                             REDUNDANT EQUATIONS.
                              LP1 = L + 1
                              RECALC = .TRUE.
                              LB = MIN0(L+ME-KRANK,N)
                              I = LP1
                              IP1 = I + 1
C                          ...EXIT
                              GO TO 240
C
C                             UPDATE-COL-SS-AND-FIND-PIVOT-COL
  150                         CONTINUE
                              ASSIGN 160 TO IGO993
C        .....................EXIT
                              GO TO 510
C
C                             PERFORM-COL-INTERCHANGE
C
C                             ELIMINATE ELEMENTS IN THE I-TH COL.
  160                         CONTINUE
                              J = ME
  170                         IF (J .LE. IR) GO TO 190
                                 JM1 = J - 1
                                 IF (W(J,I) .EQ. ZERO) GO TO 180
                                    CALL DROTMG(SCALE(JM1),SCALE(J),
     *                                          W(JM1,I),W(J,I),SPARAM)
                                    W(J,I) = ZERO
                                    CALL DROTM(NP1-I,W(JM1,IP1),MDW,
     *                                         W(J,IP1),MDW,SPARAM)
  180                            CONTINUE
                                 J = JM1
                              GO TO 170
  190                         CONTINUE
C
C                             SET IC=I=COL BEING ELIMINATED
                              IC = I
                              ASSIGN 200 TO IGO990
C           ..................EXIT
                              GO TO 470
C
C                             TEST-INDEP-OF-INCOMING-COL
  200                         CONTINUE
                              IF (INDEP) GO TO 230
C
C                                REMOVE ANY REDUNDANT OR DEPENDENT
C                                EQUALITY CONSTRAINTS.
                                 JJ = IR
  210                            IF (IR .GT. ME) GO TO 220
                                    W(IR,1) = ZERO
                                    CALL DCOPY(N,W(IR,1),0,W(IR,1),MDW)
                                    RNORM = RNORM
     *                                      + SCALE(IR)*W(IR,NP1)**2
     *                                        /ALSQ
                                    W(IR,NP1) = ZERO
                                    SCALE(IR) = ONE
                                    IR = IR + 1
                                 GO TO 210
  220                            CONTINUE
C
C                                REDUCE ME TO REFLECT ANY DISCOVERED
C                                DEPENDENT EQUALITY CONSTRAINTS.
                                 ME = JJ - 1
                                 MEP1 = ME + 1
C                       .........EXIT
                                 GO TO 250
  230                         CONTINUE
                              I = IP1
                              IP1 = IP1 + 1
  240                      CONTINUE
C                       ...EXIT
                           IF (I .GT. LB) GO TO 250
                           IR = KRANK + I - L
                           LEND = N
                           MEND = ME
                           ASSIGN 150 TO IGO996
C              ............EXIT
                           GO TO 410
  250                   CONTINUE
C                 ......EXIT
                        IF (KRANK .GE. L1) GO TO 380
C
C                       TRY TO DETERMINE THE VARIABLES KRANK+1 THROUGH
C                       L1 FROM THE LEAST SQUARES EQUATIONS.  CONTINUE
C                       THE TRIANGULARIZATION WITH PIVOT ELEMENT
C                       W(MEP1,I).
C
                        RECALC = .TRUE.
C
C                       SET FACTOR=ALSQ TO REMOVE EFFECT OF HEAVY WEIGHT
C                       FROM TEST FOR COL INDEPENDENCE.
                        FACTOR = ALSQ
                        KK = KRP1
                        I = KK
                        IP1 = I + 1
C                    ...EXIT
                        GO TO 370
C
C                       UPDATE-COL-SS-AND-FIND-PIVOT-COL
  260                   CONTINUE
                        ASSIGN 270 TO IGO993
C        ...............EXIT
                        GO TO 510
C
C                       PERFORM-COL-INTERCHANGE
C
C                       ELIMINATE I-TH COL BELOW THE IR-TH ELEMENT.
  270                   CONTINUE
                        IRP1 = IR + 1
                        J = M
                        DO 290 JJ = IRP1, M
                           JM1 = J - 1
                           IF (W(J,I) .EQ. ZERO) GO TO 280
                              CALL DROTMG(SCALE(JM1),SCALE(J),W(JM1,I),
     *                                    W(J,I),SPARAM)
                              W(J,I) = ZERO
                              CALL DROTM(NP1-I,W(JM1,IP1),MDW,W(J,IP1),
     *                                   MDW,SPARAM)
  280                      CONTINUE
                           J = JM1
  290                   CONTINUE
C
C                       TEST IF NEW PIVOT ELEMENT IS NONZERO. IF NOT,
C                       THE COL IS DEPENDENT.
                        T = SCALE(IR)*W(IR,I)**2
                        IF (FAC*T + EPENSQ .NE. EPENSQ) GO TO 300
                           INDEP = .FALSE.
                        GO TO 330
  300                   CONTINUE
C
C                          COL TEST PASSED. NOW MUST PASS ROW NORM TEST
C                          TO BE CLASSIFIED AS INDEPENDENT.
                           RN = ZERO
                           DO 320 I1 = IR, M
                              DO 310 J1 = IP1, N
                                 RN = DMAX1(RN,SCALE(I1)*W(I1,J1)**2)
  310                         CONTINUE
  320                      CONTINUE
                           INDEP = FAC*T + RN .NE. RN
  330                   CONTINUE
C
C                       IF INDEPENDENT, SWAP THE IR-TH AND KRP1-ST ROWS
C                       TO MAINTAIN THE TRIANGULAR FORM.  UPDATE THE
C                       RANK INDICATOR KRANK AND THE EQUALITY CONSTRAINT
C                       POINTER ME.
                        IF (.NOT.INDEP) GO TO 350
                           IF (IR .EQ. KRP1) GO TO 340
                              CALL DSWAP(NP1,W(KRP1,1),MDW,W(IR,1),MDW)
                              CALL DSWAP(1,SCALE(KRP1),1,SCALE(IR),1)
                              ITEMP = ITYPE(KRP1)
                              ITYPE(KRP1) = ITYPE(IR)
                              ITYPE(IR) = ITEMP
  340                      CONTINUE
                           ME = MEP1
                           MEP1 = ME + 1
                           KRANK = KRP1
                           KRP1 = KRANK + 1
                        GO TO 360
  350                   CONTINUE
C                 .........EXIT
                           GO TO 380
  360                   CONTINUE
                        I = IP1
                        IP1 = IP1 + 1
  370                CONTINUE
C                 ...EXIT
                     IF (I .GT. L1) GO TO 380
C
C                    SET IR TO POINT TO THE MEP1-ST ROW.
                     IR = MEP1
                     LEND = L
                     MEND = M
                     ASSIGN 260 TO IGO996
C              ......EXIT
                     GO TO 410
  380             CONTINUE
C
C                 IF PSEUDORANK IS LESS THAN L, APPLY HOUSEHOLDER TRANS.
C                 FROM RIGHT.
                  IF (KRANK .GE. L) GO TO 400
                     DO 390 I = 1, KRANK
                        J = KRP1 - I
                        CALL DH12(1,J,KRP1,L,W(J,1),MDW,H(J),W,MDW,1,
     *                            J-1)
  390                CONTINUE
  400             CONTINUE
                  NIV = KRANK + NSOLN - L
                  NIV1 = NIV + 1
                  IF (L .EQ. N) DONE = .TRUE.
C
C                 END OF INITIAL TRIANGULARIZATION.
                  IDOPE(1) = ME
                  IDOPE(2) = MEP1
                  IDOPE(3) = KRANK
                  IDOPE(4) = KRP1
                  IDOPE(5) = NSOLN
                  IDOPE(6) = NIV
                  IDOPE(7) = NIV1
                  IDOPE(8) = L1
C     ............EXIT
                  GO TO 530
  410          CONTINUE
C
C              TO UPDATE-COL-SS-AND-FIND-PIVOT-COL
C
C              THE COL SS VECTOR WILL BE UPDATED AT EACH STEP. WHEN
C              NUMERICALLY NECESSARY, THESE VALUES WILL BE RECOMPUTED.
C
               IF (IR .EQ. 1 .OR. RECALC) GO TO 430
C                 UPDATE COL SS =SUM OF SQUARES.
                  DO 420 J = I, LEND
                     H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
  420             CONTINUE
C
C                 TEST FOR NUMERICAL ACCURACY.
                  MAX = IDAMAX(LEND-I+1,H(I),1) + I - 1
                  RECALC = HBAR + TENM3*H(MAX) .EQ. HBAR
  430          CONTINUE
C
C              IF REQUIRED, RECALCULATE COL SS, USING ROWS IR THROUGH
C              MEND.
               IF (RECALC) GO TO 440
                  GO TO IGO996, (10,150,260)
C     ............EXIT
                  GO TO 530
  440          CONTINUE
               DO 460 J = I, LEND
                  H(J) = ZERO
                  DO 450 K = IR, MEND
                     H(J) = H(J) + SCALE(K)*W(K,J)**2
  450             CONTINUE
  460          CONTINUE
C
C              FIND COL WITH LARGEST SS.
               MAX = IDAMAX(LEND-I+1,H(I),1) + I - 1
               HBAR = H(MAX)
               GO TO IGO996, (10,150,260)
C     .........EXIT
               GO TO 530
  470       CONTINUE
C
C           TO TEST-INDEP-OF-INCOMING-COL
C
C           TEST THE COL IC TO DETERMINE IF IT IS LINEARLY INDEPENDENT
C           OF THE COLS ALREADY IN THE BASIS.  IN THE INIT TRI
C           STEP, WE USUALLY WANT THE HEAVY WEIGHT ALAMDA TO
C           BE INCLUDED IN THE TEST FOR INDEPENDENCE.  IN THIS CASE THE
C           VALUE OF FACTOR WILL HAVE BEEN SET TO 1.E0 BEFORE THIS
C           PROCEDURE IS INVOKED.  IN THE POTENTIALLY RANK DEFICIENT
C           PROBLEM, THE VALUE OF FACTOR WILL HAVE BEEN
C           SET TO ALSQ=ALAMDA**2 TO REMOVE THE EFFECT OF THE HEAVY
C           WEIGHT FROM THE TEST FOR INDEPENDENCE.
C
C           WRITE NEW COL AS PARTITIONED VECTOR
C                   (A1)  NUMBER OF COMPONENTS IN SOLN SO FAR = NIV
C                   (A2)  M-NIV COMPONENTS
C           AND COMPUTE  SN = INVERSE WEIGHTED LENGTH OF A1
C                        RN = INVERSE WEIGHTED LENGTH OF A2
C           CALL THE COL INDEPENDENT WHEN RN .GT. TAU*SN
            SN = ZERO
            RN = ZERO
            DO 500 J = 1, MEND
               T = SCALE(J)
               IF (J .LE. ME) T = T/FACTOR
               T = T*W(J,IC)**2
               IF (J .GE. IR) GO TO 480
                  SN = SN + T
               GO TO 490
  480          CONTINUE
                  RN = RN + T
  490          CONTINUE
  500       CONTINUE
            INDEP = RN .GT. TAU**2*SN
            GO TO IGO990, (30,200)
            GO TO IGO996, (10,150,260)
C     ......EXIT
            GO TO 530
  510    CONTINUE
C
C        TO PERFORM-COL-INTERCHANGE
C
         IF (MAX .EQ. I) GO TO 520
C           EXCHANGE ELEMENTS OF PERMUTED INDEX VECTOR AND PERFORM COL
C           INTERCHANGES.
            ITEMP = IPIVOT(I)
            IPIVOT(I) = IPIVOT(MAX)
            IPIVOT(MAX) = ITEMP
            CALL DSWAP(M,W(1,MAX),1,W(1,I),1)
            H(MAX) = H(I)
  520    CONTINUE
         GO TO IGO993, (20,160,270)
         GO TO IGO990, (30,200)
         GO TO IGO996, (10,150,260)
  530 CONTINUE
      RETURN
      END
*DECK DWNLSM
      SUBROUTINE DWNLSM(W,MDW,MME,MA,N,L,PRGOPT,X,RNORM,MODE,IPIVOT,
     *   ITYPE,WD,H,SCALE,Z,TEMP,D)
C***BEGIN PROLOGUE  DWNLSM
C***REFER TO  DWNNLS
C***ROUTINES CALLED  D1MACH,DASUM,DAXPY,DCOPY,DH12,DNRM2,DROTM,DROTMG,
C                    DSCAL,DSWAP,DWNLIT,IDAMAX,XERROR
C***DESCRIPTION
C
C   **** Double Precision version of WNLSM ****
C
C     This is a companion subprogram to DWNNLS( ).
C     The documentation for DWNNLS( ) has more complete
C     usage instructions.
C
C     Written by Karen H. Haskell, Sandia Laboratories,
C     with the help of R.J. Hanson, Sandia Laboratories,
C     December 1976 - January 1978.
C     Revised March 4, 1982.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine DWNNLS, the following work arrays are used in
C     subroutine DWNLSM  (they are passed through the calling
C     sequence from DWNNLS for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C***END PROLOGUE  DWNLSM
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (BEGIN CHANGES AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SROTMG/DROTMG/,
C     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/SROTM/DROTM/,/AMAX1/DMAX1/,
C     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/E0/D0/,/SSWAP/DSWAP/,
C     /ISAMAX/IDAMAX/,/SRELPR/DRELPR/
C
C     SUBROUTINE DWNLSM (W,MDW,MME,MA,N,L,PRGOPT,X,RNORM,MODE,
C    1                  IPIVOT,ITYPE,WD,H,SCALE,Z,TEMP,D)
C++
      DOUBLE PRECISION W(MDW,1), X(*), WD(*), H(*), SCALE(*), DOPE(4)
      DOUBLE PRECISION Z(*), TEMP(*), PRGOPT(*), D(*), SPARAM(5)
      DOUBLE PRECISION ALAMDA, ALPHA, ALSQ, AMAX, BNORM, EANORM
      DOUBLE PRECISION DRELPR, FAC, ONE, BLOWUP
      DOUBLE PRECISION RNORM, SM, T, TAU, TWO, WMAX, ZERO, ZZ, Z2
      DOUBLE PRECISION DMAX1, DSQRT, DNRM2, DASUM, DMIN1, D1MACH
      INTEGER IPIVOT(*), ITYPE(*), IDAMAX, IDOPE(8)
      LOGICAL HITCON, FEASBL, DONE, POS
      DATA ZERO /0.D0/, ONE /1.D0/, TWO /2.D0/, DRELPR /0.D0/
C
C     INITIALIZE-VARIABLES
C***FIRST EXECUTABLE STATEMENT  DWNLSM
      ASSIGN 10 TO IGO998
      GO TO 180
C
C     PERFORM INITIAL TRIANGULARIZATION IN THE SUBMATRIX
C     CORRESPONDING TO THE UNCONSTRAINED VARIABLES USING
C     THE PROCEDURE INITIALLY-TRIANGULARIZE.
   10 ASSIGN 20 TO IGO995
      GO TO 280
C
C     PERFORM DWNNLS ALGORITHM USING THE FOLLOWING STEPS.
C
C     UNTIL(DONE)
C
C        COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
C        WHEN (HITCON) ADD-CONSTRAINTS
C
C        ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
C
C        FIN
C
C     COMPUTE-FINAL-SOLUTION
C
   20 IF (DONE) GO TO 80
C
      ASSIGN 30 TO IGO991
      GO TO 300
C
C     COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
   30 IF (.NOT.(HITCON)) GO TO 50
      ASSIGN 40 TO IGO986
      GO TO 370
   40 GO TO 70
C
C     WHEN (HITCON) ADD-CONSTRAINTS
C
   50 ASSIGN 60 TO IGO983
      GO TO 640
   60 CONTINUE
C
C     ELSE PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
C
   70 GO TO 20
C
   80 ASSIGN 90 TO IGO980
      GO TO 1000
C
C     COMPUTE-FINAL-SOLUTION
C
   90 RETURN
  100 CONTINUE
C
C     TO PROCESS-OPTION-VECTOR
      FAC = 1.E-4
C
C     THE NOMINAL TOLERANCE USED IN THE CODE,
      TAU = DSQRT(DRELPR)
C
C     THE NOMINAL BLOW-UP FACTOR USED IN THE CODE.
      BLOWUP = TAU
C
C     THE NOMINAL COLUMN SCALING USED IN THE CODE IS
C     THE IDENTITY SCALING.
      D(1) = ONE
      CALL DCOPY(N, D, 0, D, 1)
C
C     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE.
      NOPT = 1000
C
C     DEFINE BOUND FOR POSITIVE VALUE OF LINK.
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (.NOT.(LINK.LE.0 .OR. LINK.GT.NLINK)) GO TO 110
      NERR = 3
      IOPT = 1
      CALL XERROR( 'DWNNLS( ) THE OPTION VECTOR IS UNDEFINED', 39, NERR,
     1 IOPT)
      MODE = 2
      RETURN
  110 IF (.NOT.(LINK.GT.1)) GO TO 160
      NTIMES = NTIMES + 1
      IF (.NOT.(NTIMES.GT.NOPT)) GO TO 120
      NERR = 3
      IOPT = 1
      CALL XERROR('DWNNLS( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING
     1.', 53,     NERR, IOPT)
      MODE = 2
      RETURN
  120 KEY = PRGOPT(LAST+1)
      IF (.NOT.(KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.ZERO)) GO TO 140
      DO 130 J=1,N
        T = DNRM2(M,W(1,J),1)
        IF (T.NE.ZERO) T = ONE/T
        D(J) = T
  130 CONTINUE
  140 IF (KEY.EQ.7) CALL DCOPY(N, PRGOPT(LAST+2), 1, D, 1)
      IF (KEY.EQ.8) TAU = DMAX1(DRELPR,PRGOPT(LAST+2))
      IF (KEY.EQ.9) BLOWUP = DMAX1(DRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.NLINK)) GO TO 150
      NERR = 3
      IOPT = 1
      CALL XERROR( 'DWNNLS( ) THE OPTION VECTOR IS UNDEFINED', 39, NERR,
     1 IOPT)
      MODE = 2
      RETURN
  150 LAST = LINK
      LINK = NEXT
      GO TO 110
  160 DO 170 J=1,N
        CALL DSCAL(M, D(J), W(1,J), 1)
  170 CONTINUE
      GO TO 1260
  180 CONTINUE
C
C     TO INITIALIZE-VARIABLES
C
C     DRELPR IS THE PRECISION FOR THE PARTICULAR MACHINE
C     BEING USED.  THIS LOGIC AVOIDS RECOMPUTING IT EVERY ENTRY.
      IF (.NOT.(DRELPR.EQ.ZERO)) GO TO 210
      DRELPR = ONE
  190 IF (ONE+DRELPR.EQ.ONE) GO TO 200
      DRELPR = DRELPR/TWO
      GO TO 190
  200 DRELPR = DRELPR*TWO
  210 M = MA + MME
      ME = MME
      MEP1 = ME + 1
      ASSIGN 220 TO IGO977
      GO TO 100
C
C     PROCESS-OPTION-VECTOR
  220 DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      LP1 = L + 1
      NSOLN = L
      NSP1 = NSOLN + 1
      NP1 = N + 1
      NM1 = N - 1
      L1 = MIN0(M,L)
C
C     COMPUTE SCALE FACTOR TO APPLY TO EQUAL. CONSTRAINT EQUAS.
      DO 230 J=1,N
        WD(J) = DASUM(M,W(1,J),1)
  230 CONTINUE
      IMAX = IDAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = DASUM(M,W(1,NP1),1)
      ALAMDA = EANORM/(DRELPR*FAC)
C     TO KEEP FROM OVERFLOWING ON A VAX 11/780
      ALAMDA = DMIN1(ALAMDA,DSQRT(D1MACH(2)))
C
C     DEFINE SCALING DIAG MATRIX FOR MOD GIVENS USAGE AND
C     CLASSIFY EQUATION TYPES.
      ALSQ = ALAMDA**2
      DO 260 I=1,M
C
C     WHEN EQU I IS HEAVILY WEIGHTED ITYPE(I)=0, ELSE ITYPE(I)=1.
        IF (.NOT.(I.LE.ME)) GO TO 240
        T = ALSQ
        ITEMP = 0
        GO TO 250
  240   T = ONE
        ITEMP = 1
  250   SCALE(I) = T
        ITYPE(I) = ITEMP
  260 CONTINUE
C
C     SET THE SOLN VECTOR X(*) TO ZERO AND THE COL INTERCHANGE
C     MATRIX TO THE IDENTITY.
      X(1) = ZERO
      CALL DCOPY(N, X, 0, X, 1)
      DO 270 I=1,N
        IPIVOT(I) = I
  270 CONTINUE
      GO TO 1230
  280 CONTINUE
C
C     TO INITIALLY-TRIANGULARIZE
C
C     SET FIRST L COMPS. OF DUAL VECTOR TO ZERO BECAUSE
C     THESE CORRESPOND TO THE UNCONSTRAINED VARIABLES.
      IF (.NOT.(L.GT.0)) GO TO 290
      WD(1) = ZERO
      CALL DCOPY(L, WD, 0, WD, 1)
C
C     THE ARRAYS IDOPE(*) AND DOPE(*) ARE USED TO PASS
C     INFORMATION TO DWNLIT().  THIS WAS DONE TO AVOID
C     A LONG CALLING SEQUENCE OR THE USE OF COMMON.
  290 IDOPE(1) = ME
      IDOPE(2) = MEP1
      IDOPE(3) = 0
      IDOPE(4) = 1
      IDOPE(5) = NSOLN
      IDOPE(6) = 0
      IDOPE(7) = 1
      IDOPE(8) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = FAC
      DOPE(4) = TAU
      CALL DWNLIT(W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     1 IDOPE, DOPE, DONE)
      ME = IDOPE(1)
      MEP1 = IDOPE(2)
      KRANK = IDOPE(3)
      KRP1 = IDOPE(4)
      NSOLN = IDOPE(5)
      NIV = IDOPE(6)
      NIV1 = IDOPE(7)
      L1 = IDOPE(8)
      GO TO 1240
  300 CONTINUE
C
C     TO COMPUTE-SEARCH-DIRECTION-AND-FEASIBLE-POINT
C
C     SOLVE THE TRIANGULAR SYSTEM OF CURRENTLY NON-ACTIVE
C     VARIABLES AND STORE THE SOLUTION IN Z(*).
C
C     SOLVE-SYSTEM
      ASSIGN 310 TO IGO958
      GO TO 1110
C
C     INCREMENT ITERATION COUNTER AND CHECK AGAINST MAX. NUMBER
C     OF ITERATIONS.
  310 ITER = ITER + 1
      IF (.NOT.(ITER.GT.ITMAX)) GO TO 320
      MODE = 1
      DONE = .TRUE.
C
C     CHECK TO SEE IF ANY CONSTRAINTS HAVE BECOME ACTIVE.
C     IF SO, CALCULATE AN INTERPOLATION FACTOR SO THAT ALL
C     ACTIVE CONSTRAINTS ARE REMOVED FROM THE BASIS.
  320 ALPHA = TWO
      HITCON = .FALSE.
      IF (.NOT.(L.LT.NSOLN)) GO TO 360
      DO 350 J=LP1,NSOLN
        ZZ = Z(J)
        IF (.NOT.(ZZ.LE.ZERO)) GO TO 340
        T = X(J)/(X(J)-ZZ)
        IF (.NOT.(T.LT.ALPHA)) GO TO 330
        ALPHA = T
        JCON = J
  330   HITCON = .TRUE.
  340   CONTINUE
  350 CONTINUE
  360 GO TO 1220
  370 CONTINUE
C
C     TO ADD-CONSTRAINTS
C
C     USE COMPUTED ALPHA TO INTERPOLATE BETWEEN LAST
C     FEASIBLE SOLUTION X(*) AND CURRENT UNCONSTRAINED
C     (AND INFEASIBLE) SOLUTION Z(*).
      IF (.NOT.(LP1.LE.NSOLN)) GO TO 390
      DO 380 J=LP1,NSOLN
        X(J) = X(J) + ALPHA*(Z(J)-X(J))
  380 CONTINUE
  390 FEASBL = .FALSE.
      GO TO 410
  400 IF (FEASBL) GO TO 610
C
C     REMOVE COL JCON AND SHIFT COLS JCON+1 THROUGH N TO THE
C     LEFT. SWAP COL JCON INTO THE N-TH POSITION.  THIS ACHIEVES
C     UPPER HESSENBERG FORM FOR THE NONACTIVE CONSTRAINTS AND
C     LEAVES AN UPPER HESSENBERG MATRIX TO RETRIANGULARIZE.
  410 DO 420 I=1,M
        T = W(I,JCON)
        CALL DCOPY(N-JCON, W(I,JCON+1), MDW, W(I,JCON), MDW)
        W(I,N) = T
  420 CONTINUE
C
C     UPDATE PERMUTED INDEX VECTOR TO REFLECT THIS SHIFT AND SWAP.
      ITEMP = IPIVOT(JCON)
      IF (.NOT.(JCON.LT.N)) GO TO 440
      DO 430 I=JCON,NM1
        IPIVOT(I) = IPIVOT(I+1)
  430 CONTINUE
  440 IPIVOT(N) = ITEMP
C
C     SIMILARLY REPERMUTE X(*) VECTOR.
      CALL DCOPY(N-JCON, X(JCON+1), 1, X(JCON), 1)
      X(N) = ZERO
      NSP1 = NSOLN
      NSOLN = NSOLN - 1
      NIV1 = NIV
      NIV = NIV - 1
C
C     RETRIANGULARIZE UPPER HESSENBERG MATRIX AFTER ADDING CONSTRAINTS.
      J = JCON
      I = KRANK + JCON - L
  450 IF (.NOT.(J.LE.NSOLN)) GO TO 570
      IF (.NOT.(ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0)) GO TO 470
      ASSIGN 460 TO IGO938
      GO TO 620
C
C     (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) ZERO-IP1-TO-I-IN-COL-J
  460 GO TO 560
  470 IF (.NOT.(ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1)) GO TO 490
      ASSIGN 480 TO IGO938
      GO TO 620
C
C     (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) ZERO-IP1-TO-I-IN-COL-J
  480 GO TO 560
  490 IF (.NOT.(ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0)) GO TO 510
      CALL DSWAP(NP1, W(I,1), MDW, W(I+1,1), MDW)
      CALL DSWAP(1, SCALE(I), 1, SCALE(I+1), 1)
      ITEMP = ITYPE(I+1)
      ITYPE(I+1) = ITYPE(I)
      ITYPE(I) = ITEMP
C
C     SWAPPED ROW WAS FORMERLY A PIVOT ELT., SO IT WILL
C     BE LARGE ENOUGH TO PERFORM ELIM.
      ASSIGN 500 TO IGO938
      GO TO 620
C
C     ZERO-IP1-TO-I-IN-COL-J
  500 GO TO 560
  510 IF (.NOT.(ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1)) GO TO 550
      T = SCALE(I)*W(I,J)**2/ALSQ
      IF (.NOT.(T.GT.TAU**2*EANORM**2)) GO TO 530
      ASSIGN 520 TO IGO938
      GO TO 620
  520 GO TO 540
  530 CALL DSWAP(NP1, W(I,1), MDW, W(I+1,1), MDW)
      CALL DSWAP(1, SCALE(I), 1, SCALE(I+1), 1)
      ITEMP = ITYPE(I+1)
      ITYPE(I+1) = ITYPE(I)
      ITYPE(I) = ITEMP
      W(I+1,J) = ZERO
  540 CONTINUE
  550 CONTINUE
  560 I = I + 1
      J = J + 1
      GO TO 450
C
C     SEE IF THE REMAINING COEFFS IN THE SOLN SET ARE FEASIBLE.  THEY
C     SHOULD BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.  IF ANY ARE
C     INFEASIBLE IT IS DUE TO ROUNDOFF ERROR.  ANY THAT ARE NON-
C     POSITIVE WILL BE SET TO ZERO AND REMOVED FROM THE SOLN SET.
  570 IF (.NOT.(LP1.LE.NSOLN)) GO TO 590
      DO 580 JCON=LP1,NSOLN
        IF (X(JCON).LE.ZERO) GO TO 600
  580 CONTINUE
  590 FEASBL = .TRUE.
  600 CONTINUE
      GO TO 400
  610 GO TO 1200
  620 CONTINUE
C
C     TO ZERO-IP1-TO-I-IN-COL-J
      IF (.NOT.(W(I+1,J).NE.ZERO)) GO TO 630
      CALL DROTMG(SCALE(I), SCALE(I+1), W(I,J), W(I+1,J), SPARAM)
      W(I+1,J) = ZERO
      CALL DROTM(NP1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, SPARAM)
  630 GO TO 1290
  640 CONTINUE
C
C     TO PERFORM-MULTIPLIER-TEST-AND-DROP-A-CONSTRAINT
      CALL DCOPY(NSOLN, Z, 1, X, 1)
      IF (.NOT.(NSOLN.LT.N)) GO TO 650
      X(NSP1) = ZERO
      CALL DCOPY(N-NSOLN, X(NSP1), 0, X(NSP1), 1)
  650 I = NIV1
  660 IF (.NOT.(I.LE.ME)) GO TO 690
C
C     RECLASSIFY LEAST SQUARES EQATIONS AS EQUALITIES AS
C     NECESSARY.
      IF (.NOT.(ITYPE(I).EQ.0)) GO TO 670
      I = I + 1
      GO TO 680
  670 CALL DSWAP(NP1, W(I,1), MDW, W(ME,1), MDW)
      CALL DSWAP(1, SCALE(I), 1, SCALE(ME), 1)
      ITEMP = ITYPE(I)
      ITYPE(I) = ITYPE(ME)
      ITYPE(ME) = ITEMP
      MEP1 = ME
      ME = ME - 1
  680 GO TO 660
C
C     FORM INNER PRODUCT VECTOR WD(*) OF DUAL COEFFS.
  690 IF (.NOT.(NSP1.LE.N)) GO TO 730
      DO 720 J=NSP1,N
        SM = ZERO
        IF (.NOT.(NSOLN.LT.M)) GO TO 710
        DO 700 I=NSP1,M
          SM = SM + SCALE(I)*W(I,J)*W(I,NP1)
  700   CONTINUE
  710   WD(J) = SM
  720 CONTINUE
  730 GO TO 750
  740 IF (POS .OR. DONE) GO TO 970
C
C     FIND J SUCH THAT WD(J)=WMAX IS MAXIMUM.  THIS DETERMINES
C     THAT THE INCOMING COL J WILL REDUCE THE RESIDUAL VECTOR
C     AND BE POSITIVE.
  750 WMAX = ZERO
      IWMAX = NSP1
      IF (.NOT.(NSP1.LE.N)) GO TO 780
      DO 770 J=NSP1,N
        IF (.NOT.(WD(J).GT.WMAX)) GO TO 760
        WMAX = WD(J)
        IWMAX = J
  760   CONTINUE
  770 CONTINUE
  780 IF (.NOT.(WMAX.LE.ZERO)) GO TO 790
      DONE = .TRUE.
      GO TO 960
C
C     SET DUAL COEFF TO ZERO FOR INCOMING COL.
  790 WD(IWMAX) = ZERO
C
C     WMAX .GT. ZERO, SO OKAY TO MOVE COL IWMAX TO SOLN SET.
C     PERFORM TRANSFORMATION TO RETRIANGULARIZE, AND TEST
C     FOR NEAR LINEAR DEPENDENCE.
C     SWAP COL IWMAX INTO NSOLN-TH POSITION TO MAINTAIN UPPER
C     HESSENBERG FORM OF ADJACENT COLS, AND ADD NEW COL TO
C     TRIANGULAR DECOMPOSITION.
      NSOLN = NSP1
      NSP1 = NSOLN + 1
      NIV = NIV1
      NIV1 = NIV + 1
      IF (.NOT.(NSOLN.NE.IWMAX)) GO TO 800
      CALL DSWAP(M, W(1,NSOLN), 1, W(1,IWMAX), 1)
      WD(IWMAX) = WD(NSOLN)
      WD(NSOLN) = ZERO
      ITEMP = IPIVOT(NSOLN)
      IPIVOT(NSOLN) = IPIVOT(IWMAX)
      IPIVOT(IWMAX) = ITEMP
C
C     REDUCE COL NSOLN SO THAT THE MATRIX OF NONACTIVE
C     CONSTRAINTS VARIABLES IS TRIANGULAR.
  800 J = M
  810 IF (.NOT.(J.GT.NIV)) GO TO 870
      JM1 = J - 1
      JP = JM1
C
C     WHEN OPERATING NEAR THE ME LINE, TEST TO SEE IF THE PIVOT ELT.
C     IS NEAR ZERO.  IF SO, USE THE LARGEST ELT. ABOVE IT AS THE PIVOT.
C     THIS IS TO MAINTAIN THE SHARP INTERFACE BETWEEN WEIGHTED AND
C     NON-WEIGHTED ROWS IN ALL CASES.
      IF (.NOT.(J.EQ.MEP1)) GO TO 850
      IMAX = ME
      AMAX = SCALE(ME)*W(ME,NSOLN)**2
  820 IF (.NOT.(JP.GE.NIV)) GO TO 840
      T = SCALE(JP)*W(JP,NSOLN)**2
      IF (.NOT.(T.GT.AMAX)) GO TO 830
      IMAX = JP
      AMAX = T
  830 JP = JP - 1
      GO TO 820
  840 JP = IMAX
  850 IF (.NOT.(W(J,NSOLN).NE.ZERO)) GO TO 860
      CALL DROTMG(SCALE(JP), SCALE(J), W(JP,NSOLN), W(J,NSOLN), SPARAM)
      W(J,NSOLN) = ZERO
      CALL DROTM(NP1-NSOLN, W(JP,NSP1), MDW, W(J,NSP1), MDW, SPARAM)
  860 J = JM1
      GO TO 810
C
C     SOLVE FOR Z(NSOLN)=PROPOSED NEW VALUE FOR X(NSOLN).
C     TEST IF THIS IS NONPOSITIVE OR TOO LARGE.
C     IF THIS WAS TRUE OR IF THE PIVOT TERM WAS ZERO REJECT
C     THE COL AS DEPENDENT.
  870 IF (.NOT.(W(NIV,NSOLN).NE.ZERO)) GO TO 890
      ISOL = NIV
      ASSIGN 880 TO IGO897
      GO TO 980
C
C     TEST-PROPOSED-NEW-COMPONENT
  880 GO TO 940
  890 IF (.NOT.(NIV.LE.ME .AND. W(MEP1,NSOLN).NE.ZERO)) GO TO 920
C
C     TRY TO ADD ROW MEP1 AS AN ADDITIONAL EQUALITY CONSTRAINT.
C     CHECK SIZE OF PROPOSED NEW SOLN COMPONENT.
C     REJECT IT IF IT IS TOO LARGE.
      ISOL = MEP1
      ASSIGN 900 TO IGO897
      GO TO 980
C
C     TEST-PROPOSED-NEW-COMPONENT
  900 IF (.NOT.(POS)) GO TO 910
C
C     SWAP ROWS MEP1 AND NIV, AND SCALE FACTORS FOR THESE ROWS.
      CALL DSWAP(NP1, W(MEP1,1), MDW, W(NIV,1), MDW)
      CALL DSWAP(1, SCALE(MEP1), 1, SCALE(NIV), 1)
      ITEMP = ITYPE(MEP1)
      ITYPE(MEP1) = ITYPE(NIV)
      ITYPE(NIV) = ITEMP
      ME = MEP1
      MEP1 = ME + 1
  910 GO TO 930
  920 POS = .FALSE.
  930 CONTINUE
  940 IF (POS) GO TO 950
      NSP1 = NSOLN
      NSOLN = NSOLN - 1
      NIV1 = NIV
      NIV = NIV - 1
  950 CONTINUE
  960 GO TO 740
  970 GO TO 1250
  980 CONTINUE
C
C     TO TEST-PROPOSED-NEW-COMPONENT
      Z2 = W(ISOL,NP1)/W(ISOL,NSOLN)
      Z(NSOLN) = Z2
      POS = Z2.GT.ZERO
      IF (.NOT.(Z2*EANORM.GE.BNORM .AND. POS)) GO TO 990
      POS = .NOT.(BLOWUP*Z2*EANORM.GE.BNORM)
  990 GO TO 1280
 1000 CONTINUE
C     TO COMPUTE-FINAL-SOLUTION
C
C     SOLVE SYSTEM, STORE RESULTS IN X(*).
C
      ASSIGN 1010 TO IGO958
      GO TO 1110
C     SOLVE-SYSTEM
 1010 CALL DCOPY(NSOLN, Z, 1, X, 1)
C
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO X(*) IF KRANK.LT.L
      IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.L)) GO TO 1030
      DO 1020 I=1,KRANK
        CALL DH12(2, I, KRP1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
 1020 CONTINUE
C
C     FILL IN TRAILING ZEROES FOR CONSTRAINED VARIABLES NOT IN SOLN.
 1030 IF (.NOT.(NSOLN.LT.N)) GO TO 1040
      X(NSP1) = ZERO
      CALL DCOPY(N-NSOLN, X(NSP1), 0, X(NSP1), 1)
C
C     REPERMUTE SOLN VECTOR TO NATURAL ORDER.
 1040 DO 1070 I=1,N
        J = I
 1050   IF (IPIVOT(J).EQ.I) GO TO 1060
        J = J + 1
        GO TO 1050
 1060   IPIVOT(J) = IPIVOT(I)
        IPIVOT(I) = J
        CALL DSWAP(1, X(J), 1, X(I), 1)
 1070 CONTINUE
C
C     RESCALE THE SOLN USING THE COL SCALING.
      DO 1080 J=1,N
        X(J) = X(J)*D(J)
 1080 CONTINUE
      IF (.NOT.(NSOLN.LT.M)) GO TO 1100
      DO 1090 I=NSP1,M
        T = W(I,NP1)
        IF (I.LE.ME) T = T/ALAMDA
        T = (SCALE(I)*T)*T
        RNORM = RNORM + T
 1090 CONTINUE
 1100 RNORM = DSQRT(RNORM)
      GO TO 1210
C
C     TO SOLVE-SYSTEM
C
 1110 CONTINUE
      IF (.NOT.(DONE)) GO TO 1120
      ISOL = 1
      GO TO 1130
 1120 ISOL = LP1
 1130 IF (.NOT.(NSOLN.GE.ISOL)) GO TO 1190
C
C     COPY RT. HAND SIDE INTO TEMP VECTOR TO USE OVERWRITING METHOD.
      CALL DCOPY(NIV, W(1,NP1), 1, TEMP, 1)
      DO 1180 JJ=ISOL,NSOLN
        J = NSOLN - JJ + ISOL
        IF (.NOT.(J.GT.KRANK)) GO TO 1140
        I = NIV - JJ + ISOL
        GO TO 1150
 1140   I = J
 1150   IF (.NOT.(J.GT.KRANK .AND. J.LE.L)) GO TO 1160
        Z(J) = ZERO
        GO TO 1170
 1160   Z(J) = TEMP(I)/W(I,J)
        CALL DAXPY(I-1, -Z(J), W(1,J), 1, TEMP, 1)
 1170   CONTINUE
 1180 CONTINUE
 1190 GO TO 1270
 1200 GO TO IGO986, (40)
 1210 GO TO IGO980, (90)
 1220 GO TO IGO991, (30)
 1230 GO TO IGO998, (10)
 1240 GO TO IGO995, (20)
 1250 GO TO IGO983, (60)
 1260 GO TO IGO977, (220)
 1270 GO TO IGO958, (310, 1010)
 1280 GO TO IGO897, (880, 900)
 1290 GO TO IGO938, (460, 480, 500, 520)
      END
*DECK DWNNLS
      SUBROUTINE DWNNLS(W,MDW,ME,MA,N,L,PRGOPT,X,RNORM,MODE,IWORK,WORK)
C***BEGIN PROLOGUE  DWNNLS
C***DATE WRITTEN   790701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  K1A2A
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(WNNLS-S DWNNLS-D),
C             CONSTRAINED LEAST SQUARES,CURVE FITTING,DATA FITTING,
C             EQUALITY CONSTRAINTS,INEQUALITY CONSTRAINTS,
C             NONNEGATIVITY CONSTRAINTS,QUADRATIC PROGRAMMING
C***AUTHOR  HANSON, R. J., (SNLA)
C           HASKELL, K. H., (SNLA)
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality constraints and nonnegativity constraints on
C            selected variables.
C***DESCRIPTION
C
C       *** Double Precision version of WNNLS *****
C
C     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),IWORK(M+N),WORK(M+5*N)
C
C     Written by Karen H. Haskell, Sandia Laboratories,
C     and R.J. Hanson, Sandia Laboratories.
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem.  Suppose there are given matrices E and A of
C     respective dimensions ME by N and MA by N, and vectors F
C     and B of respective lengths ME and MA.  This subroutine
C     solves the problem
C
C               EX = F, (equations to be exactly satisfied)
C
C               AX = B, (equations to be approximately satisfied,
C                        in the least squares sense)
C
C               subject to components L+1,...,N nonnegative
C
C     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
C
C     The problem is reposed as problem DWNNLS
C
C               (WT*E)X = (WT*F)
C               (   A)    (   B), (least squares)
C               subject to components L+1,...,N nonnegative.
C
C     The subprogram chooses the heavy weight (or penalty parameter) WT.
C
C     The parameters for DWNNLS are
C
C     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     W(*,*),MDW,  The array W(*,*) is double subscripted with first
C     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
C                  discussion let us call M = ME + MA.  Then MDW
C                  must satisfy MDW.GE.M.  The condition MDW.LT.M
C                  is an error.
C
C                  The array W(*,*) contains the matrices and vectors
C
C                       (E  F)
C                       (A  B)
C
C                  in rows and columns 1,...,M and 1,...,N+1
C                  respectively.  Columns 1,...,L correspond to
C                  unconstrained variables X(1),...,X(L).  The
C                  remaining variables are constrained to be
C                  nonnegative. The condition L.LT.0 or L.GT.N is
C                  an error.
C
C     PRGOPT(*)    This double precision array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case LINK=1 and the values KEY and DATA SET
C                  are not referenced. The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=LINK1 (link to first entry of next group)
C               .  PRGOPT(2)=KEY1 (key to the option change)
C               .  PRGOPT(3)=DATA VALUE (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2)=DATA VALUE
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK.GT.NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000 an error
C                  message is printed and the subprogram returns.
C
C                  OPTIONS..
C
C                  KEY=6
C                         Scale the nonzero columns of the
C                  entire data matrix
C                  (E)
C                  (A)
C                  to have length one. The DATA SET for
C                  this option is a single value.  It must
C                  be nonzero if unit length column scaling is
C                  desired.
C
C                  KEY=7
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  with a user-provided diagonal matrix.
C                  The DATA SET for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=8
C                         Change the rank determination tolerance from
C                  the nominal value of SQRT(SRELPR).  This quantity
C                  can be no smaller than SRELPR, The arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The DATA SET for this option
C                  is the new tolerance.
C
C                  KEY=9
C                         Change the blow-up parameter from the
C                  nominal value of SQRT(SRELPR).  The reciprocal of
C                  this parameter is used in rejecting solution
C                  components as too large when a variable is
C                  first brought into the active set.  Too large
C                  means that the proposed component times the
C                  reciprocal of the parameter is not less than
C                  the ratio of the norms of the right-side
C                  vector and the data matrix.
C                  This parameter can be no smaller than SRELPR,
C                  the arithmetic-storage precision.
C
C                  For example, suppose we want to provide
C                  a diagonal matrix to scale the problem
C                  matrix and change the tolerance used for
C                  determining linear dependence of dropped col
C                  vectors.  For these options the dimensions of
C                  PRGOPT(*) must be at least N+6.  The FORTRAN
C                  statements defining these options would
C                  be as follows.
C
C                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
C                  PRGOPT(2)=7 (user-provided scaling key)
C
C                  CALL DCOPY(N,D,1,PRGOPT(3),1) (copy the N
C                  scaling factors from a user array called D(*)
C                  into PRGOPT(3)-PRGOPT(N+2))
C
C                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
C                  PRGOPT(N+4)=8 (linear dependence tolerance key)
C                  PRGOPT(N+5)=... (new value of the tolerance)
C
C                  PRGOPT(N+6)=1 (no more options to change)
C
C
C     IWORK(1),    The amounts of working storage actually allocated
C     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
C                  respectively.  These quantities are compared with
C                  the actual amounts of storage needed for DWNNLS( ).
C                  Insufficient storage allocated for either WORK(*)
C                  or IWORK(*) is considered an error.  This feature
C                  was included in DWNNLS( ) because miscalculating
C                  the storage formulas for WORK(*) and IWORK(*)
C                  might very well lead to subtle and hard-to-find
C                  execution errors.
C
C                  The length of WORK(*) must be at least
C
C                  LW = ME+MA+5*N
C                  This test will not be made if IWORK(1).LE.0.
C
C                  The length of IWORK(*) must be at least
C
C                  LIW = ME+MA+N
C                  This test will not be made if IWORK(2).LE.0.
C
C     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     X(*)         An array dimensioned at least N, which will
C                  contain the N components of the solution vector
C                  on output.
C
C     RNORM        The residual norm of the solution.  The value of
C                  RNORM contains the residual vector length of the
C                  equality constraints and least squares equations.
C
C     MODE         The value of MODE indicates the success or failure
C                  of the subprogram.
C
C                  MODE = 0  Subprogram completed successfully.
C
C                       = 1  Max. number of iterations (equal to
C                            3*(N-L)) exceeded. Nearly all problems
C                            should complete in fewer than this
C                            number of iterations. An approximate
C                            solution and its corresponding residual
C                            vector length are in X(*) and RNORM.
C
C                       = 2  Usage error occurred.  The offending
C                            condition is noted with the error
C                            processing subprogram, XERROR( ).
C
C     User-designated
C     Working arrays..
C
C     WORK(*)      A double precision working array of length at least
C                  M + 5*N.
C
C     IWORK(*)     An integer-valued working array of length at least
C                  M+N.
C***REFERENCES  K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, SAND77-0552, JUNE 1978.
C               K.H. HASKELL AND R.J. HANSON, *SELECTED ALGORITHMS FOR
C                 THE LINEARLY CONSTRAINED LEAST SQUARES PROBLEM--
C                 A USERS GUIDE*, SAND78-1290, AUGUST 1979.
C               K.H. HASKELL AND R.H. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, MATH. PROG. 21 (1981),
C                 PP. 98-118.
C               R.J. HANSON AND K.H. HASKELL, *TWO ALGORITHMS FOR THE
C                 LINEARLY CONSTRAINED LEAST SQUARES PROBLEM*, ACM
C                 TRANS. ON MATH. SOFTWARE, SEPT. 1982.
C***ROUTINES CALLED  DWNLSM,XERROR,XERRWV
C***END PROLOGUE  DWNNLS
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /DOUBLE PRECISION (12 BLANKS)/DOUBLE PRECISION/,/,
C     DUMMY/,SNGL(DUMMY)/
C
      INTEGER IOPT, IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW,
     *     ME, MODE, N, NERR
      DOUBLE PRECISION DUMMY, PRGOPT(*), RNORM, W(MDW,*), WORK(*),
     *     X(*)
C
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 60
C        BEGIN BLOCK PERMITTING ...EXITS TO 20
C           BEGIN BLOCK PERMITTING ...EXITS TO 10
C***FIRST EXECUTABLE STATEMENT  DWNNLS
               MODE = 0
C     .........EXIT
               IF (MA + ME .LE. 0 .OR. N .LE. 0) GO TO 60
C           ...EXIT
               IF (IWORK(1) .LE. 0) GO TO 10
               LW = ME + MA + 5*N
C           ...EXIT
               IF (IWORK(1) .GE. LW) GO TO 10
               NERR = 2
               IOPT = 1
               CALL XERRWV(
     *         '  DWNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR WORK(*),
     * NEED LW=I1 BELOW',70,NERR,IOPT,1,LW,0,0,DUMMY,DUMMY)
               MODE = 2
C     .........EXIT
               GO TO 60
   10       CONTINUE
C        ...EXIT
            IF (IWORK(2) .LE. 0) GO TO 20
            LIW = ME + MA + N
C        ...EXIT
            IF (IWORK(2) .GE. LIW) GO TO 20
            NERR = 2
            IOPT = 1
            CALL XERRWV(
     *      '  DWNNLS( ), INSUFFICIENT STORAGE ALLOCATED FOR IWORK(*), N
     *EED LIW=I1 BELOW',72,NERR,IOPT,1,LIW,0,0,DUMMY,DUMMY)
            MODE = 2
C     ......EXIT
            GO TO 60
   20    CONTINUE
         IF (MDW .GE. ME + MA) GO TO 30
            NERR = 1
            IOPT = 1
            CALL XERROR('  DWNNLS( ), THE VALUE MDW.LT.ME+MA IS AN ERROR
     *      ',44,NERR,IOPT)
            MODE = 2
         GO TO 50
   30    CONTINUE
         IF (0 .LE. L .AND. L .LE. N) GO TO 40
            NERR = 2
            IOPT = 1
            CALL XERROR('  DWNNLS( ), L.LE.0.AND.L.LE.N IS REQUIRED',39,
     *                  NERR,IOPT)
            MODE = 2
         GO TO 50
   40    CONTINUE
C
C           THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
C           WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
C           REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
C
            L1 = N + 1
            L2 = L1 + N
            L3 = L2 + ME + MA
            L4 = L3 + N
            L5 = L4 + N
C
            CALL DWNLSM(W,MDW,ME,MA,N,L,PRGOPT,X,RNORM,MODE,IWORK,
     *                  IWORK(L1),WORK(1),WORK(L1),WORK(L2),WORK(L3),
     *                  WORK(L4),WORK(L5))
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(FDUMP-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A2
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(ISAMAX-S IDAMAX-D ICAMAX-C),
C             LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find the smallest index of that component of a d.p. vector
C            having the maximum magnitude.
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  IDAMAX
C
      DOUBLE PRECISION DX(*),DMAX,XMAG
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
*DECK J4SAVE
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C***ROUTINES CALLED  (NONE)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  1 August 1985
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERABT
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERABT-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Abort program execution and print error message.
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  1 August 1982
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
*DECK XERCTL
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERCTL-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allow user control over handling of errors.
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
*DECK XERPRT
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERPRT-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Print error messages.
C***DESCRIPTION
C
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  1 August 1985
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      SAVE F, G, LA, LCOM, LBLANK
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
*DECK XERROR
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERROR-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Process an error (diagnostic) message.
C***DESCRIPTION
C
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C    1                43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
*DECK XERRWV
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERRWV-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Process an error message allowing 2 integer and 2 real
C            values to be included in the message.
C***DESCRIPTION
C
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  1 August 1985
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
*DECK XERSAV
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XERSAV-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Record that an error has occurred.
C***DESCRIPTION
C
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  1 August 1985
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
*DECK XGETUA
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(XGETUA-A),ERROR
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., 'XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE', SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
*DECK DLSEI
      SUBROUTINE DLSEI(W,MDW,ME,MA,MG,N,PRGOPT,X,RNORME,RNORML,MODE,WS,
     *   IP)
C***BEGIN PROLOGUE  DLSEI
C***DATE WRITTEN   790701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  K1A2A,D9
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(LSEI-S DLSEI-D),
C             CONSTRAINED LEAST SQUARES,CURVE FITTING,DATA FITTING,
C             EQUALITY CONSTRAINTS,INEQUALITY CONSTRAINTS,
C             QUADRATIC PROGRAMMING
C***AUTHOR  HANSON, R. J., (SNLA)
C           HASKELL, K. H., (SNLA)
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality and inequality constraints, and optionally compute
C            a covariance matrix.
C***DESCRIPTION
C
C       **** Double Precision version of LSEI ****
C
C     DIMENSION W(MDW,N+1),PRGOPT(*),X(N),
C     WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     above, K=MAX(MA+MG,N).
C
C     Written by R. J. Hanson and K. H. Haskell.  For further math.
C     and algorithmic details see Sandia Laboratories Tech. Repts.
C     SAND-77-0552, (1978), and SAND-78-1290, (1979), or Math. Prog.,
C     Vol. 21 (1981) pp. 98-118 and ACM Trans. on Math. Software,
C     Sept. 1982.
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem with both equality and inequality constraints, and, if the
C     user requests, obtains a covariance matrix of the solution
C     parameters.
C
C     Suppose there are given matrices E, A and G of respective
C     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
C     respective lengths ME, MA and MG.  This subroutine solves the
C     linearly constrained least squares problem
C
C                   EX = F, (E ME by N) (equations to be exactly
C                                       satisfied)
C                   AX = B, (A MA by N) (equations to be
C                                       approximately satisfied,
C                                       least squares sense)
C                   GX .GE. H,(G MG by N) (inequality constraints)
C
C     The inequalities GX .GE. H mean that every component of the
C     product GX must be .GE. the corresponding component of H.
C
C     In case the equality constraints cannot be satisfied, a
C     generalized inverse solution residual vector length is obtained
C     for F-EX.  This is the minimal length possible for F-EX.
C
C
C     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
C     rank of the matrix E is estimated during the computation.  We call
C     this value KRANKE.  It is an output parameter in IP(1) defined
C     below.  Using a generalized inverse solution of EX=F, a reduced
C     least squares problem with inequality constraints is obtained.
C     The tolerances used in these tests for determining the rank
C     of E and the rank of the reduced least squares problem are
C     given in Sandia Tech. Rept. SAND-78-1290.  They can be
C     modified by the user if new values are provided in
C     the option list of the array PRGOPT(*).
C
C     The user must dimension all arrays appearing in the call list..
C     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     where K=MAX(MA+MG,N).  This allows for a solution of a range of
C     problems in the given working space.  The dimension of WS(*)
C     given is a necessary overestimate.  Once a particular problem
C     has been run, the output parameter IP(3) gives the actual
C     dimension required for that problem.
C
C     The parameters for DLSEI( ) are
C
C     Input.. All TYPE REAL variables are DOUBLE PRECISION
C
C     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
C     ME,MA,MG,N    first dimensioning parameter equal to MDW.
C                   For this discussion let us call M = ME+MA+MG.  Then
C                   MDW must satisfy MDW .GE. M.  The condition
C                   MDW .LT. M is an error.
C
C                   The array W(*,*) contains the matrices and vectors
C
C                                  (E  F)
C                                  (A  B)
C                                  (G  H)
C
C                   in rows and columns 1,...,M and 1,...,N+1
C                   respectively.
C
C                   The integers ME, MA, and MG are the
C                   respective matrix row dimensions
C                   of E, A and G.  Each matrix has N columns.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case, LINK=1 and the values KEY and DATA SET
C                  are not referenced.  The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=link1 (link to first entry of next group)
C               .  PRGOPT(2)=key1 (key to the option change)
C               .  PRGOPT(3)=data value (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=link2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=key2 (key to the option change)
C               .  PRGOPT(LINK1+2)=data value
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK .GT. NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array, a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000, an error
C                  message is printed and the subprogram returns.
C
C                  Options..
C
C                  KEY=1
C                         Compute in W(*,*) the N by N
C                  covariance matrix of the solution variables
C                  as an output parameter.  Nominally the
C                  covariance matrix will not be computed.
C                  (This requires no user input.)
C                  The data set for this option is a single value.
C                  It must be nonzero when the covariance matrix
C                  is desired.  If it is zero, the covariance
C                  matrix is not computed.  When the covariance matrix
C                  is computed, the first dimensioning parameter
C                  of the array W(*,*) must satisfy MDW .GE. MAX0(M,N).
C
C                  KEY=10
C                         Suppress scaling of the inverse of the
C                  normal matrix by the scale factor RNORM**2/
C                  MAX(1, no. of degrees of freedom).  This option
C                  only applies when the option for computing the
C                  covariance matrix (KEY=1) is used.  With KEY=1 and
C                  KEY=10 used as options the unscaled inverse of the
C                  normal matrix is returned in W(*,*).
C                  The data set for this option is a single value.
C                  When it is nonzero no scaling is done.  When it is
C                  zero scaling is done.  The nominal case is to do
C                  scaling so if option (KEY=1) is used alone, the
C                  matrix will be scaled on output.
C
C                  KEY=2
C                         Scale the nonzero columns of the
C                         entire data matrix.
C                  (E)
C                  (A)
C                  (G)
C
C                  to have length one.  The data set for this
C                  option is a single value.  It must be
C                  nonzero if unit length column scaling
C                  is desired.
C
C                  KEY=3
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  (G)
C
C                  with a user-provided diagonal matrix.
C                  The data set for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=4
C                         Change the rank determination tolerance for
C                  the equality constraint equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity SRELPR is the
C                  largest positive number such that T=1.+SRELPR
C                  satisfies T .EQ. 1.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  KEY=5
C                         Change the rank determination tolerance for
C                  the reduced least squares equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  For example, suppose we want to change
C                  the tolerance for the reduced least squares
C                  problem, compute the covariance matrix of
C                  the solution parameters, and provide
C                  column scaling for the data matrix.  For
C                  these options the dimension of PRGOPT(*)
C                  must be at least N+9.  The Fortran statements
C                  defining these options would be as follows:
C
C                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
C                  PRGOPT(2)=1 (covariance matrix key)
C                  PRGOPT(3)=1 (covariance matrix wanted)
C
C                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
C                  PRGOPT(5)=5 (least squares equas.  tolerance key)
C                  PRGOPT(6)=... (new value of the tolerance)
C
C                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
C                  PRGOPT(8)=3 (user-provided column scaling key)
C
C                  CALL DCOPY(N,D,1,PRGOPT(9),1)  (Copy the N
C                    scaling factors from the user array D(*)
C                    to PRGOPT(9)-PRGOPT(N+8))
C
C                  PRGOPT(N+9)=1 (no more options to change)
C
C                  The contents of PRGOPT(*) are not modified
C                  by the subprogram.
C                  The options for DWNNLS( ) can also be included
C                  in this array.  The values of KEY recognized
C                  by DWNNLS( ) are 6, 7 and 8.  Their functions
C                  are documented in the usage instructions for
C                  subroutine DWNNLS( ).  Normally these options
C                  do not need to be modified when using DLSEI( ).
C
C     IP(1),       The amounts of working storage actually
C     IP(2)        allocated for the working arrays WS(*) and
C                  IP(*), respectively.  These quantities are
C                  compared with the actual amounts of storage
C                  needed by DLSEI( ).  Insufficient storage
C                  allocated for either WS(*) or IP(*) is an
C                  error.  This feature was included in DLSEI( )
C                  because miscalculating the storage formulas
C                  for WS(*) and IP(*) might very well lead to
C                  subtle and hard-to-find execution errors.
C
C                  The length of WS(*) must be at least
C
C                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
C
C                  where K = max(MA+MG,N)
C                  This test will not be made if IP(1).LE.0.
C
C                  The length of IP(*) must be at least
C
C                  LIP = MG+2*N+2
C                  This test will not be made if IP(2).LE.0.
C
C     Output..All TYPE REAL variables are DOUBLE PRECISION
C
C     X(*),RNORME,  The array X(*) contains the solution parameters
C     RNORML        if the integer output flag MODE = 0 or 1.
C                   The definition of MODE is given directly below.
C                   When MODE = 0 or 1, RNORME and RNORML
C                   respectively contain the residual vector
C                   Euclidean lengths of F - EX and B - AX.  When
C                   MODE=1 the equality constraint equations EX=F
C                   are contradictory, so RNORME .NE. 0.  The residual
C                   vector F-EX has minimal Euclidean length.  For
C                   MODE .GE. 2, none of these parameters are
C                   defined.
C
C     MODE          Integer flag that indicates the subprogram
C                   status after completion.  If MODE .GE. 2, no
C                   solution has been computed.
C
C                   MODE =
C
C                   0  Both equality and inequality constraints
C                      are compatible and have been satisfied.
C
C                   1  Equality constraints are contradictory.
C                      A generalized inverse solution of EX=F was used
C                      to minimize the residual vector length F-EX.
C                      In this sense, the solution is still meaningful.
C
C                   2  Inequality constraints are contradictory.
C
C                   3  Both equality and inequality constraints
C                      are contradictory.
C
C                   The following interpretation of
C                   MODE=1,2 or 3 must be made.  The
C                   sets consisting of all solutions
C                   of the equality constraints EX=F
C                   and all vectors satisfying GX .GE. H
C                   have no points in common.  (In
C                   particular this does not say that
C                   each individual set has no points
C                   at all, although this could be the
C                   case.)
C
C                   4  Usage error occurred.  The value
C                      of MDW is .LT. ME+MA+MG, MDW is
C                      .LT. N and a covariance matrix is
C                      requested, or the option vector
C                      PRGOPT(*) is not properly defined,
C                      or the lengths of the working arrays
C                      WS(*) and IP(*), when specified in
C                      IP(1) and IP(2) respectively, are not
C                      long enough.
C
C     W(*,*)        The array W(*,*) contains the N by N symmetric
C                   covariance matrix of the solution parameters,
C                   provided this was requested on input with
C                   the option vector PRGOPT(*) and the output
C                   flag is returned with MODE = 0 or 1.
C
C     IP(*)         The integer working array has three entries
C                   that provide rank and working array length
C                   information after completion.
C
C                      IP(1) = rank of equality constraint
C                              matrix.  Define this quantity
C                              as KRANKE.
C
C                      IP(2) = rank of reduced least squares
C                              problem.
C
C                      IP(3) = the amount of storage in the
C                              working array WS(*) that was
C                              actually used by the subprogram.
C                              The formula given above for the length
C                              of WS(*) is a necessary overestimate.
C                              If exactly the same problem matrices
C                              are used in subsequent executions,
C                              the declared dimension of WS(*) can
C                              be reduced to this output value.
C     User Designated
C     Working Arrays..
C
C     WS(*),IP(*)              These are respectively type real
C                              and type integer working arrays.
C                              Their required minimal lengths are
C                              given above.
C***REFERENCES  K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, SAND77-0552, JUNE 1978.
C               K.H. HASKELL AND R.J. HANSON, *SELECTED ALGORITHMS FOR
C                 THE LINEARLY CONSTRAINED LEAST SQUARES PROBLEM--
C                 A USERS GUIDE*, SAND78-1290, AUGUST 1979.
C               K.H. HASKELL AND R.J. HANSON, *AN ALGORITHM FOR
C                 LINEAR LEAST SQUARES PROBLEMS WITH EQUALITY AND
C                 NONNEGATIVITY CONSTRAINTS*, MATH. PROG. 21 (1981),
C                 PP. 98-118.
C               R.J. HANSON AND K.H. HASKELL, *TWO ALGORITHMS FOR THE
C                 LINEARLY CONSTRAINED LEAST SQUARES PROBLEM*, ACM
C                 TRANS. ON MATH. SOFTWARE, SEPT. 1982.
C***ROUTINES CALLED  DASUM,DAXPY,DCOPY,DDOT,DH12,DLSI,DNRM2,DSCAL,DSWAP,
C                    XERROR,XERRWV
C***END PROLOGUE  DLSEI
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (START EDITING AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/SASUM/DASUM/,/SDOT/DDOT/,
C     /SNRM2/DNRM2/,/ SQRT/ DSQRT/,/ ABS/ DABS/,/AMAX1/DMAX1/,
C     /SCOPY/DCOPY/,/SSCAL/DSCAL/,/SAXPY/DAXPY/,/SSWAP/DSWAP/,/E0/D0/,
C     /, DUMMY/,SNGL(DUMMY)/,/SRELPR/DRELPR/
C
C
C     SUBROUTINES CALLED
C
C     DLSI           PART OF THIS PACKAGE.  SOLVES A
C                   CONSTRAINED LEAST SQUARES PROBLEM WITH
C                   INEQUALITY CONSTRAINTS.
C
C++
C     DDOT,DSCAL,   SUBROUTINES FROM THE BLAS PACKAGE.
C     DAXPY,DASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     DCOPY,DNRM2,
C     DSWAP
C
C     DH12           SUBROUTINE TO CONSTRUCT AND APPLY A
C                    HOUSEHOLDER TRANSFORMATION.
C
C     XERROR        FROM SLATEC ERROR PROCESSING PACKAGE.
C                   THIS IS DOCUMENTED IN SANDIA TECH. REPT.,
C                   SAND78-1189.
C
C     SUBROUTINE DLSEI(W,MDW,ME,MA,MG,N,PRGOPT,X,RNORME,RNORML,MODE,WS,
C    1 IP)
C
C     REVISED MAY, 1982.
C
      DOUBLE PRECISION W(MDW,*), PRGOPT(*), X(*), WS(*), RNORME, RNORML
      INTEGER IP(3)
      DOUBLE PRECISION DUMMY, ENORM, DRELPR, FNORM, GAM, HALF, ONE, RB,
     1RN, RNMAX, SIZE, SN, SNMAX, T, TAU, UJ, UP, VJ, XNORM, XNRME, ZERO
      DOUBLE PRECISION DASUM, DDOT, DNRM2, DSQRT, DABS, DMAX1
      LOGICAL COV
      DATA ZERO /0.D0/, DRELPR /0.D0/, HALF /0.5D0/, ONE /1.D0/
C
C     CHECK THAT ENOUGH STORAGE WAS ALLOCATED IN WS(*) AND IP(*).
C***FIRST EXECUTABLE STATEMENT  DLSEI
      IF (.NOT.(IP(1).GT.0)) GO TO 20
      LCHK = 2*(ME+N) + MAX0(MA+MG,N) + (MG+2)*(N+7)
      IF (.NOT.(IP(1).LT.LCHK)) GO TO 10
      MODE = 4
      NERR = 2
      IOPT = 1
      CALL XERRWV('DLSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR WS(*), N
     1EED LW=I1 BELOW', 67, NERR, IOPT, 1, LCHK, 0, 0,
     2 SNGL(DUMMY),SNGL(DUMMY))
      RETURN
   10 CONTINUE
   20 IF (.NOT.(IP(2).GT.0)) GO TO 40
      LCHK = MG + 2*N + 2
      IF (.NOT.(IP(2).LT.LCHK)) GO TO 30
      MODE = 4
      NERR = 2
      IOPT = 1
      CALL XERRWV('DLSEI( ), INSUFFICIENT STORAGE ALLOCATED FOR IP(*), N
     1EED LIP=I1 BELOW', 68, NERR, IOPT, 1, LCHK, 0, 0,
     2SNGL(DUMMY),SNGL(DUMMY))
      RETURN
   30 CONTINUE
C
C     COMPUTE MACHINE PRECISION=DRELPR ONLY WHEN NECESSARY.
   40 IF (.NOT.(DRELPR.EQ.ZERO)) GO TO 70
      DRELPR = ONE
   50 IF (ONE+DRELPR.EQ.ONE) GO TO 60
      DRELPR = DRELPR*HALF
      GO TO 50
   60 DRELPR = DRELPR + DRELPR
C
C     COMPUTE NUMBER OF POSSIBLE RIGHT MULTIPLYING HOUSEHOLDER
C     TRANSFORMATIONS.
   70 M = ME + MA + MG
      MODE = 0
      IF (N.LE.0 .OR. M.LE.0) RETURN
      IF (.NOT.(MDW.LT.M)) GO TO 80
      NERR = 1
      IOPT = 1
      CALL XERROR('DLSEI( ), MDW.LT.ME+MA+MG IS AN ERROR', 36, NERR,
     1 IOPT)
      MODE = 4
      RETURN
   80 NP1 = N + 1
      KRANKE = MIN0(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
C
C     PROCESS-OPTION-VECTOR
      ASSIGN 90 TO IGO990
      GO TO 480
   90 IF (.NOT.(COV .AND. MDW.LT.N)) GO TO 100
      NERR = 2
      IOPT = 1
      CALL XERROR('DLSEI( ), MDW.LT.N, WHEN COV MATRIX NEEDED, IS AN ERR
     1OR', 54,    NERR, IOPT)
      MODE = 4
      RETURN
  100 L = KRANKE
C
C     COMPUTE NORM OF EQUALITY CONSTRAINT MATRIX AND RT SIDE.
      ENORM = ZERO
      DO 110 J=1,N
        ENORM = DMAX1(ENORM,DASUM(ME,W(1,J),1))
  110 CONTINUE
      FNORM = DASUM(ME,W(1,NP1),1)
      IF (.NOT.(L.GT.0)) GO TO 190
      SNMAX = ZERO
      RNMAX = ZERO
      DO 180 I=1,L
C
C     COMPUTE MAXIMUM RATIO OF VECTOR LENGTHS. PARTITION
C     IS AT COL. I.
        DO 150 K=I,ME
          SN = DDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
          RN = DDOT(I-1,W(K,1),MDW,W(K,1),MDW)
          IF (.NOT.(RN.EQ.ZERO .AND. SN.GT.SNMAX)) GO TO 120
          SNMAX = SN
          IMAX = K
          GO TO 140
  120     IF (.NOT.(K.EQ.I .OR. (SN*RNMAX.GT.RN*SNMAX))) GO TO 130
          SNMAX = SN
          RNMAX = RN
          IMAX = K
  130     CONTINUE
  140     CONTINUE
  150   CONTINUE
C
C     INTERCHANGE ROWS IF NECESSARY.
        IF (I.NE.IMAX) CALL DSWAP(NP1, W(I,1), MDW, W(IMAX,1), MDW)
        IF (.NOT.(SNMAX.GT.TAU**2*RNMAX)) GO TO 160
C
C     ELIMINATE ELEMS I+1,...,N IN ROW I.
        CALL DH12(1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW, 1,
     1   M-I)
        GO TO 170
  160   KRANKE = I - 1
        GO TO 200
  170   CONTINUE
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
C     SAVE DIAG. TERMS OF LOWER TRAP. MATRIX.
      CALL DCOPY(KRANKE, W, MDW+1, WS(KRANKE+1), 1)
C
C     USE HOUSEHOLDER TRANS FROM LEFT TO ACHIEVE KRANKE BY KRANKE UPPER
C     TRIANGULAR FORM.
      IF (.NOT.(KRANKE.GT.0 .AND. KRANKE.LT.ME)) GO TO 220
      DO 210 KK=1,KRANKE
        K = KRANKE + 1 - KK
C
C     APPLY TRANFORMATION TO MATRIX COLS. 1,...,K-1.
        CALL DH12(1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW, K-1)
C
C     APPLY TO RT SIDE VECTOR.
        CALL DH12(2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1), 1, 1, 1)
  210 CONTINUE
  220 IF (.NOT.(KRANKE.GT.0)) GO TO 240
C
C     SOLVE FOR VARIABLES 1,...,KRANKE IN NEW COORDINATES.
      CALL DCOPY(KRANKE, W(1,NP1), 1, X, 1)
      DO 230 I=1,KRANKE
        X(I) = (X(I)-DDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  230 CONTINUE
C
C     COMPUTE RESIDUALS FOR REDUCED PROBLEM.
  240 MEP1 = ME + 1
      RNORML = ZERO
      IF (.NOT.(ME.LT.M)) GO TO 270
      DO 260 I=MEP1,M
        W(I,NP1) = W(I,NP1) - DDOT(KRANKE,W(I,1),MDW,X,1)
        SN = DDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
        RN = DDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
        IF (.NOT.(RN.LE.TAU**2*SN .AND. KRANKE.LT.N)) GO TO 250
        W(I,KRANKE+1) = ZERO
        CALL DCOPY(N-KRANKE, W(I,KRANKE+1), 0, W(I,KRANKE+1), MDW)
  250   CONTINUE
  260 CONTINUE
C
C     COMPUTE EQUAL. CONSTRAINT EQUAS. RESIDUAL LENGTH.
  270 RNORME = DNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
C
C     MOVE REDUCED PROBLEM DATA UPWARD IF KRANKE.LT.ME.
      IF (.NOT.(KRANKE.LT.ME)) GO TO 290
      DO 280 J=1,NP1
        CALL DCOPY(M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  280 CONTINUE
C
C     COMPUTE SOLN OF REDUCED PROBLEM.
  290 CALL DLSI(W(KRANKE+1,KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,
     1 X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
      IF (.NOT.(ME.GT.0)) GO TO 330
C
C     TEST FOR CONSISTENCY OF EQUALITY CONSTRAINTS.
      MDEQC = 0
      XNRME = DASUM(KRANKE,W(1,NP1),1)
      IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
      MODE = MODE + MDEQC
C
C     CHECK IF SOLN TO EQUAL. CONSTRAINTS SATISFIES INEQUAL.
C     CONSTRAINTS WHEN THERE ARE NO DEGREES OF FREEDOM LEFT.
      IF (.NOT.(KRANKE.EQ.N .AND. MG.GT.0)) GO TO 320
      XNORM = DASUM(N,X,1)
      MAPKE1 = MA + KRANKE + 1
      MEND = MA + KRANKE + MG
      DO 310 I=MAPKE1,MEND
        SIZE = DASUM(N,W(I,1),MDW)*XNORM + DABS(W(I,NP1))
        IF (.NOT.(W(I,NP1).GT.TAU*SIZE)) GO TO 300
        MODE = MODE + 2
        GO TO 450
  300   CONTINUE
  310 CONTINUE
  320 CONTINUE
  330 IF (.NOT.(KRANKE.GT.0)) GO TO 420
C
C     REPLACE DIAG. TERMS OF LOWER TRAP. MATRIX.
      CALL DCOPY(KRANKE, WS(KRANKE+1), 1, W, MDW+1)
C
C     REAPPLY TRANS TO PUT SOLN IN ORIGINAL COORDINATES.
      DO 340 II=1,KRANKE
        I = KRANKE + 1 - II
        CALL DH12(2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  340 CONTINUE
C
C     COMPUTE COV MATRIX OF EQUAL. CONSTRAINED PROBLEM.
      IF (.NOT.(COV)) GO TO 410
      DO 400 JJ=1,KRANKE
        J = KRANKE + 1 - JJ
        IF (.NOT.(J.LT.N)) GO TO 390
        RB = WS(J)*W(J,J)
        IF (RB.NE.ZERO) RB = ONE/RB
        JP1 = J + 1
        DO 350 I=JP1,N
          W(I,J) = DDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)*RB
  350   CONTINUE
        GAM = DDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)*RB
        GAM = HALF*GAM
        CALL DAXPY(N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
        DO 370 I=JP1,N
          DO 360 K=I,N
            W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
            W(K,I) = W(I,K)
  360     CONTINUE
  370   CONTINUE
        UJ = WS(J)
        VJ = GAM*UJ
        W(J,J) = UJ*VJ + UJ*VJ
        DO 380 I=JP1,N
          W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  380   CONTINUE
        CALL DCOPY(N-J, W(J,JP1), MDW, W(JP1,J), 1)
  390   CONTINUE
  400 CONTINUE
  410 CONTINUE
C
C     APPLY THE SCALING TO THE COVARIANCE MATRIX.
  420 IF (.NOT.(COV)) GO TO 440
      DO 430 I=1,N
        L = N1 + I
        CALL DSCAL(N, WS(L-1), W(I,1), MDW)
        CALL DSCAL(N, WS(L-1), W(1,I), 1)
  430 CONTINUE
  440 CONTINUE
  450 CONTINUE
C
C     RESCALE SOLN. VECTOR.
      IF (.NOT.(MODE.LE.1)) GO TO 470
      DO 460 J=1,N
        L = N1 + J
        X(J) = X(J)*WS(L-1)
  460 CONTINUE
  470 IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
  480 CONTINUE
C     TO PROCESS-OPTION-VECTOR
C
C     THE NOMINAL TOLERANCE USED IN THE CODE
C     FOR THE EQUALITY CONSTRAINT EQUATIONS.
      TAU = DSQRT(DRELPR)
C
C     THE NOMINAL COLUMN SCALING USED IN THE CODE IS
C     THE IDENTITY SCALING.
      WS(N1) = ONE
      CALL DCOPY(N, WS(N1), 0, WS(N1), 1)
C
C     NO COVARIANCE MATRIX IS NOMINALLY COMPUTED.
      COV = .FALSE.
C
C     DEFINE BOUND FOR NUMBER OF OPTIONS TO CHANGE.
      NOPT = 1000
      NTIMES = 0
C
C     DEFINE BOUND FOR POSITIVE VALUES OF LINK.
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (.NOT.(LINK.LE.0 .OR. LINK.GT.NLINK)) GO TO 490
      NERR = 3
      IOPT = 1
      CALL XERROR('DLSEI( ) THE OPTION VECTOR IS UNDEFINED', 38, NERR,
     1 IOPT)
      MODE = 4
      RETURN
  490 IF (.NOT.(LINK.GT.1)) GO TO 540
      NTIMES = NTIMES + 1
      IF (.NOT.(NTIMES.GT.NOPT)) GO TO 500
      NERR = 3
      IOPT = 1
      CALL XERROR('DLSEI( ). THE LINKS IN THE OPTION VECTOR ARE CYCLING.
     1', 52,      NERR, IOPT)
      MODE = 4
      RETURN
  500 KEY = PRGOPT(LAST+1)
      IF (KEY.EQ.1) COV = PRGOPT(LAST+2).NE.ZERO
      IF (.NOT.(KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.ZERO)) GO TO 520
      DO 510 J=1,N
        T = DNRM2(M,W(1,J),1)
        IF (T.NE.ZERO) T = ONE/T
        L = N1 + J
        WS(L-1) = T
  510 CONTINUE
  520 IF (KEY.EQ.3) CALL DCOPY(N, PRGOPT(LAST+2), 1, WS(N1), 1)
      IF (KEY.EQ.4) TAU = DMAX1(DRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      IF (.NOT.(NEXT.LE.0 .OR. NEXT.GT.NLINK)) GO TO 530
      NERR = 3
      IOPT = 1
      CALL XERROR('DLSEI( ) THE OPTION VECTOR IS UNDEFINED', 38, NERR,
     1 IOPT)
      MODE = 4
      RETURN
  530 LAST = LINK
      LINK = NEXT
      GO TO 490
  540 DO 550 J=1,N
        L = N1 + J
        CALL DSCAL(M, WS(L-1), W(1,J), 1)
  550 CONTINUE
      GO TO 560
  560 GO TO IGO990, (90)
      END
*DECK DLSI
      SUBROUTINE DLSI(W,MDW,MA,MG,N,PRGOPT,X,RNORM,MODE,WS,IP)
C***BEGIN PROLOGUE  DLSI
C***REFER TO  DLSEI
C***ROUTINES CALLED  DASUM,DAXPY,DCOPY,DDOT,DH12,DHFTI,DLPDP,DSCAL,
C                    DSWAP
C***DESCRIPTION
C
C  **** Double Precision version of LSI ****
C     This is a companion subprogram to DLSEI( ).
C     The documentation for DLSEI( ) has more complete
C     usage instructions.
C     Written by R. J. Hanson, SLA.
C
C     Solve..
C              AX = B,  A  MA by N  (least squares equations)
C     subject to..
C
C              GX.GE.H, G  MG by N  (inequality constraints)
C
C     Input..
C
C      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
C                       (G H)
C
C     MDW,MA,MG,N
C              contain (resp) var. dimension of W(*,*),
C              and matrix dimensions.
C
C     PRGOPT(*),
C              Program option vector.
C
C     OUTPUT..
C
C      X(*),RNORM
C
C              Solution vector(unless MODE=2), length of AX-B.
C
C      MODE
C              =0   Inequality constraints are compatible.
C              =2   Inequality constraints contradictory.
C
C      WS(*),
C              Working storage of dimension K+N+(MG+2)*(N+7),
C              where K=MAX(MA+MG,N).
C      IP(MG+2*N+1)
C              Integer working storage
C
C     Revised April 16, 1982.
C***END PROLOGUE  DLSI
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     (START EDITING AT LINE WITH C++ IN COLS. 1-3.)
C     /REAL (12 BLANKS)/DOUBLE PRECISION/,/DASUM/DASUM/,/DDOT/DDOT/,
C     / SQRT/ DSQRT/,/AMAX1/DMAX1/,/DSWAP/DSWAP/,
C     /DCOPY/DCOPY/,/DSCAL/DSCAL/,/DAXPY/DAXPY/,/E0/D0/,/SRELPR/DRELPR/
C
C     SUBROUTINES CALLED
C
C     DLPDP          THIS SUBPROGRAM MINIMIZES A SUM OF SQUARES
C                   OF UNKNOWNS SUBJECT TO LINEAR INEQUALITY
C                   CONSTRAINTS.  PART OF THIS PACKAGE.
C
C++
C     DDOT,DSCAL    SUBROUTINES FROM THE BLAS PACKAGE.
C     DAXPY,DASUM,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     DCOPY,DSWAP
C
C     DHFTI          SOLVES AN UNCONSTRAINED LINEAR LEAST SQUARES
C                   PROBLEM.  PART OF THIS PACKAGE.
C
C     DH12           SUBROUTINE TO CONSTRUCT AND APPLY A HOUSEHOLDER
C                   TRANSFORMATION.
C
C     SUBROUTINE DLSI(W,MDW,MA,MG,N,PRGOPT,X,RNORM,MODE,WS,IP)
C
      DOUBLE PRECISION W(MDW,*), PRGOPT(*), RNORM, WS(*), X(*)
      INTEGER IP(*)
      DOUBLE PRECISION ANORM, DRELPR, FAC, GAM, HALF, ONE, RB, TAU, TOL,
     1 XNORM, ZERO
      DOUBLE PRECISION DASUM, DDOT, DSQRT, DMAX1
      LOGICAL COV, SCLCOV
C
      DATA ZERO /0.D0/, DRELPR /0.D0/, ONE /1.D0/, HALF /.5D0/
C
C     COMPUTE MACHINE PRECISION=DRELPR ONLY WHEN NECESSARY.
C***FIRST EXECUTABLE STATEMENT  DLSI
      IF (.NOT.(DRELPR.EQ.ZERO)) GO TO 30
      DRELPR = ONE
   10 IF (ONE+DRELPR.EQ.ONE) GO TO 20
      DRELPR = DRELPR*HALF
      GO TO 10
   20 DRELPR = DRELPR + DRELPR
   30 MODE = 0
      RNORM = ZERO
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 70
      ASSIGN 40 TO IGO994
      GO TO 500
C
C     PROCESS-OPTION-VECTOR
C
C     COMPUTE MATRIX NORM OF LEAST SQUARES EQUAS.
   40 ANORM = ZERO
      DO 50 J=1,N
        ANORM = DMAX1(ANORM,DASUM(MA,W(1,J),1))
   50 CONTINUE
C
C     SET TOL FOR DHFTI( ) RANK TEST.
      TAU = TOL*ANORM
C
C     COMPUTE HOUSEHOLDER ORTHOGONAL DECOMP OF MATRIX.
      IF (N.GT.0) WS(1) = ZERO
      CALL DCOPY(N, WS, 0, WS, 1)
      CALL DCOPY(MA, W(1,NP1), 1, WS, 1)
      K = MAX0(M,N)
      MINMAN = MIN0(MA,N)
      N1 = K + 1
      N2 = N1 + N
      CALL DHFTI(W, MDW, MA, N, WS, 1, 1, TAU, KRANK, RNORM, WS(N2),
     1 WS(N1), IP)
      FAC = ONE
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
      ASSIGN 60 TO IGO990
      GO TO 80
C
C     REDUCE-TO-DLPDP-AND-SOLVE
   60 CONTINUE
   70 IP(1) = KRANK
      IP(2) = N + MAX0(M,N) + (MG+2)*(N+7)
      RETURN
   80 CONTINUE
C
C     TO REDUCE-TO-DLPDP-AND-SOLVE
      MAP1 = MA + 1
C
C     COMPUTE INEQ. RT-HAND SIDE FOR DLPDP.
      IF (.NOT.(MA.LT.M)) GO TO 260
      IF (.NOT.(MINMAN.GT.0)) GO TO 160
      DO 90 I=MAP1,M
        W(I,NP1) = W(I,NP1) - DDOT(N,W(I,1),MDW,WS,1)
   90 CONTINUE
      DO 100 I=1,MINMAN
        J = IP(I)
C
C     APPLY PERMUTATIONS TO COLS OF INEQ. CONSTRAINT MATRIX.
        CALL DSWAP(MG, W(MAP1,I), 1, W(MAP1,J), 1)
  100 CONTINUE
C
C     APPLY HOUSEHOLDER TRANSFORMATIONS TO CONSTRAINT MATRIX.
      IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.N)) GO TO 120
      DO 110 II=1,KRANK
        I = KRANK + 1 - II
        L = N1 + I
        CALL DH12(2, I, KRANK+1, N, W(I,1), MDW, WS(L-1), W(MAP1,1),
     1   MDW, 1, MG)
  110 CONTINUE
C
C     COMPUTE PERMUTED INEQ. CONSTR. MATRIX TIMES R-INVERSE.
  120 DO 150 I=MAP1,M
        IF (.NOT.(0.LT.KRANK)) GO TO 140
        DO 130 J=1,KRANK
          W(I,J) = (W(I,J)-DDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  130   CONTINUE
  140   CONTINUE
  150 CONTINUE
C
C     SOLVE THE REDUCED PROBLEM WITH DLPDP ALGORITHM,
C     THE LEAST PROJECTED DISTANCE PROBLEM.
  160 CALL DLPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X, XNORM,
     1 MDLPDP, WS(N2), IP(N+1))
      IF (.NOT.(MDLPDP.EQ.1)) GO TO 240
      IF (.NOT.(KRANK.GT.0)) GO TO 180
C
C     COMPUTE SOLN IN ORIGINAL COORDINATES.
      DO 170 II=1,KRANK
        I = KRANK + 1 - II
        X(I) = (X(I)-DDOT(II-1,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170 CONTINUE
C
C     APPLY HOUSEHOLDER TRANS. TO SOLN VECTOR.
  180 IF (.NOT.(0.LT.KRANK .AND. KRANK.LT.N)) GO TO 200
      DO 190 I=1,KRANK
        L = N1 + I
        CALL DH12(2, I, KRANK+1, N, W(I,1), MDW, WS(L-1), X, 1, 1, 1)
  190 CONTINUE
  200 IF (.NOT.(MINMAN.GT.0)) GO TO 230
C
C     REPERMUTE VARIABLES TO THEIR INPUT ORDER.
      DO 210 II=1,MINMAN
        I = MINMAN + 1 - II
        J = IP(I)
        CALL DSWAP(1, X(I), 1, X(J), 1)
  210 CONTINUE
C
C     VARIABLES ARE NOW IN ORIG. COORDINATES.
C     ADD SOLN OF UNSCONSTRAINED PROB.
      DO 220 I=1,N
        X(I) = X(I) + WS(I)
  220 CONTINUE
C
C     COMPUTE THE RESIDUAL VECTOR NORM.
      RNORM = DSQRT(RNORM**2+XNORM**2)
  230 GO TO 250
  240 MODE = 2
  250 GO TO 270
  260 CALL DCOPY(N, WS, 1, X, 1)
  270 IF (.NOT.(COV .AND. KRANK.GT.0)) GO TO 490
C
C     COMPUTE COVARIANCE MATRIX BASED ON THE ORTHOGONAL DECOMP.
C     FROM DHFTI( ).
C
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
C
C     COPY DIAG. TERMS TO WORKING ARRAY.
      CALL DCOPY(KRANK, W, MDW+1, WS(N2), 1)
C
C     RECIPROCATE DIAG. TERMS.
      DO 280 J=1,KRANK
        W(J,J) = ONE/W(J,J)
  280 CONTINUE
      IF (.NOT.(KRANK.GT.1)) GO TO 310
C
C     INVERT THE UPPER TRIANGULAR QR FACTOR ON ITSELF.
      DO 300 I=1,KRM1
        IP1 = I + 1
        DO 290 J=IP1,KRANK
          W(I,J) = -DDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  290   CONTINUE
  300 CONTINUE
C
C     COMPUTE THE INVERTED FACTOR TIMES ITS TRANSPOSE.
  310 DO 330 I=1,KRANK
        DO 320 J=I,KRANK
          W(I,J) = DDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  320   CONTINUE
  330 CONTINUE
      IF (.NOT.(KRANK.LT.N)) GO TO 450
C
C     ZERO OUT LOWER TRAPEZOIDAL PART.
C     COPY UPPER TRI. TO LOWER TRI. PART.
      DO 340 J=1,KRANK
        CALL DCOPY(J, W(1,J), 1, W(J,1), MDW)
  340 CONTINUE
      DO 350 I=KRP1,N
        W(I,1) = ZERO
        CALL DCOPY(I, W(I,1), 0, W(I,1), MDW)
  350 CONTINUE
C
C     APPLY RIGHT SIDE TRANSFORMATIONS TO LOWER TRI.
      N3 = N2 + KRP1
      DO 430 I=1,KRANK
        L = N1 + I
        K = N2 + I
        RB = WS(L-1)*WS(K-1)
        IF (.NOT.(RB.LT.ZERO)) GO TO 420
C
C     IF RB.GE.ZERO, TRANSFORMATION CAN BE REGARDED AS ZERO.
        RB = ONE/RB
C
C     STORE UNSCALED RANK-ONE HOUSEHOLDER UPDATE IN WORK ARRAY.
        WS(N3) = ZERO
        CALL DCOPY(N, WS(N3), 0, WS(N3), 1)
        L = N1 + I
        K = N3 + I
        WS(K-1) = WS(L-1)
        DO 360 J=KRP1,N
          K = N3 + J
          WS(K-1) = W(I,J)
  360   CONTINUE
        DO 370 J=1,N
          L = N3 + I
          K = N3 + J
          WS(J) = DDOT(J-I,W(J,I),MDW,WS(L-1),1) + DDOT(N-J+1,W(J,J),1,
     1     WS(K-1),1)
          WS(J) = WS(J)*RB
  370   CONTINUE
        L = N3 + I
        GAM = DDOT(N-I+1,WS(L-1),1,WS(I),1)*RB
        GAM = GAM*HALF
        CALL DAXPY(N-I+1, GAM, WS(L-1), 1, WS(I), 1)
        DO 410 J=I,N
          IF (.NOT.(I.GT.1)) GO TO 390
          IM1 = I - 1
          K = N3 + J
          DO 380 L=1,IM1
            W(J,L) = W(J,L) + WS(K-1)*WS(L)
  380     CONTINUE
  390     K = N3 + J
          DO 400 L=I,J
            IL = N3 + L
            W(J,L) = W(J,L) + WS(J)*WS(IL-1) + WS(L)*WS(K-1)
  400     CONTINUE
  410   CONTINUE
  420   CONTINUE
  430 CONTINUE
C
C     COPY LOWER TRI. TO UPPER TRI. TO SYMMETRIZE THE COVARIANCE MATRIX.
      DO 440 I=1,N
        CALL DCOPY(I, W(I,1), MDW, W(1,I), 1)
  440 CONTINUE
C
C     REPERMUTE ROWS AND COLS.
  450 DO 470 II=1,MINMAN
        I = MINMAN + 1 - II
        K = IP(I)
        IF (.NOT.(I.NE.K)) GO TO 460
        CALL DSWAP(1, W(I,I), 1, W(K,K), 1)
        CALL DSWAP(I-1, W(1,I), 1, W(1,K), 1)
        CALL DSWAP(K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
        CALL DSWAP(N-K, W(I,K+1), MDW, W(K,K+1), MDW)
  460   CONTINUE
  470 CONTINUE
C
C     PUT IN NORMALIZED RESIDUAL SUM OF SQUARES SCALE FACTOR
C     AND SYMMETRIZE THE RESULTING COVARIANCE MARIX.
      DO 480 J=1,N
        CALL DSCAL(J, FAC, W(1,J), 1)
        CALL DCOPY(J, W(1,J), 1, W(J,1), MDW)
  480 CONTINUE
  490 GO TO 540
  500 CONTINUE
C
C     TO PROCESS-OPTION-VECTOR
C
C     THE NOMINAL TOLERANCE USED IN THE CODE,
      TOL = DSQRT(DRELPR)
      COV = .FALSE.
      SCLCOV = .TRUE.
      LAST = 1
      LINK = PRGOPT(1)
  510 IF (.NOT.(LINK.GT.1)) GO TO 520
      KEY = PRGOPT(LAST+1)
      IF (KEY.EQ.1) COV = PRGOPT(LAST+2).NE.ZERO
      IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2).EQ.ZERO
      IF (KEY.EQ.5) TOL = DMAX1(DRELPR,PRGOPT(LAST+2))
      NEXT = PRGOPT(LINK)
      LAST = LINK
      LINK = NEXT
      GO TO 510
  520 GO TO 530
  530 GO TO IGO994, (40)
  540 GO TO IGO990, (60)
      END
*DECK DLPDP
      SUBROUTINE DLPDP(A,MDA,M,N1,N2,PRGOPT,X,WNORM,MODE,WS,IS)
C***BEGIN PROLOGUE  DLPDP
C***REFER TO  DLSEI
C***ROUTINES CALLED  DCOPY,DDOT,DNRM2,DSCAL,DWNNLS
C***DESCRIPTION
C
C  **** Double Precision version of LPDP ****
C     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
C     where N=N1+N2.  This is a slight overestimate for WS(*).
C
C     Written by R. J. Hanson and K. H. Haskell, Sandia Labs
C     Revised Oct. 1, 1981.
C
C     Determine an N1-vector W, and
C               an N2-vector Z
C     which minimizes the Euclidean length of W
C     subject to G*W+H*Z .GE. Y.
C     This is the least projected distance problem, LPDP.
C     The matrices G and H are of respective
C     dimensions M by N1 and M by N2.
C
C     Called by subprogram DLSI( ).
C
C     The matrix
C                (G H Y)
C
C     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
C
C     The solution (W) is returned in X(*).
C                  (Z)
C
C     The value of MODE indicates the status of
C     the computation after returning to the user.
C
C          MODE=1  The solution was successfully obtained.
C
C          MODE=2  The inequalities are inconsistent.
C***END PROLOGUE  DLPDP
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /DOUBLE PRECISION (12 BLANKS)/DOUBLE
C     PRECISION/,/DNRM2/DNRM2/,/DDOT/DDOT/,
C     /DCOPY/DCOPY/,/DSCAL/DSCAL/,/DABS(/DABS(/, DABS/, DABS,/,/E0/D0/
C
      INTEGER I, IS(*), IW, IX, J, L, M, MDA, MODE, MODEW, N, N1, N2,
     *     NP1
      DOUBLE PRECISION A(MDA,*), DABS, DDOT, DNRM2, FAC, ONE,
     *     PRGOPT(*), RNORM, SC, WNORM, WS(*), X(*), YNORM, ZERO
      DATA ZERO,ONE /0.0D0,1.0D0/, FAC /0.1D0/
C***FIRST EXECUTABLE STATEMENT  DLPDP
      N = N1 + N2
      MODE = 1
      IF (M .GT. 0) GO TO 20
         IF (N .LE. 0) GO TO 10
            X(1) = ZERO
            CALL DCOPY(N,X,0,X,1)
   10    CONTINUE
         WNORM = ZERO
      GO TO 200
   20 CONTINUE
C        BEGIN BLOCK PERMITTING ...EXITS TO 190
            NP1 = N + 1
C
C           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
            DO 40 I = 1, M
               SC = DNRM2(N,A(I,1),MDA)
               IF (SC .EQ. ZERO) GO TO 30
                  SC = ONE/SC
                  CALL DSCAL(NP1,SC,A(I,1),MDA)
   30          CONTINUE
   40       CONTINUE
C
C           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
            YNORM = DNRM2(M,A(1,NP1),1)
            IF (YNORM .EQ. ZERO) GO TO 50
               SC = ONE/YNORM
               CALL DSCAL(M,SC,A(1,NP1),1)
   50       CONTINUE
C
C           SCALE COLS OF MATRIX H.
            J = N1 + 1
   60       IF (J .GT. N) GO TO 70
               SC = DNRM2(M,A(1,J),1)
               IF (SC .NE. ZERO) SC = ONE/SC
               CALL DSCAL(M,SC,A(1,J),1)
               X(J) = SC
               J = J + 1
            GO TO 60
   70       CONTINUE
            IF (N1 .LE. 0) GO TO 130
C
C              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
               IW = 0
               DO 80 I = 1, M
C
C                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
                  CALL DCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
C
C                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
                  CALL DCOPY(N1,A(I,1),MDA,WS(IW+1),1)
                  IW = IW + N1
C
C                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
   80          CONTINUE
               WS(IW+1) = ZERO
               CALL DCOPY(N,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N
               WS(IW+1) = ONE
               IW = IW + 1
C
C              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
C              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
C              F = TRANSPOSE OF (0,...,0,1).
               IX = IW + 1
               IW = IW + M
C
C              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
C              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,NP1,N2,NP1-N2,M,0,PRGOPT,WS(IX),RNORM,
     *                     MODEW,IS,WS(IW+1))
C
C              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
               SC = ONE - DDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*DABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)
     *            GO TO 110
                  SC = ONE/SC
                  DO 90 J = 1, N1
                     X(J) = SC*DDOT(M,A(1,J),1,WS(IX),1)
   90             CONTINUE
C
C                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
C                 VECTOR.
                  DO 100 I = 1, M
                     A(I,NP1) = A(I,NP1) - DDOT(N1,A(I,1),MDA,X,1)
  100             CONTINUE
               GO TO 120
  110          CONTINUE
                  MODE = 2
C        .........EXIT
                  GO TO 190
  120          CONTINUE
  130       CONTINUE
            IF (N2 .LE. 0) GO TO 180
C
C              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
               IW = 0
               DO 140 I = 1, M
                  CALL DCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
  140          CONTINUE
               WS(IW+1) = ZERO
               CALL DCOPY(N2,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N2
               WS(IW+1) = ONE
               IW = IW + 1
               IX = IW + 1
               IW = IW + M
C
C              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
C              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
C              OF (0,...,0,1)).
C
C              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
C              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,N2+1,0,N2+1,M,0,PRGOPT,WS(IX),RNORM,MODEW,
     *                     IS,WS(IW+1))
C
C              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
               SC = ONE - DDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*DABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)
     *            GO TO 160
                  SC = ONE/SC
                  DO 150 J = 1, N2
                     L = N1 + J
                     X(L) = SC*DDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150             CONTINUE
               GO TO 170
  160          CONTINUE
                  MODE = 2
C        .........EXIT
                  GO TO 190
  170          CONTINUE
  180       CONTINUE
C
C           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
            CALL DSCAL(N,YNORM,X,1)
            WNORM = DNRM2(N1,X,1)
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE ZEROIN(F,B,C,RE,AE,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2646
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980
C                   *************************
C                   *       ISSUED BY       *
C                   *  SANDIA LABORATORIES, *
C                   *   A PRIME CONTRACTOR  *
C                   ********     TO THE     *
C                          *  UNITED STATES *
C                          *   DEPARTMENT   *
C                          *       OF       *
C                          *     ENERGY     *
C      *********************  ---NOTICE---  *********************
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *
C      *               NOR ANY OF THEIR EMPLOYEES,              *
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *
C      *          **********    ACCURACY,   **********          *
C      *          *        *  COMPLETENESS  *        *          *
C      *          *        *  OR USEFULNESS *        *          *
C      *          *        *     OF ANY     *        *          *
C      *          *        *  INFORMATION,  *        *          *
C      *          *        *   APPARATUS,   *        *          *
C      *       ****        *     PRODUCT    *        ****       *
C      *       *           *   OR PROCESS   *           *       *
C      *       *           *   DISCLOSED,   *           *       *
C      *       *           *  OR REPRESENTS *           *       *
C      *       *          **    THAT ITS    **          *       *
C      *       *          **  USE WOULD NOT **          *       *
C      *********          **    INFRINGE    **          *********
C                         **    PRIVATELY   **
C                         **      OWNED     **
C                         **     RIGHTS.    **
C                         **                **
C                         **                **
C                         **                **
C                         ********************
C
C     BASED ON A METHOD BY T J DEKKER
C     WRITTEN BY L F SHAMPINE AND H A WATTS
C     MODIFIED FOR THE MATH LIBRARY BY C B BAILEY
C
C     ABSTRACT
C        ZEROIN SEARCHES FOR A ZERO OF A FUNCTION F(X) BETWEEN
C        THE GIVEN VALUES B AND C UNTIL THE WIDTH OF THE INTERVAL
C        (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE SPECIFIED BY
C        THE STOPPING CRITERION, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C        THE METHOD USED IS AN EFFICIENT COMBINATION OF BISECTION AND
C        THE SECANT RULE.  IN ORDER TO INSURE THAT ZEROIN WILL CONVERGE
C        TO A ZERO, THE USER SHOULD PICK VALUES FOR B AND C AT WHICH
C        THE FUNCTION DIFFERS IN SIGN.
C
C     DESCRIPTION OF ARGUMENTS
C     F,B,C,RE AND AE ARE INPUT PARAMETERS
C     B,C AND IFLAG ARE OUTPUT PARAMETERS
C        F     - NAME OF THE REAL VALUED EXTERNAL FUNCTION.  THIS NAME
C                MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING
C                PROGRAM.  F MUST BE A FUNCTION OF ONE REAL ARGUMENT.
C        B     - ONE END OF THE INTERVAL (B,C).  THE VALUE RETURNED FOR
C                B USUALLY IS THE BETTER APPROXIMATION TO A ZERO OF F.
C        C     - THE OTHER END OF THE INTERVAL (B,C)
C        RE    - RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
C                IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
C                THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
C        AE    - ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
C                GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
C                NONZERO VALUE SHOULD BE CHOSEN FOR AE.
C        IFLAG - A STATUS CODE.  USER MUST CHECK IFLAG AFTER EACH CALL.
C                CONTROL RETURNS TO THE USER FROM ZEROIN IN ALL CASES.
C                XERROR DOES NOT PROCESS DIAGNOSTICS IN THESE CASES.
C                 1 B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C), AND
C                   F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C                 2 F(B) = 0.  HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
C                   COLLAPSED TO THE REQUESTED TOLERANCE.
C                 3 B MAY BE NEAR A SINGULAR POINT OF F(X).
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE AND THE FUNCTION CHANGES SIGN IN (B,C) BUT
C                   F(X) INCREASED IN MAGNITUDE AS (B,C) COLLAPSED,I.E.
C                     ABS(F(B OUT)) .GT. MAX(ABS(F(B IN)),ABS(F(C IN)))
C                 4 NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
C                   INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
C                   THE USER MUST EXAMINE THIS CASE AND DECIDE WHETHER
C                   B IS NEAR A LOCAL MINIMUM OF F(X), OR B IS NEAR A
C                   ZERO OF EVEN MULTIPLICITY, OR NEITHER OF THESE.
C                 5 TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.
C
C     REFERENCES
C       1.  L F SHAMPINE AND H A WATTS, ZEROIN, A ROOT-SOLVING CODE,
C           SC-TM-70-631, SEPT 1970.
C       2.  T J DEKKER, FINDING A ZERO BY MEANS OF SUCCESSIVE LINEAR
C           INTERPOLATION, *CONSTRUCTIVE ASPECTS OF THE FUNDAMENTAL
C           THEOREM OF ALGEBRA*, EDITED BY B DEJON AND P HENRICI, 1969.
C
C
C     ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS
C     DEFINED HERE BY THE FUNCTION D1MACH.
C
      ER = 2.0D0 * D1MACH(4)
C
C     INITIALIZE
      RW=DMAX1(RE,ER)
      AW=DMAX1(AE,0.0D0)
      IC=0
      ACBS=DABS(B-C)
      A=C
      T=A
      FA=F(T)
      T=B
      FB=F(T)
      FC=FA
      KOUNT=2
      FX=DMAX1(DABS(FB),DABS(FC))
C
    1 IF (DABS(FC) .GE. DABS(FB)) GO TO 2
C     PERFORM INTERCHANGE
      A=B
      FA=FB
      B=C
      FB=FC
      C=A
      FC=FA
C
    2 IF (FB .EQ. 0.0D0) GO TO 11
      CMB=0.5D0*(C-B)
      ACMB=DABS(CMB)
      TOL=RW*DABS(B)+AW
C
C     TEST STOPPING CRITERION
      IF (ACMB .LE. TOL) GO TO 10
C
C     CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
C     WHERE WE ARRANGE P .GE. 0.
C     THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
      P=(B-A)*FB
      Q=FA-FB
      IF (P .GE. 0.0D0) GO TO 3
      P=-P
      Q=-Q
C
C     UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
C     IN THE SIZE OF OUR BOUNDING INTERVAL.
    3 A=B
      FA=FB
      IC=IC+1
      IF (IC .LT. 4) GO TO 4
      IF (8.0D0*ACMB .GE. ACBS) GO TO 6
      IC=0
      ACBS=ACMB
C
C     TEST FOR TOO SMALL A CHANGE
    4 IF (P .GT. DABS(Q)*TOL) GO TO 5
C
C     INCREMENT BY TOLERANCE
      B=B+DSIGN(TOL,CMB)
      GO TO 7
C
C     ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.0D0
    5 IF (P .GE. CMB*Q) GO TO 6
C
C     INTERPOLATE
      B=B+P/Q
      GO TO 7
C
    6 B=0.5D0*(C+B)
C     BISECT
C
C     HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
    7 T=B
      FB=F(T)
      IF (FB .EQ. 0.0D0) GO TO 11
C
C     DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
      IF (DSIGN(1.0D0,FB) .NE. DSIGN(1.0D0,FC)) GO TO 8
      C=A
      FC=FA
    8 KOUNT=KOUNT+1
      IF (KOUNT .GT. 500) GO TO 15
      GO TO 1
C
C
C     FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
C
   10 IF (DSIGN(1.0D0,FB) .EQ. DSIGN(1.0D0,FC)) GO TO 13
      IF (DABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN
   11 IFLAG = 2
      RETURN
   12 IFLAG = 3
      RETURN
   13 IFLAG = 4
      RETURN
   15 IFLAG = 5
      RETURN
      END
