#ifdef __DEBUG__
#define __DECL_CLASS_UNUSED_THIS__ \
  LOGICAL :: DUMMY_LOGICAL=.TRUE.;
#else
#define __DECL_CLASS_UNUSED_THIS__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_CLASS_UNUSED_THIS__ \
  DUMMY_LOGICAL=THIS%DUMMY_LOGICAL;
#else
#define __SUPPRESS_CLASS_UNUSED_THIS__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_CLASS_UNUSED__(C) \
 C%DUMMY_LOGICAL=.TRUE.
#else
#define __SUPPRESS_CLASS_UNUSED__(C)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_FILE_READER__ \
  INTEGER :: DUMMY_INTEGER_FR;
#else
#define __DECL_UNUSED_FILE_READER__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_FILE_READER__(FR) \
  DUMMY_INTEGER_FR=FR%getFileID();
#else
#define __SUPPRESS_UNUSED_FILE_READER__(FR)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_CHARACTER__ \
  CHARACTER :: DUMMY_CHARACTER;
#else
#define __DECL_UNUSED_CHARACTER__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_STRING__(S) \
DUMMY_CHARACTER=S(1:1);
#else
#define __SUPPRESS_UNUSED_STRING__(S)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_REAL__ \
  REAL(k_real) :: DUMMY_REAL;
#else
#define __DECL_UNUSED_REAL__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_REAL__(R) \
  DUMMY_REAL=R;
#else
#define __SUPPRESS_UNUSED_REAL__(R)
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_REAL_OUT__(R) \
  R=1.;
#else
#define __SUPPRESS_UNUSED_REAL_OUT__(R)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_INTEGER__ \
  INTEGER :: DUMMY_INT;
#else
#define __DECL_UNUSED_INTEGER__
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_LOGICAL__ \
  LOGICAL :: DUMMY_LOGICAL_2;
#else
#define __DECL_UNUSED_LOGICAL__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_INTEGER__(I) \
  DUMMY_INT=I;
#else
#define __SUPPRESS_UNUSED_INTEGER__(I)
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_LOGICAL__(I) \
  DUMMY_LOGICAL_2=I;
#else
#define __SUPPRESS_UNUSED_LOGICAL__(I)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_MATRIX_PTR__ \
  real(k_real), pointer, dimension(:,:) :: UNUSED_MATRIX_PTR=>null()
#else
#define __DECL_UNUSED_MATRIX_PTR__
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_MATRIX_WARNING__(M) \
  UNUSED_MATRIX_PTR=>M;
#else
#define __SUPPRESS_UNUSED_MATRIX_WARNING__(M)
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_VECTOR_PTR__ \
  real(k_real), pointer, dimension(:) :: UNUSED_VECTOR_PTR=>null()
#else
#define __DECL_UNUSED_VECTOR_PTR__
#endif

#ifdef __DEBUG__
#define __DECL_UNUSED_VECTOR__(N) \
  real(k_real) :: UNUSED_VECTOR(N)
#else
#define __DECL_UNUSED_VECTOR__(N)
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_VECTOR__(V) \
  UNUSED_VECTOR=V;
#else
#define __SUPPRESS_UNUSED_VECTOR__(V)
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_VECTOR_WARNING__(V) \
  UNUSED_VECTOR_PTR=>V;
#else
#define __SUPPRESS_UNUSED_VECTOR_WARNING__(V)
#endif

#ifdef __DEBUG__
#define __SUPPRESS_UNUSED_LOGICAL_OUT__(L) \
  L=.TRUE.;
#else
#define __SUPPRESS_UNUSED_LOGICAL_OUT__(L)
#endif