# Meschach matrix library (vendored; build originally by Petter Strandmark)
SET(MESCHACH_INTERNAL_SRC
  bdfactor.c
  bkpfacto.c
  chfactor.c
  copy.c
  err.c
  extras.c
  fft.c
  givens.c
  hessen.c
  hsehldr.c
  init.c
  iter0.c
  iternsym.c
  itersym.c
  ivecop.c
  lufactor.c
  machine.c
  matlab.c
  matop.c
  matrixio.c
  meminfo.c
  memory.c
  memstat.c
  mfunc.c
  norm.c
  otherio.c
  pxop.c
  qrfactor.c
  schur.c
  solve.c
  sparse.c
  sparseio.c
  spbkp.c
  spchfctr.c
  splufctr.c
  sprow.c
  spswap.c
  submat.c
  svd.c
  symmeig.c
  update.c
  vecop.c
  version.c
#  zcopy.c
#  zfunc.c
#  zgivens.c
#  zhessen.c
#  zhsehldr.c
#  zlufctr.c
#  zmachine.c
#  zmatio.c
#  zmatlab.c
#  zmatop.c
#  zmemory.c
#  znorm.c
#  zqrfctr.c
#  zschur.c
#  zsolve.c
#  ztorture.c
#  zvecop.c
)

add_library(meschach STATIC ${MESCHACH_INTERNAL_SRC})
target_include_directories(meschach PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(meschach PUBLIC THREADSAFE)
target_link_libraries(meschach PUBLIC mfreex_compat)
