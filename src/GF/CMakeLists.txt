FILE(GLOB GF_LIB_SRC
  GF_ELECTRONS.f90
  GREENS_FUNCTIONS.f90
  )
ADD_LIBRARY(GF_LIB OBJECT  ${GF_LIB_SRC})


