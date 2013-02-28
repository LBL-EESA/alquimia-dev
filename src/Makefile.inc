#
# setup the external library dependancies
#

include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

CFLAGS += -g -W -Wall -std=c99 -pedantic -Werror
CXXFLAGS += -g -W -Wall
FCFLAGS += -g -W -Wall -std=f2003 -pedantic -Wno-unused-variable -Wno-unused-parameter


DEPEND_LIBS = 

ifdef PFLOTRAN_DIR
      CFLAGS += -DHAVE_PFLOTRAN
      CXXFLAGS += ${PETSC_CC_INCLUDES}
      FCFLAGS += -I$(PFLOTRAN_DIR)/src/pflotran ${PETSC_FC_INCLUDES}
      LDFLAGS += -L$(PFLOTRAN_DIR)/src/pflotran -lpflotranchem 
      LDFLAGS += ${PETSC_LIB} ${LIB}
#	DEPEND_LIBS += \
#		$(PFLOTRAN_DIR)/src/pflotran/libpflotranchem.a
endif

ifdef CRUNCH_DIR
	CXXFLAGS += -DHAVE_CRUNCH
#	DEPEND_LIBS += \
#		$(CRUNCH_DIR)/src/libcrunch_chem.a
#       LDFLAGS += -L$(CRUNCH_DIR)/src -lcrunch_chem
endif

ifdef DEBUG
      CFLAGS += -g
      CXXFLAGS += -g
      FCFLAGS += -g
endif
