#
# setup the external library dependancies
#
DEPEND_LIBS = 

ifdef PFLOTRAN_DIR
	CXXFLAGS += -DHAVE_PFLOTRAN
#	DEPEND_LIBS += \
#		$(PFLOTRAN_DIR)/src/pflotran/libpflotran_chem.a
#       LDFLAGS += -L$(PFLOTRAN_DIR)/src/pflotran -lpflotran_chem
endif

ifdef CRUNCH_DIR
	CXXFLAGS += -DHAVE_CRUNCH
#	DEPEND_LIBS += \
#		$(CRUNCH_DIR)/src/libcrunch_chem.a
#       LDFLAGS += -L$(CRUNCH_DIR)/src -lcrunch_chem
endif

ifdef DEBUG
      CXXFLAGS += -g
      FCFLAGS += -g
endif
