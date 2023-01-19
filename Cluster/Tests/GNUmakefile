GNUSTEP_MAKEFILES = /usr/GNUstep/System/Library/Makefiles

include $(GNUSTEP_MAKEFILES)/common.make

TOOL_NAME = LogTest
LogTest_OBJC_FILES = source.m
LogTest_TOOL_LIBS += -lphaseobjc
LogTest_LIB_DIRS += -L/usr/local/lib
LogTest_INCLUDE_DIRS += -I../../src

include $(GNUSTEP_MAKEFILES)/tool.make
