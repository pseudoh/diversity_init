######################################################################
## DESCRIPTON ########################################################
# This Makefile supports the following builds for c++ projects:
# - release: builds the project with optimization enabled.
# - debug: builds the project with debug information.
# - test: builds and runs all tests.
# - doc: builds Doxygen documentation.
# - clean: deletes Debug, Documentation, Release and Test directories.
#
# It will automatically detect sources and calculate dependencies of ,
# provided that:
# (1) files have a .cpp extension iff they are c++ source files.
# (2) no source file has a double underscore "__" in its full
#       pathname.
# (3) the source files containing the main functions for the release,
#       debug and test builds are found directly in the ./src
#       directory.
# (4) all source files required for the release and debug builds,
#       aside those covered in (3), or those in linked libraries, are
#       found in the "src" subdirectory of either the project's root
#       directory or a directory listed in the LOCAL_LIBS variable.
# (5) all source files required for the test build, aside those
#       covered in (3) and (4), or those in linked libraries, are
#       found in the "test" subdirectory of either the project's root
#       directory or a directory listed in the LOCAL_LIBS variable.
# (6) Google's test and mocking frameworks, gtest and gmock, are
#       available. If their location is other than the root directory
#       of the project (as shown in the example directory structure
#       below), the variables GMOCK_DIR and GTEST_DIR must be
#       modified appropriately.
#
# The following is an example directory conforming to the above. Note
# that duplicate file names in different directories are permitted.
#
# +-- _googlemock
# |   +-- _include
# |   +-- _src
# +-- _googletest
# |   +-- _include
# |   +-- _src
# +-- _mylib1
# |   +-- _src
# |   |   +-- some_file.cpp
# |   +-- _test
# |   |   +-- some_file_test.cpp
# +-- _mylib2
# |   +-- _sublib1
# |   |   +-- _src
# |   |   |   +-- some_file.cpp
# |   |   +-- _test
# |   |   |   +-- some_file_test.cpp
# |   +-- _sublib2
# |   |   +-- _src
# |   |   |   +-- some_file.cpp
# |   |   +-- _test
# |   |   |   +-- some_file_test.cpp
# +-- _src
# |   +-- main.cpp
# |   +-- test_main.cpp
# |   +-- some_file.cpp
# +-- _test
# |   +-- some_file_test.cpp
# +-- Makefile
#
# In this example we would set
#   LOCAL_LIBS := mylib1 mylib2/sublib1 mylib2/sublib2

######################################################################
## OPTIONS ###########################################################
CXX := g++
CXXFLAGS :=  -Wextra -std=c++1y -Wno-unused-parameter #-Werror -Wall
LDFLAGS :=
LDLIBS :=
# Points to the root of Google Test and Google Mock, relative to one
# level up from the project root. Remember to tweak this if you move
# this file.
GTEST_DIR := ../googletest
GMOCK_DIR := ../googlemock
# LOCAL_LIBS is a list of paths from project base dir to library base
# dir. Each library is assumed to have:
#   - A src dir containing all relevant .cpp files
#   - A test dir containing the librariy's tests
LOCAL_LIBS := # TODO Complete
TARGET_FILENAME := us.o
RELEASE_MAIN := main.cpp
DEBUG_MAIN := main.cpp
TEST_MAIN := test_main.cpp

######################################################################
## UTILITY FUNCTIONS #################################################
# Creates a subdirectory and then makes a target in that directory
#  1 - make target
#  2 - directory name
make_target_in_subdir = mkdir -p $(2) &&\
	cp Makefile $(2)/Makefile &&\
	cd $(2) &&\
	make $(1)

#Transforms a path string to an object name, replacing "/" with "__".
#  1 - path string
path_to_object = $(patsubst ..__%.cpp, %.o, $(subst /,__,$(1)))

#Transforms an object name to a path string, replacing "__" with "/".
#  1 - object name
object_to_path = $(patsubst %.o, ../%.cpp, $(subst __,/,$(1)))

######################################################################
## BUILD TYPES #######################################################

## COMMON ############################################################
SRC_DIRS := ../src \
  $(addsuffix /src/, $(addprefix ../, $(LOCAL_LIBS)))
# Prevent multiple inclusions of main()
EXCLUDE_MAIN := \
  -not -path ../src/$(SRC_MAIN)\
  -not -path ../src/$(DEBUG_MAIN)\
  -not -path ../src/$(TEST_MAIN)
SRC_FILES := \
  $(shell find $(SRC_DIRS) $(EXCLUDE_MAIN) -name "*.cpp" 2>/dev/null)
OBJ_FILES := $(call path_to_object,$(SRC_FILES))

TEST_SRC_DIRS := ../test \
  $(addsuffix /test/, $(addprefix ../, $(LOCAL_LIBS)))
TEST_SRC_FILES := \
  $(shell find $(TEST_SRC_DIRS) -name "*.cpp" 2>/dev/null)
TEST_OBJ_FILES := $(call path_to_object,$(TEST_SRC_FILES))

.PHONY: all
all: release

## RELEASE ###########################################################
RELEASE_DIR := Release
RELEASE_CXXFLAGS := -O3
RELEASE_LDFLAGS :=
RELEASE_LDLIBS :=

.PHONY: release
release:
	$(call make_target_in_subdir,release_build,$(RELEASE_DIR))

.PHONY: release_build
release_build: MAIN_SRC := ../src/$(RELEASE_MAIN)
release_build: CXXFLAGS += $(RELEASE_CXXFLAGS)
release_build: LDFLAGS += $(RELEASE_LDFLAGS)
release_build: LDLIBS += $(RELEASE_LDLIBS)
release_build: $(TARGET_FILENAME)

## DEBUG #############################################################
DEBUG_DIR := Debug
DEBUG_CXXFLAGS := -g
DEBUG_LDFLAGS := -g
DEBUG_LDLIBS :=

.PHONY: debug
debug:
	$(call make_target_in_subdir,debug_build,$(DEBUG_DIR))

.PHONY: debug_build
debug_build: MAIN_SRC := ../src/$(DEBUG_MAIN)
debug_build: CXXFLAGS += $(DEBUG_CXXFLAGS)
debug_build: LDFLAGS += $(DEBUG_LDFLAGS)
debug_build: LDLIBS += $(DEBUG_LDLIBS)
debug_build: $(TARGET_FILENAME)

## GTEST #############################################################

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
GMOCK_CPPFLAGS := $(CXX_FLAGS) -isystem $(GTEST_DIR)/include -isystem \
	          $(GMOCK_DIR)/include

# Flags passed to the C++ compiler.
GMOCK_CXXFLAGS := -g -Wall -Wextra -pthread

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS := $(GTEST_DIR)/include/gtest/*.h \
                 $(GTEST_DIR)/include/gtest/internal/*.h

# All Google Mock headers. Note that all Google Test headers are
# included here too, as they are #included by Google Mock headers.
# Usually you shouldn't change this definition.
GMOCK_HEADERS = $(GMOCK_DIR)/include/gmock/*.h \
                $(GMOCK_DIR)/include/gmock/internal/*.h \
                $(GTEST_HEADERS)

# Builds gmock.a and gmock_main.a.  These libraries contain both
# Google Mock and Google Test.  A test should link with either gmock.a
# or gmock_main.a, depending on whether it defines its own main()
# function.  It's fine if your test only uses features from Google
# Test (and not Google Mock).

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
GMOCK_SRCS_ = $(GMOCK_DIR)/src/*.cc $(GMOCK_HEADERS)

# For simplicity and to avoid depending on implementation details of
# Google Mock and Google Test, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Mock and
# Google Test compile fast and for ordinary users their source rarely
# changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(GMOCK_CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) \
	$(GMOCK_CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gmock-all.o : $(GMOCK_SRCS_)
	$(CXX) $(GMOCK_CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) \
	 $(GMOCK_CXXFLAGS) -c $(GMOCK_DIR)/src/gmock-all.cc

gmock_main.o : $(GMOCK_SRCS_)
	$(CXX) $(GMOCK_CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) \
	$(GMOCK_CXXFLAGS) -c $(GMOCK_DIR)/src/gmock_main.cc

gmock.a : gmock-all.o gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gmock_main.a : gmock-all.o gtest-all.o gmock_main.o
	$(AR) $(ARFLAGS) $@ $^

## TEST ##############################################################
TEST_DIR := Test
TEST_CXXFLAGS := -g $(GMOCK_CPPFLAGS) -pthread
TEST_LDFLAGS := -g -I$(GTEST_DIR) -I$(GMOCK_DIR)
TEST_LDLIBS := -lpthread

.PHONY: test
test:
	$(call make_target_in_subdir,test_build,$(TEST_DIR))

.PHONY: test_build
test_build: MAIN_SRC := ../src/$(TEST_MAIN)
test_build: CXXFLAGS += $(TEST_CXXFLAGS)
test_build: LDFLAGS += $(TEST_LDFLAGS)
test_build: LDLIBS += $(TEST_LDLIBS)
test_build: test.out
	@- ./test.out > test_summary.txt
	./test.out

test.out: main.o $(OBJ_FILES) $(TEST_OBJ_FILES) gmock_main.a
	$(CXX) main.o $(OBJ_FILES) $(TEST_OBJ_FILES) gmock_main.a \
	-o $@ $(LDFLAGS) $(LDLIBS)


## DOC ###############################################################
DOC_DIR := Documentation

.PHONY: doc
doc: Doxyfile
	(cat Doxyfile ; echo "OUTPUT_DIRECTORY = $(DOC_DIR)") \
	| doxygen -

Doxyfile:
	doxygen -g

## CLEAN #############################################################
.PHONY: clean
clean:
	rm -rf $(RELEASE_DIR) $(DEBUG_DIR) $(TEST_DIR) $(DOC_DIR)

######################################################################
## DEPENDENCIES ######################################################
$(TARGET_FILENAME): main.o $(OBJ_FILES)
	$(CXX) main.o $(OBJ_FILES) -o $@ $(LDFLAGS) $(LDLIBS)

main.o: $(MAIN_SRC)
	$(CXX) $(CXXFLAGS) -c -MMD $(MAIN_SRC) -o $@

%.o: $(call object_to_path,$@)
	$(CXX) $(CXXFLAGS) -c -MMD $(call object_to_path,$@) -o $@

.PHONY: *.d
-include *.d
