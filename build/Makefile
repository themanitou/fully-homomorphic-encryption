##################################################################
#
#   Fully-Homomorphic Cryptography library,
#   based on Gentry-Halevi ideal lattice scheme.
#
#   Author: Quan Nguyen (https://github.com/themanitou)
#
#   This library is open-source software distributed under the
#   terms of the GNU Lesser General Public License (LGPL) version
#   2.1 or later.  See the file doc/copying.txt for complete
#   details on the licensing of this library.
#
###################################################################

ifndef $(CC)
CC=g++
endif

CPPFLAGS += -O2 -m64 -DFHE_TEST -DFHE_MPI_LOG -std=c++17 `pkg-config --cflags mpi-cxx`
LIBS     += -lntl `pkg-config --libs mpi-cxx`

-include src/tools/subdir.mk
-include src/subdir.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

all: ../bin/FullHomCrypt.out

../bin/FullHomCrypt.out: $(OBJS) $(USER_OBJS)
	@echo 'Building file: $@'
	mkdir -p ../bin
	@echo 'Invoking: Cross GCC Linker'
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o "$@" $(OBJS) $(LIBS)
	@echo 'Finished building: $@'
	@echo ' '

clean:
	-$(RM) $(OBJS) $(C_DEPS) ../bin/FullHomCrypt.out
	-@echo ' '
