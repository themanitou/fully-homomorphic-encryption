C_SRCS += \
../src/tools/fullhom.cpp \
../src/tools/multiplication.cpp \
../src/tools/tools.cpp \
../src/tools/fhe_mpi.cpp

OBJS += \
./src/tools/fullhom.o \
./src/tools/multiplication.o \
./src/tools/tools.o \
./src/tools/fhe_mpi.o

C_DEPS += \
./src/tools/fullhom.d \
./src/tools/multiplication.d \
./src/tools/tools.d \
./src/tools/fhe_mpi.d


src/tools/%.o: ../src/tools/%.cpp
	@echo 'Building file: $@'
	@echo 'Invoking: Cross GCC Compiler'
	$(CC) $(CPPFLAGS) $(LDFLAGS) -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $@'
	@echo ' '
