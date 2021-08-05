C_SRCS += \
../src/FullHomCrypt.cpp

OBJS += \
./src/FullHomCrypt.o

C_DEPS += \
./src/FullHomCrypt.d


src/%.o: ../src/%.cpp
	@echo 'Building file: $@'
	@echo 'Invoking: Cross GCC Compiler'
	$(CC) $(CPPFLAGS) $(LDFLAGS) -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $@'
	@echo ' '
