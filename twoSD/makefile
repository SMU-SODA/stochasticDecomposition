# Make file for Tow-stage Stochastic Decomposition 

#------------------------------------------------------------
# TODO: Replace with appropriate folder names and system descriptions
# Location of CPLEX installations and system descriptions
CPLEXDIR = /opt/ilog/cplex128/cplex
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXINCDIR   = $(CPLEXDIR)/include
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

SUBDIRS := ./src

# These utilities are from the spAlgorithms repository. You may have to
# replace the path to this repository. 
SODA_UTILS := ../../spAlgorithms/solverUtilities
SMPS_SRC := ../../spAlgorithms/smpsReader

#------------------------------------------------------------
# Compiler, Linker Selections and Definitions (CC is for c)
OBJ_SRCS := 
ASM_SRCS := 
C_SRCS := $(wildcard $(SUBDIRS)/*.c) $(SODA_UTILS)/solver.c $(SODA_UTILS)/utility.c $(SMPS_SRC)/prob.c $(SMPS_SRC)/rvgen.c $(SMPS_SRC)/smps.c
O_SRCS := 
S_UPPER_SRCS := 
EXECUTABLES := 
OBJS := $(C_SRCS:.c=.o)
C_DEPS := $(C_SRCS:.c=.d)

CC  := gcc -O3
LIBS := -lilocplex -lcplex -lpthread -lm -ldl
INCLUDES := -I$(CPLEXINCDIR) -I$(SUBDIRS) -I$(SODA_UTILS) -I$(SMPS_SRC)

#------------------------------------------------------------
# All Target
all: twoSD

# Tool invocations
twoSD: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross GCC Linker'
	gcc -L$(CPLEXLIBDIR) -o "twoSD" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o: %.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc $(INCLUDES) -O0 -g3 -Wall -c -fmessage-length=0 -fPIC -fexceptions -m64 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" -c "$<"
	@echo 'Finished building: $<'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) twoSD
	-@echo ' '

