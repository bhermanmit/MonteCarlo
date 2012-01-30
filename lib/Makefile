# Makefile for a serial random number generator.
program = libpdf.a

#===============================================================================
# User Options
#===============================================================================

COMPILER = gnu
DEBUG    = no

#===============================================================================
# Compiler Options
#===============================================================================

# Intel Fortran compiler options

ifeq ($(COMPILER), intel)
  F90 = ifort
  F90FLAGS += -fpp -warn
endif

# GNU Fortran compiler options

ifeq ($(COMPILER), gnu)
  F90 = gfortran
  F90FLAGS += -cpp -Wall
endif

# Set compiler flags for debugging

ifeq ($(DEBUG),yes)
  F90FLAGS += -g 
  LDFLAGS  += -g
  ifeq ($(COMPILER), intel)
    F90FLAGS += -traceback -ftrapuv -fp-stack-check -check all
  endif
  ifeq ($(COMPILER), gfortran)
    F90FLAGS += -pedantic -std=f2008 -fbacktrace -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
  endif
endif

RANLIB = ranlib

#===============================================================================
# Targets
#===============================================================================

OBJS = random_number_generator.o \
       pdfs.o 
      
all : $(program)

$(program): $(OBJS)
	ar ruv $(program) $(OBJS)
	$(RANLIB) $(program)

clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: clean neat

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<
