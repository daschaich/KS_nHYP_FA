CC = mpiicc
COPT = -O3 -std=c99 -Wall -Wno-strict-aliasing

SCIDAC = /g/g10/schaich2/work/lqcd/install
QMPDIR = $(SCIDAC)/qmp
QLADIR = $(SCIDAC)/qla
QIODIR = $(SCIDAC)/qio
QDPDIR = $(SCIDAC)/qdp
QOPDIR = $(SCIDAC)/qopqdp
QCDDIR = $(SCIDAC)/qcdlib

IQOP = -I$(QOPDIR)/include
LQOP = -L$(QOPDIR)/lib -lqopqdp

IQDP = -I$(QDPDIR)/include
LQDP = -L$(QDPDIR)/lib -lqdp_f -lqdp_d -lqdp_df -lqdp_f3 -lqdp_d3 -lqdp_df3 -lqdp_int -lqdp_common

IQLA = -I$(QLADIR)/include
LQLA = -L$(QLADIR)/lib -lqla_f -lqla_d -lqla_q -lqla_df -lqla_dq -lqla_f3 -lqla_d3 -lqla_q3 -lqla_df3 -lqla_dq3 -lqla_int -lqla_random -lqla_cmath

IQIO = -I$(QIODIR)/include
LQIO = -L$(QIODIR)/lib -lqio -llime

IQMP = -I$(QMPDIR)/include
LQMP = -L$(QMPDIR)/lib -lqmp

IQCD = -I$(QCDDIR)
LQCD = $(QCDDIR)/libqcd.a

ILAPACK = -DUSE_UNDERSCORES
LLAPACK = -llapack -lblas

CFLAGS = $(COPT) $(IQCD) $(IQOP) $(IQDP) $(IQLA) $(IQIO) $(IQMP) $(ILAPACK) -DQDP_PrecisionInt=2
LDLIBS = $(LQCD) $(LQOP) $(LQDP) $(LQLA) $(LQIO) $(LQMP) $(LLAPACK) -lfftw3 -lm

PROGS = KSdecomp

all: $(PROGS)

KSdecomp: KSdecomp.o

clean:
	-rm -rf $(PROGS) *.o
