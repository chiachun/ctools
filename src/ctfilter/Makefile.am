## Process this file with automake to produce Makefile.in

# Program name
bin_PROGRAMS = ctfilter

# Sources
sources = ctfilter.cpp

# Libtool library (only needed for python bindings)
noinst_LTLIBRARIES = libctfilter.la
libctfilter_la_SOURCES = $(sources)

# Program definition
ctfilter_SOURCES = main.cpp
ctfilter_LDFLAGS = $(LIB_LINK_OPTION) @LDFLAGS@
ctfilter_CPPFLAGS = @CPPFLAGS@
ctfilter_LDADD = libctfilter.la

# Parameter files
pfilesdir = $(prefix)/syspfiles
pfiles_DATA = ctfilter.par
