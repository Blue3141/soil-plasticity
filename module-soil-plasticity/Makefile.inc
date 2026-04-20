## Makefile.inc for module-soil-plasticity
## Copy this directory into mbdyn/modules/ and configure with:
##   ./configure --with-module=soil-plasticity

module_LTLIBRARIES = module-soil-plasticity.la
module_soil_plasticity_la_SOURCES = module-soil-plasticity.cc
module_soil_plasticity_la_CPPFLAGS = $(AM_CPPFLAGS)
module_soil_plasticity_la_LDFLAGS = -module -avoid-version
