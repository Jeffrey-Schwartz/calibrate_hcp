# You will likely have to change the following two lines
module_LTLIBRARIES = calibrate_hcp.la
calibrate_hcp_la_SOURCES = calibrate_hcp.c

# The rest is quite generic unless your module uses extra libraries
ACLOCAL_AMFLAGS = -I m4
moduledir = @GWYDDION_MODULE_DIR@
AM_CPPFLAGS = -I$(top_srcdir) -DG_LOG_DOMAIN=\"Module\" @GWYDDION_CFLAGS@
AM_CFLAGS = @WARNING_CFLAGS@ @HOST_CFLAGS@
AM_LDFLAGS = -avoid-version -module @HOST_LDFLAGS@ @GWYDDION_LIBS@
