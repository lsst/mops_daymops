# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re, os
import lsst.SConsUtils as scons

env = scons.makeEnv("daymops", 
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/mops/daymops/trunk/SConstruct $",
                    [["boost", "boost/test/included/unit_test.hpp"],
                     ["pex_exceptions", "lsst/pex/exceptions/Exception.h lsst/pex/exceptions/Runtime.h", "pex_exceptions:C++"],
                     ["gsl", "gsl/gsl_fit.h", "gslcblas gsl:C++"]
                     ])



env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

# the "install" target
#
Alias("install", [env.Install(env['prefix'], "python"),
                  env.InstallEups(env['prefix'] + "/ups",
                                  glob.glob("ups/*.table"))])

scons.CleanTree(r"*~ core *.os *.o *.a")

#
# Build TAGS files
#
files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
LSST-MOPS package
""")





SConscript(['src/SConscript',
            'lib/SConscript'])
