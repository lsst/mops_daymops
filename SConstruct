# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re, os
import lsst.SConsUtils as scons


env = scons.makeEnv("mops_daymops", 
                   r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/mops/daymops/trunk/SConstruct $",
                   [["boost", "boost/version.hpp", "boost_system:C++"],
                    ["boost", "boost/version.hpp", "boost_filesystem:C++"],
                    ["boost", "boost/regex.hpp", "boost_regex:C++"],
                    ["boost", "boost/filesystem.hpp", "boost_system:C++"],
                    ["boost", "boost/serialization/base_object.hpp", "boost_serialization:C++"],
                    ["boost", "boost/test/included/unit_test.hpp"],
                    ["pex_exceptions", "lsst/pex/exceptions/Exception.h lsst/pex/exceptions/Runtime.h", "pex_exceptions:C++"],
                    ["gsl", "gsl/gsl_fit.h", "gslcblas gsl:C++"],
                    ["pal", "pal.h palmac.h", "libpal:C++"],
                    ["eigen", "Eigen/Core.h"],
                    #["utils", "lsst/tr1/unordered_map.h", "utils:C++"],
                    #["daf_base", "lsst/daf/base/Citizen.h lsst/daf/base/Persistable.h", "daf_base:C++"],
                    #["pex_policy", "lsst/pex/policy.h", "pex_policy:C++"],
                    #["daf_persistence", "lsst/daf/persistence.h", "daf_persistence:C++"],
                    ["python", "Python.h"]
                    ])

env.libs["mops_daymops"] += env.getlibs("boost pex_exceptions gsl python pal")

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

# the "install" target
#
Alias("install", [ env.Install(env['prefix'], "include"),
                  env.InstallEups(env['prefix'] + "/ups",
                                  glob.glob("ups/*.table"))])

scons.CleanTree(r"*~ core *.os *.o *.a *.so")

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
            'lib/SConscript',
            'doc/SConscript'])

SConscript('doc/SConscript')

 
