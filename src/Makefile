# jmyers march 27 2012 - updated ljones Dec 19, 2012
# makefile to get things building on Gordon or elsewhere (without scons)

# -DNOPEX: build our own Exceptions which derive from
# std::exception, not LSST's pex_exception. This way we avoid needing
# any LSST packages at all!

# choose your appropriate gcc compiler
GCC = c++

# if you compiled PAL as a shared library, use this OPT line
OPT=-O3 -DNOPEX -fPIC -shared
# if you compiled PAL as a static library, use this OPT line
#OPT=-O3 -DNOPEX 
#if not using EUPS, then you'll probably uncomment these (otherwise just use eups to set them up first)
#THIRD=~/thirdPartyPrecompiled/
#EIGEN_DIR=${THIRD}eigen/
#BOOST_DIR=${THIRD}boost/
#GSL_DIR=/opt/gsl/gnu/
PAL_DIR=${PAL_DIR}

EXTINCLUDES=-I ${EIGEN_DIR}/include/ -I ${BOOST_DIR}/include/ -I ${PAL_DIR}/include/ -I ${GSL_DIR}/include
EXTLIBDIRS=-L ${BOOST_DIR}/lib -L ${PAL_DIR}/lib -L ${GSL_DIR}/lib
EXTLIBS=-lgsl -lgslcblas  -lpal

# we are being really suboptimal here and just recompiling EVERYTHING if any header changes. Fix this someday if it gets too painful
BASEINC=-I ../include/
MOPSINC=../include/lsst/mops/
MOPSHEADERS= ${MOPSINC}*.h ${MOPSINC}daymops/findTracklets/*.h ${MOPSINC}daymops/collapseTrackletsAndPostfilters/*.h 

all:	../bin/findTracklets ../bin/collapseTracklets ../bin/purifyTracklets ../bin/removeSubsets ../bin/linkTracklets 

# ../bin/findTrackletsOMP   ../bin/collapseTrackletsOMP  ../bin/purifyTrackletsOMP   ../bin/removeSubsetsOMP   ../bin/linkTrackletsOMP 

clean: 
	rm -f *.o */*.o ../bin/findTracklets ../bin/collapseTracklets ../bin/purifyTracklets ../bin/removeSubsets ../bin/findTrackletsOMP ../bin/linkTracklets   ../bin/collapseTrackletsOMP  ../bin/purifyTrackletsOMP   ../bin/removeSubsetsOMP   ../bin/linkTrackletsOMP 

MopsDetection.o: MopsDetection.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c  MopsDetection.cc ${EXTINCLUDES} ${BASEINC}

PointAndValue.o: PointAndValue.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c  PointAndValue.cc ${EXTINCLUDES} ${BASEINC}

common.o: common.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c  common.cc ${EXTINCLUDES} ${BASEINC}

Tracklet.o: Tracklet.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c  Tracklet.cc ${EXTINCLUDES} ${BASEINC}

Track.o: Track.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c Track.cc ${EXTINCLUDES} ${BASEINC}

TrackletVector.o: TrackletVector.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c TrackletVector.cc ${EXTINCLUDES} ${BASEINC}

fileUtils.o: fileUtils.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c fileUtils.cc ${EXTINCLUDES} ${BASEINC}

rmsLineFit.o: rmsLineFit.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c rmsLineFit.cc ${EXTINCLUDES} ${BASEINC}

TrackletTree.o: linkTracklets/TrackletTree.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c linkTracklets/TrackletTree.cc ${EXTINCLUDES} ${BASEINC}

TrackletTreeNode.o: linkTracklets/TrackletTreeNode.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c linkTracklets/TrackletTreeNode.cc ${EXTINCLUDES} ${BASEINC}

TrackSet.o: TrackSet.cc ${MOPSHEADERS}
	${GCC} ${OPT} -c TrackSet.cc ${EXTINCLUDES} ${BASEINC}

../bin/findTracklets: findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o \
findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc ${EXTLIBS} -o ../bin/findTracklets

../bin/findTrackletsOMP: findTracklets/findTrackletsOMP.cc findTracklets/findTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o \
-lgomp -fopenmp \
findTracklets/findTrackletsOMP.cc findTracklets/findTrackletsMain.cc ${EXTLIBS} -o ../bin/findTrackletsOMP

../bin/collapseTracklets: collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc ${EXTLIBS} -o ../bin/collapseTracklets

../bin/collapseTrackletsOMP: collapseTrackletsAndPostfilters/collapseTrackletsOMP.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
collapseTrackletsAndPostfilters/collapseTrackletsOMP.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc ${EXTLIBS} -o ../bin/collapseTrackletsOMP

../bin/purifyTracklets: collapseTrackletsAndPostfilters/purifyTracklets.cc  MopsDetection.o common.o Tracklet.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
collapseTrackletsAndPostfilters/purifyTracklets.cc ${EXTLIBS} -o ../bin/purifyTracklets

../bin/purifyTrackletsOMP: collapseTrackletsAndPostfilters/purifyTrackletsOMP.cc  MopsDetection.o common.o Tracklet.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
collapseTrackletsAndPostfilters/purifyTrackletsOMP.cc ${EXTLIBS} -o ../bin/purifyTrackletsOMP

../bin/removeSubsets: removeSubsets.cc removeSubsetsMain.cc MopsDetection.o common.o Tracklet.o fileUtils.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o common.o Tracklet.o fileUtils.o \
removeSubsets.cc removeSubsetsMain.cc ${EXTLIBS} -o ../bin/removeSubsets

../bin/removeSubsetsOMP: removeSubsetsOMP.cc removeSubsetsMainOMP.cc MopsDetection.o common.o Tracklet.o fileUtils.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
MopsDetection.o common.o Tracklet.o fileUtils.o \
-fopenmp -lgomp \
removeSubsetsOMP.cc removeSubsetsMainOMP.cc ${EXTLIBS} -o ../bin/removeSubsetsOMP

../bin/linkTracklets: linkTracklets/linkTracklets.cc linkTracklets/linkTrackletsMain.cc TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o \
linkTracklets/linkTracklets.cc linkTracklets/linkTrackletsMain.cc ${EXTLIBS} -o ../bin/linkTracklets

../bin/linkTrackletsOMP: linkTracklets/linkTrackletsOMP.cc linkTracklets/linkTrackletsMain.cc TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	${GCC} ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} \
TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
linkTracklets/linkTrackletsOMP.cc linkTracklets/linkTrackletsMain.cc ${EXTLIBS} -o ../bin/linkTracklets
