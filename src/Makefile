# jmyers march 27 2012
# makefile to get things building on Gordon or

# -DPANSTARRS: build our own Exceptions which derive from
# std::exception, not LSST's pex_exception. This way we avoid needing
# any LSST packages at all!
OPT=-O3 -DNOPEX
#if not using EUPS, then you'll probably uncomment these
THIRD=~/thirdPartyPrecompiled/
SLALIB_DIR=${THIRD}slalib/
EIGEN_DIR=${THIRD}eigen/3.0.0
BOOST_DIR=${THIRD}boost/1.37.0+2
GSL_DIR=/opt/gsl/gnu/

EXTINCLUDES=-I ${EIGEN_DIR}/include/ -I ${BOOST_DIR}/include/ -I ${SLALIB_DIR}/include/ -I ${GSL_DIR}/include
EXTLIBDIRS=-L ${BOOST_DIR}/lib -L ${SLALIB_DIR}/lib -L ${GSL_DIR}/lib
EXTLIBS=-lgsl -lgslcblas ${SLALIB_DIR}/src/*.o

#we are being really suboptimal here and just recompiling EVERYTHING if any header changes. Fix this someday if it gets too painful
BASEINC=-I ../include/
MOPSINC=../include/lsst/mops/
MOPSHEADERS= ${MOPSINC}*.h ${MOPSINC}daymops/findTracklets/*.h ${MOPSINC}daymops/collapseTrackletsAndPostfilters/*.h 

all:	../bin/findTracklets ../bin/collapseTracklets ../bin/purifyTracklets ../bin/removeSubsets ../bin/linkTracklets  ../bin/findTrackletsOMP   ../bin/collapseTrackletsOMP  ../bin/purifyTrackletsOMP   ../bin/removeSubsetsOMP   ../bin/linkTrackletsOMP 

clean: 
	rm -f *.o */*.o ../bin/findTracklets ../bin/collapseTracklets ../bin/purifyTracklets ../bin/removeSubsets ../bin/findTrackletsOMP ../bin/linkTracklets   ../bin/collapseTrackletsOMP  ../bin/purifyTrackletsOMP   ../bin/removeSubsetsOMP   ../bin/linkTrackletsOMP 

MopsDetection.o: MopsDetection.cc ${MOPSHEADERS}
	g++ ${OPT} -c  MopsDetection.cc ${EXTINCLUDES} ${BASEINC}

PointAndValue.o: PointAndValue.cc ${MOPSHEADERS}
	g++ ${OPT} -c  PointAndValue.cc ${EXTINCLUDES} ${BASEINC}

common.o: common.cc ${MOPSHEADERS}
	g++ ${OPT} -c  common.cc ${EXTINCLUDES} ${BASEINC}

Tracklet.o: Tracklet.cc ${MOPSHEADERS}
	g++ ${OPT} -c  Tracklet.cc ${EXTINCLUDES} ${BASEINC}

Track.o: Track.cc ${MOPSHEADERS}
	g++ ${OPT} -c Track.cc ${EXTINCLUDES} ${BASEINC}

TrackletVector.o: TrackletVector.cc ${MOPSHEADERS}
	g++ ${OPT} -c TrackletVector.cc ${EXTINCLUDES} ${BASEINC}

fileUtils.o: fileUtils.cc ${MOPSHEADERS}
	g++ ${OPT} -c fileUtils.cc ${EXTINCLUDES} ${BASEINC}

rmsLineFit.o: rmsLineFit.cc ${MOPSHEADERS}
	g++ ${OPT} -c rmsLineFit.cc ${EXTINCLUDES} ${BASEINC}

TrackletTree.o: linkTracklets/TrackletTree.cc ${MOPSHEADERS}
	g++ ${OPT} -c linkTracklets/TrackletTree.cc ${EXTINCLUDES} ${BASEINC}

TrackletTreeNode.o: linkTracklets/TrackletTreeNode.cc ${MOPSHEADERS}
	g++ ${OPT} -c linkTracklets/TrackletTreeNode.cc ${EXTINCLUDES} ${BASEINC}

TrackSet.o: TrackSet.cc ${MOPSHEADERS}
	g++ ${OPT} -c TrackSet.cc ${EXTINCLUDES} ${BASEINC}

../bin/findTracklets: findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o \
findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc -o ../bin/findTracklets

../bin/findTrackletsOMP: findTracklets/findTrackletsOMP.cc findTracklets/findTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o \
-lgomp -fopenmp \
findTracklets/findTrackletsOMP.cc findTracklets/findTrackletsMain.cc -o ../bin/findTrackletsOMP

../bin/collapseTracklets: collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc -o ../bin/collapseTracklets

../bin/collapseTrackletsOMP: collapseTrackletsAndPostfilters/collapseTrackletsOMP.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
collapseTrackletsAndPostfilters/collapseTrackletsOMP.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc -o ../bin/collapseTrackletsOMP

../bin/purifyTracklets: collapseTrackletsAndPostfilters/purifyTracklets.cc  MopsDetection.o common.o Tracklet.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
collapseTrackletsAndPostfilters/purifyTracklets.cc -o ../bin/purifyTracklets

../bin/purifyTrackletsOMP: collapseTrackletsAndPostfilters/purifyTrackletsOMP.cc  MopsDetection.o common.o Tracklet.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
collapseTrackletsAndPostfilters/purifyTrackletsOMP.cc -o ../bin/purifyTrackletsOMP

../bin/removeSubsets: removeSubsets.cc removeSubsetsMain.cc MopsDetection.o common.o Tracklet.o fileUtils.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o common.o Tracklet.o fileUtils.o \
removeSubsets.cc removeSubsetsMain.cc -o ../bin/removeSubsets

../bin/removeSubsetsOMP: removeSubsetsOMP.cc removeSubsetsMainOMP.cc MopsDetection.o common.o Tracklet.o fileUtils.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o common.o Tracklet.o fileUtils.o \
-fopenmp -lgomp \
removeSubsetsOMP.cc removeSubsetsMainOMP.cc -o ../bin/removeSubsetsOMP

../bin/linkTracklets: linkTracklets/linkTracklets.cc linkTracklets/linkTrackletsMain.cc TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o \
linkTracklets/linkTracklets.cc linkTracklets/linkTrackletsMain.cc -o ../bin/linkTracklets

../bin/linkTrackletsOMP: linkTracklets/linkTrackletsOMP.cc linkTracklets/linkTrackletsMain.cc TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
TrackletTree.o TrackletTreeNode.o TrackSet.o Tracklet.o Track.o MopsDetection.o common.o fileUtils.o rmsLineFit.o \
-fopenmp -lgomp \
linkTracklets/linkTrackletsOMP.cc linkTracklets/linkTrackletsMain.cc -o ../bin/linkTracklets