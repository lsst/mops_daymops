# jmyers march 27 2012
# makefile to get things building on Gordon or

# -DPANSTARRS: build our own Exceptions which derive from
# std::exception, not LSST's pex_exception. This way we avoid needing
# any LSST packages at all!
OPT=-O3 -DPANSTARRS
THIRD=~/thirdPartyPrecompiled/

EXTINCLUDES=-I ${THIRD}eigen/3.0.0/include/ -I ${THIRD}1.37.0+2/include/ -I ~/build/include/ -I ${THIRD}slalib/include/
EXTLIBDIRS=-L ${THIRD}boost/1.37.0+2/lib -L ~/build/lib -L ${THIRD}slalib/lib -L /opt/atlas/gnu/lib/
EXTLIBS=-lgsl -lcblas -latlas ${THIRD}slalib/src/*.o

#we are being really suboptimal here and just recompiling EVERYTHING if any header changes. Fix this someday if it gets too painful
BASEINC=-I ../include/
MOPSINC=../include/lsst/mops/
MOPSHEADERS= ${MOPSINC}*.h ${MOPSINC}daymops/findTracklets/*.h ${MOPSINC}daymops/collapseTrackletsAndPostfilters/*.h 

all:	../bin/findTracklets ../bin/collapseTracklets


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

../bin/findTracklets: findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o \
findTracklets/findTracklets.cc findTracklets/findTrackletsMain.cc -o ../bin/findTracklets

../bin/collapseTracklets: collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o ${MOPSHEADERS}
	g++ ${OPT} ${BASEINC} ${EXTINCLUDES} ${EXTLIBDIRS} ${EXTLIBS} \
MopsDetection.o PointAndValue.o common.o Tracklet.o Track.o TrackletVector.o fileUtils.o rmsLineFit.o \
collapseTrackletsAndPostfilters/collapseTracklets.cc collapseTrackletsAndPostfilters/collapseTrackletsMain.cc -o ../bin/collapseTracklets

