-- New schema for Mops tables (only)
-- To later be merged in to the full LSST schema
-- jmyers
-- started Jul 9 2010


-- "Lightweight" DiaSources are used by MOPS; they only hold the
-- information needed to find tracks and derive orbits and have a much
-- smaller footprint than the conventional LSST DiaSources.
--
-- Template for tables which hold abbreviated diaSource information
-- for ONE NIGHT.  Expect one table per every night of observation.
-- Look up the correct table name using
-- mops_nightNumToNightlyDiaSourceTableName
--
-- By convention, we will always use the same DiaSourceID as the full,
-- heavyweight DiaSource stored by the rest of LSST.

DROP TABLE IF EXISTS _mops_nightlyDiaSources; 

CREATE TABLE _mops_nightlyDiaSources
(
	diaSourceId BIGINT NOT NULL PRIMARY KEY,
	objectId BIGINT NULL,
	movingObjectId BIGINT NULL,

	ra DOUBLE NOT NULL,
	decl DOUBLE NOT NULL,
	taiMidPoint DOUBLE NOT NULL,
	mag DOUBLE NOT NULL,
	filterId TINYINT NOT NULL
) ;



-- maps night number to name of the correct mops_nightlyDiaSources 

DROP TABLE IF EXISTS mops_nightNumToNightlyDiaSourceTableName;

CREATE TABLE mops_nightNumToNightlyDiaSourceTableName
(
	nightNum INTEGER NOT NULL PRIMARY KEY,
	diaSourceTableName VARCHAR(200)
) ;




-- templates for tables which define tracklets from ONE NIGHT.  Expect
-- one table per every night of observation.  Look up the correct
-- table name using mops_nightNumToNightlyTrackletToDIASource


DROP TABLE IF EXISTS _mops_nightlyTrackletToDIASource;

CREATE TABLE _mops_nightlyTrackletToDIASource
(
	trackletId BIGINT NOT NULL,
	diaSourceId BIGINT NOT NULL,
	PRIMARY KEY (trackletId, diaSourceId),
	INDEX idx_mopsTrackletToDIASource_diaSourceId (diaSourceId ASC),
	INDEX idx_mopsTrackletToDIASource_trackletId (trackletId ASC),
	KEY (trackletId)
) ;

DROP TABLE IF EXISTS _mops_nightlyTracklet;

CREATE TABLE _mops_nightlyTracklet
(
	trackletId BIGINT NOT NULL PRIMARY KEY,
	ra0 DOUBLE,
	dec0 DOUBLE,
	raV DOUBLE,
	decV DOUBLE
) ;




-- maps night number to name of the correct
-- mops_nightlyTrackletToDIASource table and mops_Tracklet table for
-- the given night number

DROP TABLE IF EXISTS mops_nightNumToNightlyTrackletTableNames;

CREATE TABLE mops_nightNumToNightlyTrackletTableNames
(
	nightNum INTEGER NOT NULL PRIMARY KEY,
	trackletToDiaSourceTableName VARCHAR(200),
	trackletTableName VARCHAR(200)
) ;






-- templates for tables which define tracks from ONE NIGHT.  Expect one
-- table per every night of observation.  Look up the correct table
-- name using mops_nightNumToNightlyTrackletToDIASource

DROP TABLE IF EXISTS _mops_nightlyTrackToDIASource;

CREATE TABLE _mops_nightlyTrackToDIASource
(
	trackId BIGINT NOT NULL,
	diaSourceId BIGINT NOT NULL,
	PRIMARY KEY (trackId, diaSourceId),
	INDEX idx_mopsTrackToDIASource_diaSourceId (diaSourceId ASC),
	INDEX idx_mopsTrackToDIASource_TrackId (diaSourceId ASC)
) ;

DROP TABLE IF EXISTS _mops_nightlyTrack;

CREATE TABLE _mops_nightlyTrack
(
	trackId BIGINT NOT NULL PRIMARY KEY,
	ra0    DOUBLE,
	dec0   DOUBLE,
	raV    DOUBLE,
	decV   DOUBLE,
	raAcc  DOUBLE,
	decAcc DOUBLE
) ;


-- maps night number to name of the correct mops_nightlyTrackletsToDIASource

DROP TABLE IF EXISTS mops_nightNumToNightlyTrackTableNames;

CREATE TABLE mops_nightNumToNightlyTrackTableNames
(
	nightNum INTEGER NOT NULL PRIMARY KEY,
	trackToDiaSourceTableName VARCHAR(200),
	trackTableName VARCHAR(200)
) ;





-- MovingObject does not need to be partitioned (at least yet) so they are monolithic tables

DROP TABLE IF EXISTS mops_MovingObjectToDiaSource;

CREATE TABLE mops_MovingObjectToDiaSource
(
	movingObjectId BIGINT NOT NULL,
	diaSourceId BIGINT NOT NULL,
	INDEX idx_mopsMovingObjectToDiaSource_movingObjectId (movingObjectId ASC),
	INDEX idx_mopsMovingObjectToDiaSource_diaSourceId (diaSourceId ASC)
) ;

DROP TABLE IF EXISTS MovingObject;

CREATE TABLE MovingObject
(
	movingObjectId BIGINT NOT NULL,
	movingObjectVersion INT NOT NULL DEFAULT '1',
	procHistoryId INTEGER NULL,
	taxonomicTypeId SMALLINT NULL,
	ssmObjectName VARCHAR(32) NULL,
	q DOUBLE NOT NULL,
	e DOUBLE NOT NULL,
	i DOUBLE NOT NULL,
	node DOUBLE NOT NULL,
	meanAnom DOUBLE NULL,
	argPeri DOUBLE NOT NULL,
	distPeri DOUBLE NOT NULL,
	timePeri DOUBLE NOT NULL,
	epoch DOUBLE NOT NULL,
	h_v DOUBLE NOT NULL,
	g DOUBLE NULL DEFAULT 0.15,
	rotationPeriod DOUBLE NULL,
	rotationEpoch DOUBLE NULL,
	albedo DOUBLE NULL,
	poleLat DOUBLE NULL,
	poleLon DOUBLE NULL,
	d3 DOUBLE NULL,
	d4 DOUBLE NULL,
	orbFitResidual DOUBLE NOT NULL,
	orbFitChi2 DOUBLE NULL,
	classification CHAR(1) NULL,
	ssmId BIGINT NULL,
	mopsStatus CHAR(1) NULL,
	stablePass CHAR(1) NULL,
	timeCreated TIMESTAMP NULL,
	uMag DOUBLE NULL,
	uMagErr FLOAT(0) NULL,
	uAmplitude FLOAT(0) NULL,
	uPeriod FLOAT(0) NULL,
	gMag DOUBLE NULL,
	gMagErr FLOAT(0) NULL,
	gAmplitude FLOAT(0) NULL,
	gPeriod FLOAT(0) NULL,
	rMag DOUBLE NULL,
	rMagErr FLOAT(0) NULL,
	rAmplitude FLOAT(0) NULL,
	rPeriod FLOAT(0) NULL,
	iMag DOUBLE NULL,
	iMagErr FLOAT(0) NULL,
	iAmplitude FLOAT(0) NULL,
	iPeriod FLOAT(0) NULL,
	zMag DOUBLE NULL,
	zMagErr FLOAT(0) NULL,
	zAmplitude FLOAT(0) NULL,
	zPeriod FLOAT(0) NULL,
	yMag DOUBLE NULL,
	yMagErr FLOAT(0) NULL,
	yAmplitude FLOAT(0) NULL,
	yPeriod FLOAT(0) NULL,
	flag INTEGER NULL,
	src01 DOUBLE NULL,
	src02 DOUBLE NULL,
	src03 DOUBLE NULL,
	src04 DOUBLE NULL,
	src05 DOUBLE NULL,
	src06 DOUBLE NULL,
	src07 DOUBLE NULL,
	src08 DOUBLE NULL,
	src09 DOUBLE NULL,
	src10 DOUBLE NULL,
	src11 DOUBLE NULL,
	src12 DOUBLE NULL,
	src13 DOUBLE NULL,
	src14 DOUBLE NULL,
	src15 DOUBLE NULL,
	src16 DOUBLE NULL,
	src17 DOUBLE NULL,
	src18 DOUBLE NULL,
	src19 DOUBLE NULL,
	src20 DOUBLE NULL,
	src21 DOUBLE NULL,
	convCode VARCHAR(8) NULL,
	o_minus_c DOUBLE NULL,
	moid1 DOUBLE NULL,
	moidLong1 DOUBLE NULL,
	moid2 DOUBLE NULL,
	moidLong2 DOUBLE NULL,
	arcLengthDays DOUBLE NULL,
	PRIMARY KEY (movingObjectId, movingObjectVersion),
	KEY (procHistoryId),
	INDEX idx_MovingObject_taxonomicTypeId (taxonomicTypeId ASC),
	INDEX idx_MovingObject_ssmId (ssmId ASC),
	INDEX idx_MovingObject_ssmObjectName (ssmObjectName ASC),
	INDEX idx_MovingObject_status (mopsStatus ASC)
) ;
