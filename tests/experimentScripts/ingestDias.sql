CREATE TABLE mops_diaSource
(
	diaSourceId BIGINT NOT NULL PRIMARY KEY,
	objectId BIGINT NULL,
	groundTruthMovingObjectId BIGINT NULL,
	movingObjectId BIGINT NULL,

	ra DOUBLE NOT NULL,
	decl DOUBLE NOT NULL,
	taiMidPoint DOUBLE NOT NULL,
	mag DOUBLE NOT NULL,
	filterId TINYINT NOT NULL
) ;

LOAD DATA INFILE '/workspace1/jmyers/nightlyDiasAstromErr/dias_pt1_nodeep.short.astromErr.plusIds' 
 INTO TABLE mops_diaSource
 FIELDS TERMINATED BY ' '
 LINES TERMINATED BY '\n' 
  (diaSourceId, @ignored1, groundTruthMovingObjectId, ra, decl, taiMidPoint,
   mag,  @ignored2) 
;
