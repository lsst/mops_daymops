/* I think this is right... */

CREATE TABLE fullerDiaSource
(
	diaSourceId BIGINT NOT NULL PRIMARY KEY,
	opSimId BIGINT NULL,
	ssmId BIGINT NULL,

	ra DOUBLE NOT NULL,
	decl DOUBLE NOT NULL,
	taiMidPoint DOUBLE NOT NULL,
	mag DOUBLE NOT NULL,
	snr DOUBLE NOT NULL
) ;

LOAD DATA INFILE '/workspace1/jmyers/nightlyDiasAstromErr/dias_pt1_nodeep.short.astromErr.plusIds' 
 INTO TABLE mops_diaSource
 FIELDS TERMINATED BY ' '
 LINES TERMINATED BY '\n' 
  (diaSourceId, opSimId, ssmId, ra, decl, taiMidPoint,
   mag,  snr) 
;
