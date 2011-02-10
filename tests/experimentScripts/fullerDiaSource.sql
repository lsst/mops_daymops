-- jmyers feb 8 2011
-- 
-- A few of our new scripts (orbit_server input generator, some
-- analyses) require diaSource info to have SNR information. 
--
-- This is the schema they expect.
-- 
-- +-------------+------------+------+-----+---------+-------+
-- | Field       | Type       | Null | Key | Default | Extra |
-- +-------------+------------+------+-----+---------+-------+
-- | diaSourceId | bigint(20) | YES  | MUL | NULL    |       | 
-- | opSimId     | bigint(20) | YES  | MUL | NULL    |       | 
-- | ssmId       | bigint(20) | YES  |     | NULL    |       | 
-- | ra          | double     | YES  |     | NULL    |       | 
-- | decl        | double     | YES  |     | NULL    |       | 
-- | taiMidPoint | double     | YES  |     | NULL    |       | 
-- | mag         | double     | YES  |     | NULL    |       | 
-- | snr         | double     | YES  |     | NULL    |       | 
-- +-------------+------------+------+-----+---------+-------+


CREATE TABLE IF NOT EXISTS fullerDiaSource 
       ( 
       	 diaSourceId BIGINT,
	 opSimId BIGINT,
	 ssmId BIGINT,
	 ra DOUBLE,
	 decl DOUBLE,
	 taiMidPoint DOUBLE,
	 mag DOUBLE,
	 snr DOUBLE
	 );	     