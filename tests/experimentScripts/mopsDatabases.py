import MySQLdb as db


OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"

DB_USER="jmyers"
DB_PASS="jmyers"
DB_HOST="localhost"


def getCursor():
    conn = db.connect(user=DB_USER, passwd=DB_PASS, host=DB_HOST)
    curs = conn.cursor()
    return curs
