import MySQLdb as db
import MySQLdb.cursors as cursors

OPSIM_DB="opsim_3_61"
OPSIM_TABLE="output_opsim3_61"

DIAS_DB="mops_noDeepAstromError"
DIAS_TABLE="fullerDiaSource"
#DIAS_TABLE="diasUnattributedAfterSlowMoverPass_plusNoise10000PerImage"

DB_USER="jmyers"
DB_PASS="jmyers"
DB_HOST="localhost"


def getCursor(useSSCursor=False):
    # SSCursor is the "server-side" cursor and does not generate
    # massive memory bottlenecks when doing a huge select (e.g. to
    # fetch a large table into main memory, which we do frequently to
    # remove MySQL overhead.)
    if (useSSCursor):
        conn = db.connect(user=DB_USER, passwd=DB_PASS, host=DB_HOST, cursorclass=cursors.SSCursor)
    else:
        conn = db.connect(user=DB_USER, passwd=DB_PASS, host=DB_HOST)
    curs = conn.cursor()
    return curs
