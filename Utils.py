import psycopg2


def Get_Exp_Master(Vmag):
    if Vmag < 7.5:
        Exp = 0
    elif Vmag < 8.1:
        Exp = 2
    elif Vmag < 8.8:
        Exp = 3
    elif Vmag < 9.3:
        Exp = 5
    elif Vmag < 9.6:
        Exp = 10
    elif Vmag < 10.0:
        Exp = 20
    elif Vmag < 10.6:
        Exp = 30
    elif Vmag < 11.1:
        Exp = 50
    elif Vmag < 11.5:
        Exp = 80
    elif Vmag < 12.0:
        Exp = 120
    else:
        Exp = 180
    return Exp


def Check_AzEl(Coo, Az0):
    Ra, Dec = Coo.split(' ')
    Dec = Dec.split(':')
    Dec = int(Dec[0]) + int(Dec[1]) / 60. + float(Dec[2]) / 3600.
    if (Az0 >= 0) & (Az0 <= 180) & (Dec < 20):
        return 1
    else:
        return 0


def Get_Coo(Coo):
    Ra, Dec = Coo.split(' ')
    Ra = Ra.split(':')
    Ra = int(Ra[0]) + int(Ra[1]) / 60. + float(Ra[2]) / 3600.
    Ra = Ra * 15.
    Dec = Dec.split(':')
    Dec = int(Dec[0]) + int(Dec[1]) / 60. + float(Dec[2]) / 3600.
    return Ra, Dec


def Get_Times(t0, t1, t2, t3, T_Start, T_Stop):
    if t0 < T_Start.jd:
        t0 = T_Start.jd[0]
    if t1 < T_Start.jd:
        t1 = T_Start.jd[0]
    if t2 > T_Stop.jd:
        t2 = T_Stop.jd[0]
    if t3 > T_Stop.jd:
        t3 = T_Stop.jd[0]
    return t0, t1, t2, t3


def Check_El(ObsAltAz, TrAltAz):
    ObsAltAz.remove_rows(ObsAltAz['El'] < 30)
    TrAltAz.remove_rows(TrAltAz['El'] < 30)
    return ObsAltAz, TrAltAz


def Send_Task(Date_Obs, duration, alpha, delta, exptime, frame_numbers, object, filter_e='V', filter_w='R'):
    ##Master_DB params
    db_host = "192.168.240.27"
    db_name = "imdata"
    db_user = "postgres"
    db_pass = 'imdata'
    try:
        # Date_Obs = Date_Obs.split(sep='.')[0]
        ##open connection
        conn = psycopg2.connect(host=db_host, database=db_name, user=db_user, password=db_pass)
        if conn:
            cur = conn.cursor()
            cur.execute("SELECT max(task) FROM fields")
            task = cur.fetchone()[0]
            # sql = "INSERT INTO cloud_cam (timestamp, stars, maglimit, skybg, mean_ext, std_ext) \
            #                       VALUES (\'%s\', %d, %0.1f, %0.1f, %0.2f, %0.2f);" % \
            #       (Date_Obs, NStars, maglimit, skybg, mean_ext, std_ext)

            sql = f"INSERT INTO fields (task, runtime, fintime, coord2000, exptime, filter1, filter2, " +\
                  "direction, expnumb, ftype, observer, objname, focus, xbin, ybin, xstart, ystart, xstop, ystop, " + \
                  "separate, minalt,sequence) VALUES " +\
                  f"({task}, '{Date_Obs}'::timestamp, '{Date_Obs}'::timestamp+" +\
                  f"'{duration} hours'::interval, '({alpha} d,{delta} d)'::spoint, {exptime}, '{filter_e}', '{filter_w}', true, {frame_numbers}, " +\
                  f"'PHOTOMETRY', 'N.Chazov', '{object}', false::boolean, 1, 1, 0, 0, 4098, 4098, " + \
                  f"'0 minutes'::interval, 0,1)  RETURNING id; "
            cur.execute(sql)
            conn.commit()
            cur.close()
        conn.close()
    ##        print('DB ok!')
    except Exception as e:
        print('Error in db req')
        print(e)
    return
