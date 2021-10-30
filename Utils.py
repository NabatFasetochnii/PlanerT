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
