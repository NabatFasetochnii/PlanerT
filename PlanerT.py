import os
import warnings

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy.utils import iers
from TTF import TTF_Query  # (DateObs, min_depth, Target)
from Utils import Get_Exp_Master, Check_AzEl, Get_Coo, Get_Times, Check_El
from matplotlib.backends.backend_pdf import PdfPages

# mpl.rcParams["figure.figsize"] = [6.4, 4.8]

mpl.use('Agg')
iers.conf.auto_download = False
warnings.simplefilter("ignore")
#################################################################
min_depth = 3.0
max_mag = 16.
Kou = EarthLocation(lat=57.036537 * u.deg, lon=59.545735 * u.deg, height=290 * u.m)

#################################################################
Local_Time = Time.now()
# Test = TimeDelta(5*24*3600, format='sec')
# Local_Time = Local_Time + Test
TimeZone = TimeDelta(5 * 3600, format='sec')
Start = Time(np.round(Local_Time.jd), format='jd') - TimeZone

#################################################################
# get Sun and Moon
Plan = Table()
Step = TimeDelta(180, format='sec')
Epochs = Time(np.arange(Start.jd, Start.jd + 1., Step.jd), format='jd')
Plan['JD'] = Epochs

altazframe = AltAz(obstime=Plan['JD'], location=Kou)
sunaltaz = get_sun(Epochs).transform_to(altazframe)
Plan['Sun_El'] = sunaltaz.alt
Plan['Sun_Az'] = sunaltaz.az
moonaltaz = get_moon(Epochs).transform_to(altazframe)
Plan['Moon_El'] = moonaltaz.alt
Plan['Moon_Az'] = moonaltaz.az

Plan['Time'] = (Plan['JD'].jd - Plan['JD'][0].jd) * 24.  # #hours after midday

Zero = np.where(Plan['Sun_El'] <= 0)
Sunset = np.min(Plan['JD'][Zero].jd)
Sunrise = np.max(Plan['JD'][Zero].jd)
Title = 'Sunset: ' + Time(Sunset, format='jd').datetime.strftime('%Y-%m-%d %H:%M')
Title += ', '
Title += 'sunrise:' + Time(Sunrise, format='jd').datetime.strftime('%Y-%m-%d %H:%M')
Min_lim = np.min(Plan['Time'][Zero])
Max_lim = np.max(Plan['Time'][Zero])

# #get TESS target list
DT = Local_Time.datetime.strftime('%m-%d-%Y')
Tess_List = TTF_Query(DT, min_depth, max_mag, '')

figHeight = 7.5
# figHeight = 12 + 0.1 * len(Tess_List)
fig, axs = plt.subplots(2, 1, figsize=(8, figHeight), dpi=125)

num = len(Tess_List)
biba = 20
print('number of objects is', num)
Name = 'PlanT.pdf'
# AttributeError: 'bool' object has no attribute 'any'

if (Tess_List is not None) & (len(Tess_List) > 0):

    Title = Title + '\nTESS transits with depth>' + str(min_depth) + 'mmag and Vmag<' + str(max_mag)

    if os.path.isfile(Name):
        os.remove(Name)
    if os.path.isfile('comments.txt'):
        os.remove('comments.txt')

    with PdfPages(Name) as pdf:

        comments_text_file = open('comments.txt', 'a')

        for row in range(0, len(Tess_List), biba):

            axs[0].cla()
            axs[1].cla()
            try:
                Night = np.where(Plan['Sun_El'] <= -18)
                N_Start = np.min(Plan['Time'][Night])
                N_Stop = np.max(Plan['Time'][Night])
                axs[0].axvspan(N_Start, N_Stop, facecolor='k', alpha=0.3)
            except:
                pass

            try:
                Twilight = np.where(Plan['Sun_El'] < -12)
                T_Start = np.min(Plan['Time'][Twilight])
                T_Stop = np.max(Plan['Time'][Twilight])
                axs[0].axvspan(T_Start, T_Stop, facecolor='k', alpha=0.2)
                T_Start = Plan['JD'][Plan['Time'] == T_Start]
                T_Stop = Plan['JD'][Plan['Time'] == T_Stop]
            except:
                pass
            # axs[0].set_position([0.125, 0.45, 0.8, 0.49])
            axs[0].set_position([0.125, 0.55, 0.8, 0.39])
            # axs[0].set_position([0.125, figHeight/(14.+(figHeight-6)), 0.8, figHeight/(25.+(figHeight-6))])
            axs[0].hlines(-18, 0, 24, 'k', 'dashed', alpha=0.5)
            axs[0].hlines(-12, 0, 24, 'k', 'dashed', alpha=0.2)
            axs[0].hlines(30, 0, 24, 'r', 'dashed', alpha=0.3)

            axs[0].plot(Plan['Time'], Plan['Sun_El'], 'r-.', alpha=0.5, label='Sun')
            axs[0].plot(Plan['Time'], Plan['Moon_El'], 'b-.', alpha=0.5, label='Moon')

            #################################################################
            collabel = ('Target', 'Coord', 'Start(UTC)', 'Duration(h)', 'Depth(mmag)',
                        'Vmag', 'Exp(s)', '2Moon(d)', 'Priority')
            full_w = 22.5
            widths = ([3 / full_w, 5 / full_w, 4 / full_w, 2 / full_w, 2 / full_w,
                       1.5 / full_w, 1.5 / full_w, 2 / full_w, 1.5 / full_w])
            axs[1].set_position([0.09, -0.025, 0.85, 0.52])
            # axs[1].set_position([0.125, 0.05, 0.8, 0.35])
            axs[1].axis('tight')
            axs[1].axis('off')

            Data = []
            colors = []
            comments = []

            if (num - row - 1) < biba:
                delim = num - row - 1
            else:
                delim = biba

            for i in range(delim):
                # times
                t1 = Tess_List['jd_start'][row + i] + 2450000
                t2 = Tess_List['jd_end'][row + i] + 2450000
                t0 = t1 - 0.5 / 24.
                t3 = t2 + 0.5 / 24.
                # check twilight
                t0, t1, t2, t3 = Get_Times(t0, t1, t2, t3, T_Start, T_Stop)
                ObsTime = Time(np.arange(t0, t3, Step.jd), format='jd')
                TrTime = Time(np.arange(t1, t2, Step.jd), format='jd')
                if (len(TrTime) > 10) & (len(ObsTime) > 10):
                    # get position
                    Ra, Dec = Get_Coo(Tess_List['coords(J2000)'][row + i])
                    Eq = SkyCoord(ra=Ra * u.deg, dec=Dec * u.deg, frame='icrs')
                    ObsAltAz = Eq.transform_to(AltAz(obstime=ObsTime, location=Kou))
                    ObsAltAz = Table([ObsAltAz.obstime, ObsAltAz.alt, ObsAltAz.az], names=('obstime', 'El', 'Az'))
                    TrAltAz = Eq.transform_to(AltAz(obstime=TrTime, location=Kou))
                    TrAltAz = Table([TrAltAz.obstime, TrAltAz.alt, TrAltAz.az], names=('obstime', 'El', 'Az'))

                    # check elevation
                    ObsAltAz, TrAltAz = Check_El(ObsAltAz, TrAltAz)

                    # check moon
                    Moon = get_moon(Time((t0 + t3) / 2., format='jd'))
                    sep = Moon.separation(Eq).degree
                    sep = np.round(sep, 2)

                    colors.append(['w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'])
                    if Check_AzEl(Tess_List['coords(J2000)'][row + i], Tess_List['az_start'][row + i]):
                        colors[-1][1] = 'lightpink'  # bad master angle
                        # status[row] = status[row] + 1
                    Exp = Get_Exp_Master(Tess_List['V'][row + i])
                    if Exp <= 2:
                        colors[-1][6] = 'lightpink'  # short exposure
                        # status[row] = status[row] + 1
                    if Tess_List['V'][row + i] > 14:
                        colors[-1][5] = 'lightpink'  # faint star
                        # status[row] = status[row] + 1
                    if Tess_List['depth(mmag)'][row + i] < 6:
                        colors[-1][4] = 'lightpink'  # min depth
                        # status[row] = status[row] + 1
                    if sep < 90:
                        colors[-1][7] = 'lightpink'  # moon distance
                        # status[row] = status[row] + 1
                    if Tess_List['priority'][row + i] < 3:
                        colors[-1][8] = 'lightgreen'  # priority
                        # status[row] = status[row] - 1
                    elif Tess_List['priority'][row + i] > 3:
                        colors[-1][8] = 'lightpink'  # priority
                        # status[row] = status[row] + 1

                    Begin = ObsAltAz['obstime'][0].datetime.strftime('%Y-%m-%d %H:%M:%S')
                    Duration = round((ObsAltAz['obstime'][-1].jd - ObsAltAz['obstime'][0].jd) * 24, 2)
                    Data.append([Tess_List['Name'][row + i].replace('.', '_').replace(' ', ''),
                                 Tess_List['coords(J2000)'][row + i],
                                 Begin,
                                 Duration,
                                 Tess_List['depth(mmag)'][row + i],
                                 Tess_List['V'][row + i],
                                 Get_Exp_Master(Tess_List['V'][row + i]),
                                 sep,
                                 Tess_List['priority'][row + i]])
                    comments_text_file.write(str((Tess_List['Name'][row + i] + ': ' + Tess_List['comments'][row + i])
                                                 .encode('ascii', 'ignore')) + '\n')
                    D = axs[0].plot((ObsAltAz['obstime'].jd - Plan['JD'][0].jd) * 24.,
                                    ObsAltAz['El'],
                                    linestyle='dotted', linewidth=2, label='_nolegend_')

                    axs[0].plot((TrAltAz['obstime'].jd - Plan['JD'][0].jd) * 24.,
                                TrAltAz['El'], linestyle='solid',
                                label=Tess_List['Name'][row + i],
                                color=D[-1].get_color(), linewidth=2)
                #        colors[-1][0] = D[-1].get_color()

            the_table = axs[1].table(cellText=Data, colLabels=collabel,
                                     colWidths=widths, cellColours=colors,
                                     loc='upper center')
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(6)

            axs[0].set_xticks(np.linspace(0, 23, 24))
            axs[0].set_xticklabels(['07', '08', '09', '10', '11', '12', '13', '14', '15',
                                    '16', '17', '18', '19', '20', '21', '22', '23', '00',
                                    '01', '02', '03', '04', '05', '06'])
            axs[0].tick_params(axis='both', labelsize=6, direction='in')
            axs[0].set_xlim(Min_lim, Max_lim)
            axs[0].set_xlabel('UTC, ' + Start.datetime.strftime('%Y-%m-%d'), fontsize=6)
            axs[0].set_ylim(-20, None)
            axs[0].set_ylabel('Elevation (deg)', fontsize=6)
            axs[0].legend(loc=3, fontsize=6)
            axs[0].grid()

            fig.suptitle(Title, fontsize=8)
            pdf.savefig()
# ascii.write(comments, 'comments.txt', fast_writer=False,
#             overwrite=True)
    comments_text_file.close()
else:
    Title = 'List is empty for some reasons. Maybe night starts/ends at nautical twilight'
    print(Title)

#################################################################
# axs[0].set_xticks(np.linspace(0, 23, 24))
# axs[0].set_xticklabels(['07', '08', '09', '10', '11', '12', '13', '14', '15',
#                         '16', '17', '18', '19', '20', '21', '22', '23', '00',
#                         '01', '02', '03', '04', '05', '06'])
# axs[0].tick_params(axis='both', labelsize=6, direction='in')
# axs[0].set_xlim(Min_lim, Max_lim)
# axs[0].set_xlabel('UTC, ' + Start.datetime.strftime('%Y-%m-%d'), fontsize=6)
# axs[0].set_ylim(-20, None)
# axs[0].set_ylabel('Elevation (deg)', fontsize=6)
# axs[0].legend(loc=3, fontsize=6)
# axs[0].grid()
#
# fig.suptitle(Title, fontsize=8)

# Path2Save = str(DT)
# if not os.path.exists(Path2Save):
#     os.makedirs(Path2Save)

# plt.savefig(Name)
