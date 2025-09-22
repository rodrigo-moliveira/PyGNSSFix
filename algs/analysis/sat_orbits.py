""" Satellite Orbits Evaluation Script.

Script to plot the differences between precise and navigation orbits (position + velocity).
 Two tests are executed:
    1) Orbit differences: the precise and navigation orbits are evaluated at the knot points defined in the SP3 file
    2) Interpolation Error: the orbit differences are computed with interpolated precise and navigation orbits formed
                            at the defined interpolation epochs

The configurations to execute this script (path to SP3 and RINEX Nav files), and the epochs to perform the test are
defined in the following file attributes, that can be configured below

Attributes:
    sp3_file(list[str]): list of paths with the SP3 files to process
    rnx_nav_file(str): path to the RINEX NAV file
    gal_nav_type(str): select the GAL navigation message type (FNAV or INAV)
    first_epoch(str): first epoch to read from the SP3 file, in format `yyyy-mm-dd HH:MM:DD`
    last_epoch(str): last epoch to read from the SP3 file, in format `yyyy-mm-dd HH:MM:DD`
    first_interpolation_epoch(str): first epoch for the interpolation test, in format `yyyy-mm-dd HH:MM:DD`
    last_interpolation_epoch(str): last epoch for the interpolation test `yyyy-mm-dd HH:MM:DD`
    interpolation_period(int): interpolation period, in seconds
 """
import datetime
import os
import time

from src import RUNS_PATH
from src.io.rinex_parser import RinexNavReader
from src.data_mng.gnss.sat_orbit_data import SatelliteOrbits
from src.data_mng.gnss.navigation_data import NavigationData
from src.common_log import set_logs, MAIN_LOG, get_logger
from src.data_types.date import Epoch
from src.models.plots import utils

########################################
# ------- Script Configurations ------ #
# configuring sp3 and nav files
sp3_file = ["datasets/gnss/BRUX/COD0R03FIN_20190140000_01D_05M_ORB.SP3"]
rnx_nav_file = "datasets/gnss/BRUX/BRDC00IGS_R_20190140000_01D_MN.rnx"
gal_nav_type = "FNAV"  # FNAV or INAV

# Define first and last epochs for the test (optional, set to None to run the full file)
first_epoch = "2019-01-14 10:00:00"  # may be set to None
last_epoch = "2019-01-14 18:00:00"  # may be set to None

# Define interpolation interval and time period
first_interpolation_epoch = "2019-01-14 10:15:00"
last_interpolation_epoch = "2019-01-14 10:16:00"
interpolation_period = 1  # in seconds


# --- End of Script Configurations --- #
########################################


def _plot(time_vec, data_vec, x_label, y_label, title):
    ax = utils.plot_1D(time_vec, [x[0] for x in data_vec], x_label=x_label, y_label=y_label, title=title, label="X")
    utils.plot_1D(time_vec, [x[1] for x in data_vec], x_label=x_label, y_label=y_label, title=title, ax=ax, label="Y")
    utils.plot_1D(time_vec, [x[2] for x in data_vec], x_label=x_label, y_label=y_label, title=title, ax=ax, label="Z",
                  set_legend=True)
    return ax


def plots(log, plot_dir, time_vec, data_vec, product, interpolation=False):
    for sat, data in data_vec.items():
        if interpolation:
            title = f"{product[0]} Interpolation Differences (Precise Minus Navigation Interpolated Orbits) for" \
                    f" {str(sat)}"
        else:
            title = f"{product[0]} Differences (Precise Minus Navigation Orbits) for {str(sat)}"

        ax = _plot(time_vec, data, x_label="Time", y_label=f"{product[0]} error [{product[1]}]", title=title)

        plot_path = plot_dir + "\\" + f"{product[0]}_{str(sat)}" + ".png"
        log.info(f"Saving figure {plot_path}")
        ax.figure.savefig(plot_path, format='png')

    # utils.show_all()


def apply_test_interpolation(log, output_dir, sat_orbits, first_interpolation_epoch,
                             last_interpolation_epoch, interpolation_period):
    sat_list = sat_orbits.get_satellites()
    first = Epoch.strptime(first_interpolation_epoch, scale="GPST")
    last = Epoch.strptime(last_interpolation_epoch, scale="GPST")

    epoch_list = list()
    ep = first
    while ep <= last:
        epoch_list.append(ep)
        ep = ep + datetime.timedelta(seconds=interpolation_period)

    time_vec = []
    pos_vec = {}
    vel_vec = {}
    for sat in sat_list:
        pos_vec[sat] = []
        vel_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            precise_pos, precise_vel, _ = sat_orbits.get_orbit_precise(sat, epoch)
            nav_pos, nav_vel, _, _ = sat_orbits.get_orbit_broadcast(sat, epoch)

            log.info(f"[POS] {epoch}, {sat}, nav position={nav_pos}[m], precise position={precise_pos}[m] : "
                     f"diff={precise_pos - nav_pos}[m]")
            log.info(f"[VEL] {epoch}, {sat}, nav velocity={nav_vel}[m/s], precise velocity={precise_vel}[m/s] : "
                     f"diff={precise_vel - nav_vel}[m/s]")

            pos_vec[sat].append(precise_pos - nav_pos)
            vel_vec[sat].append(precise_vel - nav_vel)

    plots(log, f"{output_dir}\\interpolated_orbits", time_vec, pos_vec, ["Position", "m"], True)
    plots(log, f"{output_dir}\\interpolated_orbits", time_vec, vel_vec, ["Velocity", "m/s"], True)


def apply_test(log, output_dir, sat_orbits):
    sat_list = sat_orbits.get_satellites()
    epoch_list = list(sat_orbits.get_epochs())
    epoch_list = epoch_list[5:-5]

    time_vec = []
    pos_vec = {}
    vel_vec = {}
    for sat in sat_list:
        pos_vec[sat] = []
        vel_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            precise_pos, precise_vel, _ = sat_orbits.get_orbit_precise(sat, epoch)
            nav_pos, nav_vel, _, _ = sat_orbits.get_orbit_broadcast(sat, epoch)

            log.info(f"[POS] {epoch}, {sat}, nav position={nav_pos}[m], precise position={precise_pos}[m] : "
                     f"diff={precise_pos - nav_pos}[m]")
            log.info(f"[VEL] {epoch}, {sat}, nav velocity={nav_vel}[m/s], precise velocity={precise_vel}[m/s] : "
                     f"diff={precise_vel - nav_vel}[m/s]")

            pos_vec[sat].append(precise_pos - nav_pos)
            vel_vec[sat].append(precise_vel - nav_vel)

    plots(log, f"{output_dir}\\orbit_differences", time_vec, pos_vec, ["Position", "m"], False)
    plots(log, f"{output_dir}\\orbit_differences", time_vec, vel_vec, ["Velocity", "m/s"], False)


def create_folder():
    data_dir = str(RUNS_PATH)
    if data_dir[-1] != '//':
        data_dir = data_dir + '//'
    data_dir = data_dir + time.strftime('%Y-%m-%dT%HH%MM%SS', time.localtime()) + '_SAT_ORBITS//'
    data_dir = os.path.abspath(data_dir)

    # try to create data dir
    if not os.path.exists(data_dir):
        try:
            data_dir = os.path.abspath(data_dir)
            os.makedirs(data_dir)
            os.makedirs(f"{data_dir}\\orbit_differences")
            os.makedirs(f"{data_dir}\\interpolated_orbits")

        except:
            raise IOError(f"Cannot create dir: {data_dir}")
    return data_dir


def main():

    output_dir = create_folder()

    # initialize logger objects
    set_logs("DEBUG", f"{output_dir}\\log.txt")
    log = get_logger(MAIN_LOG)

    log.info("Executing test: Satellite Orbit Differences - SP3 vs RINEX Nav orbit differences")

    log.info("The following configurations have been selected:")
    log.info(f"RINEX CLK file: {sp3_file}\n"
             f"\t\tRINEX NAV file: {rnx_nav_file}\n"
             f"\t\tGAL Ephemeride type: {gal_nav_type}\n"
             f"\t\tFirst CLK epoch: {first_epoch}\n"
             f"\t\tLast CLK epoch: {last_epoch}\n"
             f"\t\tFirst Interpolation epoch: {first_interpolation_epoch}\n"
             f"\t\tLast Interpolation epoch: {last_interpolation_epoch}\n"
             f"\t\tInterpolation period: {interpolation_period} seconds")

    # construct data managers
    nav_data = NavigationData()
    sat_orbits = SatelliteOrbits()

    # reading the files
    log.info(f'Galileo messages selected by user are {gal_nav_type}')
    log.info('Launching RinexNavReader.')
    RinexNavReader(rnx_nav_file, nav_data, gal_nav_type)

    log.info("Launching SatelliteOrbits constructor")
    sat_orbits.init(nav_data, sp3_file, True, interp_order=9, first_epoch=first_epoch, last_epoch=last_epoch)

    log.info("Executing Satellite Orbit Test 1: (SP3 vs RINEX Nav) orbit differences")
    apply_test(log, output_dir, sat_orbits)

    log.info("Executing Satellite Orbit Test 2: (SP3 vs RINEX Nav) interpolated orbit differences")
    apply_test_interpolation(log, output_dir, sat_orbits, first_interpolation_epoch,
                             last_interpolation_epoch, interpolation_period)

    log.info("Test finished successfully!")
    log.info(f"Outputs saved to {output_dir}")


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
