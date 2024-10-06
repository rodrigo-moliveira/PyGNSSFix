""" Script to plot the differences between precise and navigation clocks.
 Two tests are executed:
    1) Clock differences: the precise and navigation clocks are evaluated at the knot points defined in the CLK file
    2) Interpolation Error: the clocks differences are computed with interpolated precise and navigation clocks formed
                            at the defined interpolation epochs

The configurations for the script (path to CLK and RINEX Nav files), and the epochs to perform the test are
defined in the main function.
 """
import datetime
import time
import os

from src import RUNS_PATH
from src.data_types.date import Epoch
from src.io.rinex_parser import RinexNavReader
from src.data_mng.gnss.sat_clock_data import SatelliteClocks
from src.data_mng.gnss.navigation_data import NavigationData
from src.common_log import set_logs, MAIN_LOG, get_logger
from src.models.plots import utils


def plots(log, time_vec, data_vec, plot_dir, interpolation=False):
    for sat, data in data_vec.items():
        if interpolation:
            title = f"Clock Interpolation Differences (Precise Minus Navigation Interpolated Clocks) for" \
                    f" {str(sat)}"
        else:
            title = f"Clock Differences (Precise Minus Navigation) for {str(sat)}"
        ax = utils.plot_1D(time_vec, data, x_label="Time", y_label="Clock error [s]",
                           title=title)

        plot_path = plot_dir + "\\" + str(sat) + ".png"
        log.info(f"Saving figure {plot_path}")
        ax.figure.savefig(plot_path, format='png')

    # utils.show_all()


def apply_test_interpolation(log, output_dir, sat_clocks, first_interpolation_epoch,
                             last_interpolation_epoch, interpolation_period):
    sat_list = sat_clocks.get_satellites()

    first = Epoch.strptime(first_interpolation_epoch, scale="GPST")
    last = Epoch.strptime(last_interpolation_epoch, scale="GPST")

    epoch_list = list()
    ep = first
    while ep <= last:
        epoch_list.append(ep)
        ep = ep + datetime.timedelta(seconds=interpolation_period)

    time_vec = []
    data_vec = {}
    for sat in sat_list:
        data_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            nav_clock, _ = sat_clocks.get_clock_broadcast(sat, epoch)
            precise_clock = sat_clocks.get_clock_precise(sat, epoch)

            # compare
            log.info(f"{epoch}, {sat}, nav_clock={nav_clock}[s], precise_clock={precise_clock}[s] : "
                     f"diff={precise_clock - nav_clock}[s]")

            data_vec[sat].append(precise_clock - nav_clock)

    plots(log, time_vec, data_vec, f"{output_dir}\\interpolated_clocks", True)


def apply_test(log, output_dir, sat_clocks):
    sat_list = sat_clocks.get_satellites()
    epoch_list = sat_clocks.get_epochs()

    time_vec = []
    data_vec = {}
    for sat in sat_list:
        data_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            nav_clock, _ = sat_clocks.get_clock_broadcast(sat, epoch)
            precise_clock = sat_clocks.get_clock_precise(sat, epoch)

            # compare
            log.info(f"{epoch}, {sat}, nav_clock={nav_clock}[s], precise_clock={precise_clock}[s] : "
                     f"diff={precise_clock - nav_clock}[s]")

            data_vec[sat].append(precise_clock - nav_clock)

    plots(log, time_vec, data_vec, f"{output_dir}\\clock_differences", False)


def create_folder():
    data_dir = str(RUNS_PATH)
    if data_dir[-1] != '//':
        data_dir = data_dir + '//'
    data_dir = data_dir + time.strftime('%Y-%m-%dT%HH%MM%SS', time.localtime()) + '_SAT_CLOCKS//'
    data_dir = os.path.abspath(data_dir)

    # try to create data dir
    if not os.path.exists(data_dir):
        try:
            data_dir = os.path.abspath(data_dir)
            os.makedirs(data_dir)
            os.makedirs(f"{data_dir}\\clock_differences")
            os.makedirs(f"{data_dir}\\interpolated_clocks")

        except:
            raise IOError(f"Cannot create dir: {data_dir}")
    return data_dir


def main():
    ########################################
    # ------- Script Configurations ------ #
    # configuring clock and nav files
    rnx_clock_file = ["datasets/gnss/BRUX/COD0R03FIN_20190140000_01D_30S_CLK.CLK"]
    rnx_nav_file = "datasets/gnss/BRUX/BRDC00IGS_R_20190140000_01D_MN.rnx"
    gal_nav_type = "FNAV"  # FNAV or INAV

    # Define first and last epochs for the test (optional, set to None to run the full file)
    first_epoch = "2019-01-14 10:00:00"  # may be set to None
    last_epoch = "2019-01-14 11:00:00"  # may be set to None

    # Define interpolation interval and time period
    first_interpolation_epoch = "2019-01-14 10:15:00"
    last_interpolation_epoch = "2019-01-14 10:16:00"
    interpolation_period = 1  # in seconds

    # --- End of Script Configurations --- #
    ########################################

    output_dir = create_folder()

    # initialize logger objects
    set_logs("DEBUG", f"{output_dir}\\log.txt")
    log = get_logger(MAIN_LOG)

    log.info("Executing test: Satellite Clock Differences - SP3 vs RINEX Nav clock differences")

    # construct data managers
    nav_data = NavigationData()
    sat_clocks = SatelliteClocks()

    # reading the files
    log.info(f'Galileo messages selected by user are {gal_nav_type}')
    log.info('Launching RinexNavReader.')
    RinexNavReader(rnx_nav_file, nav_data, gal_nav_type)

    log.info("Launching SatelliteClocks constructor")
    sat_clocks.init(nav_data, rnx_clock_file, True, interp_order=1, first_epoch=first_epoch, last_epoch=last_epoch)

    log.info("Executing Satellite Clock Test 1: (RINEX CLOCK vs RINEX NAV) clock differences")
    apply_test(log, output_dir, sat_clocks)

    log.info("Executing Satellite Clock Test 2: (RINEX CLOCK vs RINEX NAV) interpolated clock differences")
    apply_test_interpolation(log, output_dir, sat_clocks, first_interpolation_epoch, last_interpolation_epoch,
                             interpolation_period)

    log.info("Test finished successfully!")
    log.info(f"Outputs saved to {output_dir}")


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
