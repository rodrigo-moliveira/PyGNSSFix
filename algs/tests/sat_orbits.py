import datetime

from src.io.rinex_parser import RinexNavReader
from src.data_mng.gnss.sat_orbit_data import SatelliteOrbits
from src.data_mng.gnss.navigation_data import NavigationData
from src.data_types.gnss.satellite import get_satellite
from src.common_log import set_logs, MAIN_LOG, get_logger
from src.modules.performance.plot_manager import plot_1D, show_all
from src.data_types.date.date import Epoch


def plots(time_vec, data_vec):
    for sat, data in data_vec.items():
        plot_1D(time_vec, [x[0] for x in data], x_label="Time", y_label="Position error [m]",
                title=f"Orbit differences in X for {str(sat)}")
        plot_1D(time_vec, [x[1] for x in data], x_label="Time", y_label="Position error [m]",
                title=f"Orbit differences in Y for {str(sat)}")
        plot_1D(time_vec, [x[2] for x in data], x_label="Time", y_label="Position error [m]",
                title=f"Orbit differences in Z for {str(sat)}")
    show_all()


def apply_test_interpolation(log, sat_orbits, sats_to_plot):
    first_epoch = Epoch.strptime("2019-01-14 10:00:00", scale="GPST")
    sat_list = sat_orbits.get_satellites()
    epoch_list = [first_epoch + datetime.timedelta(seconds=60*i) for i in range(100)]

    time_vec = []
    data_vec = {}
    for sat in sats_to_plot:
        data_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            precise_orbit = sat_orbits.get_orbit_precise(sat, epoch)
            nav_orbit = sat_orbits.get_orbit_broadcast(sat, epoch)

            # compare
            #log.info(f"{epoch}, {sat}, nav_orbit={nav_orbit}[m], precise_orbit={precise_orbit}[m] : "
            #         f"diff={precise_orbit - nav_orbit}[m]")

            if sat in sats_to_plot:
                data_vec[sat].append(precise_orbit - nav_orbit)

    plots(time_vec, data_vec)

def apply_test(log, sat_orbits, sats_to_plot):
    sat_list = sat_orbits.get_satellites()
    epoch_list = sat_orbits.get_epochs()

    time_vec = []
    data_vec = {}
    for sat in sats_to_plot:
        data_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            precise_orbit = sat_orbits.get_orbit_precise(sat, epoch)
            nav_orbit = sat_orbits.get_orbit_broadcast(sat, epoch)

            # compare
            #log.info(f"{epoch}, {sat}, nav_orbit={nav_orbit}[m], precise_orbit={precise_orbit}[m] : "
            #         f"diff={precise_orbit - nav_orbit}[m]")

            if sat in sats_to_plot:
                data_vec[sat].append(precise_orbit - nav_orbit)

    plots(time_vec, data_vec)


def main():
    # initialize logger objects
    set_logs("DEBUG")
    log = get_logger(MAIN_LOG)

    log.info("Satellite Orbits Tests - (SP3 vs RINEX NAV orbit differences")

    # configuring sp3 and nav files
    sp3_file = ["datasets/gnss/BRUX/COD0R03FIN_20190140000_01D_05M_ORB.SP3"]
    rnx_nav_file = "datasets/gnss/BRUX/BRDC00IGS_R_20190140000_01D_MN.rnx"
    gal_nav_type = "FNAV"  # FNAV or INAV
    sat_lst = ["E01", "E02",
               "G01", "G02"]
    sats_to_plot = [get_satellite(x) for x in sat_lst]  # select satellites to plot

    # construct data managers
    nav_data = NavigationData()
    sat_orbits = SatelliteOrbits()

    # reading the files
    log.info(f'Galileo messages selected by user are {gal_nav_type}')
    log.info('Launching RinexNavReader.')
    RinexNavReader(rnx_nav_file, nav_data, gal_nav_type)

    log.info("Launching SatelliteOrbits constructor")
    sat_orbits.init(nav_data, sp3_file, True, interp_order=9, first_epoch="2019-01-14 02:00:00",
                    last_epoch="2019-01-14 22:00:00")

    apply_test_interpolation(log, sat_orbits, sats_to_plot)

    log.info("Test finished successfully")


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
