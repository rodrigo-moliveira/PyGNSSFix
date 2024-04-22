from src.io.rinex_parser import RinexNavReader
from src.data_mng.gnss.sat_orbit_data import SatelliteOrbits
from src.data_mng.gnss.navigation_data import NavigationData
from src.data_types.gnss.satellite import get_satellite
from src.common_log import set_logs, MAIN_LOG, get_logger
from src.modules.performance.plot_manager import plot_1D, show_all


def plots(time_vec, data_vec):
    for sat, data in data_vec.items():
        plot_1D(time_vec, data, x_label="Time", y_label="Clock error [s]",
                title=f"Clock differences for {str(sat)}")
    show_all()


def apply_test(log, sat_clocks, sats_to_plot):
    sat_list = sat_clocks.get_satellites()
    epoch_list = sat_clocks.get_epochs()

    time_vec = []
    data_vec = {}
    for sat in sats_to_plot:
        data_vec[sat] = []

    for epoch in epoch_list:
        time_vec.append(epoch)
        for sat in sat_list:
            nav_clock = sat_clocks.get_clock_broadcast(sat, epoch)
            precise_clock = sat_clocks.get_clock_precise(sat, epoch)

            # compare
            log.info(f"{epoch}, {sat}, nav_clock={nav_clock}[s], precise_clock={precise_clock}[s] : "
                     f"diff={precise_clock - nav_clock}[s]")

            if sat in sats_to_plot:
                data_vec[sat].append(precise_clock - nav_clock)

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
    sat_lst = ["E01"]#, "E02", "E03", "E04", "E05", "E30", "E31", "E33",
               #"G01", "G02", "G03", "G04", "G05", "G06", "G30", "G31", "G32"]
    sats_to_plot = [get_satellite(x) for x in sat_lst]  # select satellites to plot

    # construct data managers
    nav_data = NavigationData()
    sat_orbits = SatelliteOrbits()

    # reading the files
    log.info(f'Galileo messages selected by user are {gal_nav_type}')
    log.info('Launching RinexNavReader.')
    #RinexNavReader(rnx_nav_file, nav_data, gal_nav_type)

    log.info("Launching SatelliteOrbits constructor")
    sat_orbits.init(nav_data, sp3_file, True, "2019-01-14 02:00:00", "2019-01-14 20:00:00")

    apply_test(log, sat_orbits, sats_to_plot)

    log.info("Test finished successfully")


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
