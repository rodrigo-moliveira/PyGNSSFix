# to plot skyplots, https://astroplan.readthedocs.io/en/latest/tutorials/plots.html#plots-sky-charts
from matplotlib import pyplot as plt

# Note -> elevation is positive towards East, that is elevation = 90 [deg] points towards East
from ..math_utils.Constants import Constant


def plot_sky(trajectory, target_name, ax=None, style_kwargs=None,
             north_to_east_ccw=True, grid=True,
             warn_below_horizon=False):
    # Note: this function is adapted from https://github.com/astropy/astroplan/blob/master/astroplan/plots/sky.py
    # but I had to adapt it to this project

    # Set up axes & plot styles if needed.
    if ax is None:
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('marker', 'o')
    style_kwargs.setdefault('label', target_name)
    style_kwargs.setdefault('s', 20)

    # Grab elevation and azimuth
    # elevation is smartly plotted through the radius dimension of the polar coordinates
    altitude = [90.1 - x[1] * Constant.RAD2DEG for x in trajectory]  # plotted in degrees
    azimuth = [x[0] for x in trajectory]  # Azimuth MUST be given to plot() in radians

    # if altitude > 91.0:
    #    if warn_below_horizon:
    #        msg = 'Target "{0}" is below the horizon at time: {1}'
    #        #warnings.warn(msg, PlotBelowHorizonWarning)

    # More axes set-up.
    # Position of azimuth = 0 (data, not label).
    ax.set_theta_zero_location('N')

    # Direction of azimuth increase. Clockwise is -1
    if north_to_east_ccw is False:
        ax.set_theta_direction(-1)

    # Plot target coordinates.
    ax.scatter(azimuth, altitude, **style_kwargs)

    # Set radial limits.
    ax.set_rlim(0, 90.1)

    # Grid, ticks & labels.
    # May need to set ticks and labels AFTER plotting points.
    if grid is True:
        ax.grid(True, which='major')
    if grid is False:
        ax.grid(False)
    degree_sign = u'\N{DEGREE SIGN}'

    # For positively-increasing range (e.g., range(1, 90, 15)),
    # labels go from middle to outside.
    r_labels = [
        '90' + degree_sign,
        '',
        '60' + degree_sign,
        '',
        '30' + degree_sign,
        '',
        '0' + degree_sign + ' Alt.',
    ]

    theta_labels = []
    for chunk in range(0, 7):
        label_angle = chunk * 45.0
        while label_angle >= 360.0:
            label_angle -= 360.0
        if chunk == 0:
            theta_labels.append('N ' + '\n' + str(label_angle) + degree_sign
                                + ' Az')
        elif chunk == 2:
            theta_labels.append('E' + '\n' + str(label_angle) + degree_sign)
        elif chunk == 4:
            theta_labels.append('S' + '\n' + str(label_angle) + degree_sign)
        elif chunk == 6:
            theta_labels.append('W' + '\n' + str(label_angle) + degree_sign)
        else:
            theta_labels.append(str(label_angle) + degree_sign)
    theta_labels.append('')

    # Set ticks and labels.
    ax.set_rgrids(range(1, 106, 15), r_labels, angle=-45)
    ax.set_thetagrids(range(0, 360, 45), theta_labels)

    # Redraw the figure for interactive sessions.
    ax.figure.canvas.draw()

    return ax
