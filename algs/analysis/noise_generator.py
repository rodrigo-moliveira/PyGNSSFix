from src.models.noise.noise_configuration import NoiseModel
import numpy as np

from src.models.noise.stats import power_spectral_density, allan_variance
from src.models.plots.utils import plot_1D, show_all, loglog

###############################
# globals / control variables #
###############################
fs = 1  # sample frequency in [Hz]
t0 = 0  # first time instance of simulation in [s]
tf = 18000  # last time instance of simulation in [s]

###############################
# User Def / Configurations  ##
###############################
units = "rad/s"  # units of the measured variable
b_rw = True  # Simulate Random Walk
b_gm = True  # Simulate Gauss Markov
b_wn = True  # Simulate White Noise
b_rc = True  # Simulate Random Constant

sigma_rw = [1,1,1]        # random walk sigma
sigma_gm = [1,1,1]        # gauss markov sigma
tau_gm = 10000          # gauss markov correlation time
sigma_wn = [1,1,1]        # white noise sigma
sigma_rc = [1,1,1]        # random constant standard deviation


def plots(name, axes, trajectory, psd, allan_var):
    # unpacking data
    time, process = trajectory
    f, power = psd
    a_var, tau = allan_var

    a_dev = np.sqrt(a_var)

    # plot timeseries
    ax = None
    for i in range(axes):
        ax = plot_1D(time, process[:, i], ax=ax, label=f"Axis {i}", x_label="Time [s]",
                     y_label=f"process [{units}]",
                     title=f"Time series of {str(name)}", set_legend=True, grid=True)

    # plot power spectral density
    ax = None
    for i in range(axes):
        ax = loglog(f, power[:, i], ax=ax, label=f"Axis {i}", x_label="Frequency [Hz]",
                    y_label=f"Power ({units})^2/Hz",
                    title=f"Power Spectral Density of {str(name)}", set_legend=True)

    # plot allan deviation
    ax = None
    for i in range(axes):
        ax = loglog(tau, a_dev[:, i], ax=ax, label=f"Axis {i}", x_label="Tau [s]",
                    y_label=f"ADEV [{units}]",
                    title=f"Allan Deviation of {str(name)}", set_legend=True)



def analyze_process(time, noise_model: NoiseModel, sampling_time):
    # process series
    n = len(time)
    process = noise_model.gen(n, sampling_time)

    # power spectral density
    power, f = power_spectral_density(process, 1/sampling_time)

    # allan variance
    a_var, tau = allan_variance(process, fs)

    plots(noise_model.name, noise_model.process_gen.axis, (time, process), (f, power), (a_var, tau))

    print(f"### Statistics of the simulated realization of process {repr(noise_model.name)} ###")
    print(f"\tmean = {process.mean(axis=0)} [{units}]")
    print(f"\tstd = {process.std(axis=0)} [{units}]")
    print(f"\tmean power = {power.mean(axis=0)} [({units})^2/Hz]")
    print("")



def main():
    # get time vector
    time = np.arange(t0, tf, step=1 / fs)

    # Gauss Markov generator
    if b_gm:
        model_gm = NoiseModel("model_gauss_markov", "gauss_markov", process_noise=sigma_gm, correlation_time=tau_gm)
        analyze_process(time, model_gm, 1/fs)

    # Random Walk generator
    if b_rw:
        model_rw = NoiseModel("model_random_walk", "random_walk", process_noise=sigma_rw)
        analyze_process(time, model_rw, 1/fs)

    # White Noise generator
    if b_wn:
        model_wn = NoiseModel("model_white_noise", "white_noise", process_noise=sigma_wn)
        analyze_process(time, model_wn, 1/fs)

    # Random Constant generator
    if b_rc:
        model_rc = NoiseModel("model_random_constant", "random_constant", process_noise=sigma_rc)
        analyze_process(time, model_rc, 1/fs)

    show_all()



print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
