"""
Deprecated

import numpy as np
import allantools as at
import matplotlib.pyplot as plt

# TODO: clean all this up with https://chatgpt.com/c/ce29b0ea-20a7-42e5-a84b-27e08274bd92
#    see also here https://allantools.readthedocs.io/en/latest/functions.html


Compute the Allan deviation for a given time series.

The Allan deviation is a measure of frequency stability in time series data.
This function calculates the Allan deviation using the theoretical formula,
which involves computing the time-domain differences of the input series.

Parameters:
-----------
x_series : array-like
    The time series data, which typically represents clock estimates or phase measurements.
tau : float
    The averaging time interval in seconds.

Returns:
--------
allan_dev : float
    The Allan deviation for the specified averaging time interval tau.

Procedure:
----------
1. Ensure the input time series `x_series` has at least three elements.
   If not, raise a ValueError.
2. Initialize `variance_sum` to zero. This will accumulate the sum of squared differences.
3. Loop through the series from the first element to the third-to-last element:
   a. For each element at index `k`, compute the second difference:
      `diff = x_series[k + 2] - 2 * x_series[k + 1] + x_series[k]`
   b. Square the difference and add it to `variance_sum`.
4. After the loop, compute the Allan variance using the formula:
   `allan_var = variance_sum / (2 * tau**2 * (N - 2))`
   where `N` is the length of `x_series`.
5. The Allan deviation is the square root of the Allan variance:
   `allan_dev = np.sqrt(allan_var)`



Notes:
------
- The input time series `x_series` should be evenly spaced in time.
- The averaging time interval `tau` should be a multiple of the sampling interval.
- This function assumes that `tau` is constant across the time series.

def allan_deviation(x_series, tau):
    N = len(x_series)
    if N < 3:
        raise ValueError('Time array must have at least 3 elements.')

    variance_sum = 0.0

    for k in range(N - 2):
        variance_sum += (x_series[k + 2] - 2 * x_series[k + 1] + x_series[k])**2

    allan_var = variance_sum / (2 * tau**2 * (N - 2))
    allan_dev = np.sqrt(allan_var)

    return allan_dev

def compute_allan_deviation_from_clock_bias2(y, tau0):
    # Convert clock estimates to frequency deviations
    f = np.diff(y) / tau0

    # Determine the range of tau values
    N = len(f)
    max_m = int(np.log2(N // 2))  # Maximum m based on the length of the data
    m_values = [2 ** i for i in range(max_m)]  # Tau values in powers of 2
    tau_values = [m * tau0 for m in m_values]
    allan_devs = []

    for m in m_values:
        M = N // m  # Number of windows
        if M < 2:
            break  # Need at least two windows to compute Allan deviation

        # Compute averaged frequency deviations for each window of size m
        averaged_freq = np.array([np.mean(f[i * m:(i + 1) * m]) for i in range(M)])

        # Compute the differences of averaged frequency deviations
        diff_squared = np.diff(averaged_freq) ** 2

        # Compute Allan variance and Allan deviation
        allan_variance = np.sum(diff_squared) / (2 * (M - 1))
        allan_dev = np.sqrt(allan_variance)
        allan_devs.append(allan_dev)

    return tau_values, allan_devs

def compute_allan_deviation_from_clock_bias(y, tau0):
    # Convert clock estimates to frequency deviations
    f = np.diff(y) / tau0

    # Calculate Allan deviation
    N = len(f)
    allan_devs = []
    tau_values = []

    for i in range(1, N):
        val = np.sum((f[i:] - f[:-i]) ** 2)
        tau = i * tau0
        allan_dev = np.sqrt(val / (2 * (N - i)))
        tau_values.append(tau)
        allan_devs.append(allan_dev)

    return tau_values, allan_devs

def compute_allan_deviation(f, tau0):
    N = len(f)
    max_m = int(np.log2(N // 2))  # Maximum m based on the length of the data
    m_values = [2 ** i for i in range(max_m)]  # Tau values in powers of 2
    #tau_values = [m * tau0 for m in m_values]
    allan_devs = []

    for m in m_values:
        M = N // m  # Number of windows
        if M < 2:
            break  # Need at least two windows to compute Allan deviation

        # Compute averaged frequency deviations for each window of size m
        averaged_freq = np.array([np.mean(f[i * m:(i + 1) * m]) for i in range(M)])

        # Compute the differences of averaged frequency deviations
        diff_squared = np.diff(averaged_freq) ** 2

        # Compute Allan variance and Allan deviation
        allan_variance = np.sum(diff_squared) / (2 * (M - 1))
        allan_dev = np.sqrt(allan_variance)
        allan_devs.append(allan_dev)

    return tau_values, allan_devs

# Assume y is the array of clock estimates and tau0 is the constant epoch interval
y = np.random.normal(0,1,10000)  # Your clock estimates here
y = np.array([1.5865929980152038e-07,1.5737721278969898e-07,1.5734583077920247e-07,1.6010771253197945e-07,1.6204466850601402e-07,1.5921906382662883e-07,1.6003802554301417e-07,1.6115775473329696e-07,1.5807687573008019e-07,1.5826973449298664e-07,1.6264361769938928e-07,1.6432896455820503e-07,1.6578452085366723e-07,1.6683653734646047e-07,1.6701111449691945e-07,1.6844869489397793e-07,1.6890794807406934e-07,1.680365985594066e-07,1.6854103612118545e-07,1.6798908390363825e-07,1.6692644228761198e-07,1.6763320241549798e-07,1.6511839742663347e-07,1.6449580294725763e-07,1.6537279504891555e-07,1.6577400384995336e-07,1.642145618270083e-07,1.669257886756088e-07,1.6445462243680084e-07,1.6459056535639836e-07,1.6544369975111234e-07,1.6381781858952346e-07,1.638780125433898e-07,1.6557375516973204e-07,1.655075886702462e-07,1.6209968602001227e-07,1.5968559029838404e-07,1.612961422747736e-07,1.6090696795478124e-07,1.612235307795481e-07,1.637309897355765e-07,1.6457604039086803e-07,1.628188462132856e-07,1.6165416801738302e-07,1.6027850406964268e-07,1.634981963867046e-07,1.619101828536616e-07,1.628819358334088e-07,1.6132336592840896e-07,1.5734071438992558e-07,1.5720137460918854e-07,1.5693931886815376e-07,1.5507879659070686e-07,1.5589469114315246e-07,1.5720877469011524e-07,1.552066329710234e-07,1.5728777604994702e-07,1.5965267501409393e-07,1.6349817802060255e-07,1.62853272287543e-07,1.6594704508727004e-07,1.6431167507344843e-07,1.628022067024099e-07,1.6150275924629782e-07,1.646234997883352e-07,1.649857996684897e-07,1.648013674797908e-07,1.659514265404572e-07,1.6870554095216146e-07,1.691017364624651e-07,1.7163012338650778e-07,1.7118694997374144e-07,1.6731242614593697e-07,1.6692717474557277e-07,1.6782283244345063e-07,1.6759662207475362e-07,1.7197295405428866e-07,1.6986607447353604e-07,1.709041267611451e-07,1.6805199310557514e-07,1.6484753687225545e-07,1.6731053980457564e-07,1.7707098013386034e-07,1.7743747346908896e-07,1.7901346002742933e-07,1.7482790020257864e-07,1.7220694597145135e-07,1.6890590824128443e-07,1.719143020273411e-07,1.723772495075198e-07,1.7359338202881108e-07,1.7461474777793718e-07,1.7975691429133284e-07,1.7869046557447594e-07,1.7607466064592467e-07,1.76132864468205e-07,1.7341833417135293e-07,1.7571691259333574e-07,1.760209815836855e-07,1.7659386890428642e-07,1.7798337093882295e-07,1.7858977186224775e-07,1.7668869789082957e-07,1.7273411541482508e-07,1.724731139218339e-07,1.7911083913892224e-07,1.7915674825152165e-07,1.7811041339487645e-07,1.7641465267771203e-07,1.7498466198700677e-07,1.773819941258319e-07,1.7199292633712736e-07,1.726230771529965e-07,1.7396410479371655e-07,1.788313426484016e-07,1.7876451794887529e-07,1.767286816968129e-07,1.7810370724184385e-07,1.7843898299138935e-07,1.7453906993967719e-07,1.7500269395581444e-07,1.7368609298985538e-07,1.7329235830425663e-07,1.7467525566691054e-07,1.7178110983429497e-07,1.7193630579580128e-07,1.7265677225564834e-07,1.723628757807178e-07,1.7291664106196057e-07,1.7060210077120875e-07,1.7226862267766336e-07,1.7045496644596974e-07,1.6987166262640955e-07,1.7056980551652716e-07,1.7411480313239082e-07,1.7717537201408875e-07,1.7783351747520332e-07,1.7706858186299105e-07,1.78063967790842e-07,1.7652861970337142e-07,1.7400180974729767e-07,1.7199629288611304e-07,1.7570070407854735e-07,1.8129317656618973e-07,1.812890956934967e-07,1.7880426523818913e-07,1.7826053457697842e-07,1.768657796933148e-07,1.7458226132069074e-07,1.7334225608571502e-07,1.7272807367400607e-07,1.7059544154559866e-07,1.723343688572897e-07,1.7456425977154083e-07,1.7586399439807015e-07,1.7264601462885525e-07,1.742138391487784e-07,1.7414037682609515e-07,1.713424860487204e-07,1.717070053179044e-07,1.7148045642167008e-07,1.6935806970301348e-07,1.6331085964942624e-07,1.663291783983195e-07,1.6634346654492916e-07,1.673223437590918e-07,1.6934970061189013e-07,1.6800653307358723e-07,1.704922974843724e-07,1.715173592689232e-07,1.7049973627147869e-07,1.7024508915083233e-07,1.7910809631387324e-07,1.819980471679075e-07,1.7199185452646334e-07,1.7371905040292525e-07,1.704461506427517e-07,1.7122426264365718e-07,1.6687496170493614e-07,1.6586285552999742e-07,1.6512749307189222e-07,1.6428935713458416e-07,1.6494809721701657e-07,1.6695024469532017e-07,1.698507853107525e-07,1.7247169123425604e-07,1.7204542265303858e-07,1.7589266191411945e-07,1.7873305448163617e-07,1.741529345809751e-07,1.74054943731609e-07,1.7367089384878054e-07,1.746516646565789e-07,1.70858233729292e-07,1.7406077869361006e-07,1.7428536508340022e-07,1.7392649966278007e-07,1.738129248714572e-07,1.7388615756675246e-07,1.7348834055906198e-07,1.7254326413004358e-07,1.7111323555368263e-07,1.7125085664220718e-07,1.7134468195384442e-07,1.7326177348098372e-07,1.740346624859632e-07,1.7517814605686898e-07,1.746320547698955e-07,1.7052151291566673e-07,1.7085472848347896e-07,1.7277518186784768e-07,1.7298702002990538e-07,1.7535590004603244e-07,1.7196692085627293e-07,1.7321405549968966e-07,1.7277313687168232e-07,1.7810682013979223e-07,1.8273440645133319e-07,1.795584643488636e-07,1.8104729812554385e-07,1.7894671163639735e-07,1.7604346602185613e-07,1.777123630680312e-07,1.7213183526868407e-07,1.7214776502087649e-07,1.7477900446127063e-07,1.739402050695531e-07,1.7717524277278943e-07,1.7519720484360326e-07,1.7500331093049114e-07,1.7478873234937753e-07,1.7857026150210666e-07,1.785910378027897e-07,1.7936089930389792e-07,1.7732903246078552e-07,1.770606587761163e-07,1.801851976103564e-07,1.723556626664022e-07,1.6969130047968686e-07,1.686992565363624e-07,1.708906604477454e-07,1.71656200866777e-07,1.662471992879911e-07,1.7181140005361075e-07,1.7173268800724345e-07,1.8038826094704811e-07,1.8228701155984075e-07,1.7582084451519383e-07,1.7032929676671626e-07,1.729823924398256e-07,1.753764455733759e-07,1.7624155555056482e-07,1.7811766484478342e-07,1.7685002108099093e-07,1.8205654820663324e-07,1.7956638381610115e-07,1.8233243443192228e-07,1.7936849805585365e-07,1.816908687365468e-07,1.7884317874293333e-07,1.7668947013599901e-07,1.7216028143112276e-07,1.747624124334483e-07,1.783213876201515e-07,1.7716506276434522e-07,1.786231698154006e-07,1.761937019460745e-07,1.7259172313203662e-07,1.7194435906510658e-07,1.7420073910147372e-07,1.766529604932314e-07,1.7516772218114405e-07,1.7303417819739502e-07,1.782284710587724e-07,1.7467906287307057e-07,1.785840149138961e-07,1.7805955905019566e-07,1.7458030732115404e-07,1.7377095346101607e-07,1.765430772029909e-07,1.7903472766356874e-07,1.7490727266815092e-07,1.751932119515866e-07,1.800917691488452e-07,1.7899931113454384e-07,1.8050158526496092e-07,1.8022934325016104e-07,1.8267202790044796e-07,1.800236217287817e-07,1.777574255790349e-07,1.8307090163561204e-07,1.8574992462871552e-07,1.8746736625749396e-07,1.885211032109731e-07,1.849843262964223e-07,1.818843036347098e-07,1.827344943164647e-07,1.8034486011501425e-07,1.791292287908663e-07,1.7479844238595278e-07,1.7475794355028556e-07])
tau0 = 1  # Your constant time interval

# Convert clock estimates to frequency deviations
f = np.diff(y) / tau0

# Compute Allan deviation using allantools
taus, adevs, _, _ = at.oadev(f, rate=1/tau0, data_type='freq', taus='decade')

# Plotting
plt.loglog(taus, adevs, label="oadev")
plt.xlabel('Tau (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation of GNSS Clock Estimates')
plt.grid(True)

# Compute Allan deviation using allantools
taus, adevs, _, _ = at.adev(f, rate=1/tau0, data_type='freq', taus='decade')

# Plotting
plt.loglog(taus, adevs, label="adev")
plt.xlabel('Tau (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation of GNSS Clock Estimates')
plt.grid(True)
#plt.show()

# Compute the Allan deviation
tau_values, allan_devs = compute_allan_deviation_from_clock_bias2(f, tau0)

# Plotting
plt.loglog(tau_values, allan_devs, label="compute_allan_deviation_from_clock_bias2")
plt.xlabel('Tau (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation of GNSS Clock Estimates')
plt.grid(True)


tau_values, allan_devs = compute_allan_deviation(f, tau0)

# Plotting
plt.loglog(tau_values, allan_devs, label="compute_allan_deviation")
plt.xlabel('Tau (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation of GNSS Clock Estimates')
plt.grid(True)
print(tau_values, allan_devs)

# Calculate Allan deviation for different tau values
tau_values = []
allan_devs = []
max_tau_exponent = int(np.log2(len(y) // 2))

for i in range(max_tau_exponent):
    tau = 2**i * tau0
    downsampled_series = y[::2**i]  # Downsample the series for the given tau
    allan_dev = allan_deviation(downsampled_series, tau)
    tau_values.append(tau)
    allan_devs.append(allan_dev)

# Plotting
print(tau_values, allan_devs)
plt.loglog(tau_values, allan_devs, label='From X series')
plt.xlabel('Tau (s)')
plt.ylabel('Allan Deviation')
plt.title('Allan Deviation of GNSS Clock Estimates')
plt.grid(True)
plt.legend()
plt.show()"""