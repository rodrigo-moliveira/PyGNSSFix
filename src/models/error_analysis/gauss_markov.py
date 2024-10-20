"""import numpy as np
import matplotlib.pyplot as plt

# ver a discussao https://chat.openai.com/c/ea40bb5f-1ca1-4c7d-aec6-aaaf20e3c636
# para perceber a correta discretização do continuous time gauss markov
# comparar com as fontes que tenho atualmente..

# Parameters
beta = 0.1  # Correlation constant
sigma = 0.1  # Standard deviation of process noise
initial_value = 1.0  # Initial value
T = 1.0  # Total simulation time
dt = 0.01  # Time step
num_time_points = int(T / dt)

# Initialize arrays to store process values and time points
time_points = np.linspace(0, T, num_time_points)
process_values = np.zeros(num_time_points)
process_values[0] = initial_value

# Generate the continuous-time process
for t in range(1, num_time_points):
    dW = np.random.normal(0, np.sqrt(dt))  # Wiener process increment
    dX = -beta * process_values[t-1] * dt + sigma * dW
    process_values[t] = process_values[t-1] + dX

# Plot the process
plt.plot(time_points, process_values)
plt.xlabel("Time")
plt.ylabel("Process Value")
plt.title("Continuous-Time Gauss-Markov Process")
plt.show()"""