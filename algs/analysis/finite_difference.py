import numpy as np
import pandas as pd

from src.utils.finite_diff import NumericalDifferentiator

def main():
    time_array = np.arange(0,10,0.001)

    functions = {
        "x": lambda a: a + 10,
        "y": lambda a: 2 * a ** 2 + 10 * a + 1,
        "z": lambda a: a ** 3,
    }
    functions_derivative = {
        "x": lambda a: 1,
        "y": lambda a: 4 * a + 10,
        "z": lambda a: 3 * a ** 2,
    }

    df = pd.DataFrame({"time": time_array})
    df_der = pd.DataFrame({"time": time_array})

    for name, function in functions.items():
        df[name] = function(time_array)
    for name, function in functions_derivative.items():
        df_der[name] = function(time_array)

    diff_mng = NumericalDifferentiator("optimal")
    for time in time_array:
        numerical_der = diff_mng.compute(df, 0,[1,2,3], float(time))
        index = np.flatnonzero(time_array == time)[0]
        analytical_der = np.asarray(df_der.iloc[index][1:4])
        print(time, "numerical_der", numerical_der, "analytical_der", analytical_der, " -> [ERROR]", analytical_der-numerical_der)



if __name__ == "__main__":
    main()

