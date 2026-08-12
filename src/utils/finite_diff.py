def finite_difference(x1, x2, step):
    if not isinstance(step, float) and not isinstance(step, int):
        raise TypeError(f"step must be either of type integer of float")

    return (x2 - x1) / step
