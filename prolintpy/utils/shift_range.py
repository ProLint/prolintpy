def shift_range(values, new_min=0, new_max=1):
    """
    Shift values from one range to another.
    """
    old_min = min(values)
    old_max = max(values)

    old_range = (old_max - old_min)
    new_range = (new_max - new_min)
    new_val = []
    for val in values:
        try:
            new_value = (((val - old_min) * new_range) / old_range) + new_min
        except ZeroDivisionError:
            new_value = new_min
        new_val.append(new_value)
    return new_val
