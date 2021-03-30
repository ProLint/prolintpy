import numpy as np

FUNCS = {
    'mean' : lambda c, co, *args : c[c>co].mean(),
    'sum'  : lambda c, co, *args : c[c>co].sum(),
    'max'  : lambda c, co, *args : c[c>co].max(),
    'lnr'  : lambda c, co, t, *args : c[c>co].sum() / (t.time[-1] - t.time[0]),
    'nlnr' : lambda c, co, t, cl, *args : c[c>co].sum() / (cl * (t.time[-1] - t.time[0])),
    'occ'  : lambda c, co, t, *args :  100 * (c.sum() / t.n_frames)
}

def calculate_metric(p, func, co=0, t=None, norm=False, unit='us', *args):
    """Calculate default ProLint metrics. The function is designed such
    that it can be easily extended for other custom metrics. The second argument
    has to be a callable object, which is called using the following
    default arguments: contacts, cutoff, MDTraj.Trajectory, number_of_contacts, *args.

    By default ProLint offers the following metrics:

    ``mean`` (Mean_Duration) : The average duration of all contacts.

    ``max`` (Longest_Duration) : The longest duration (averaged if more than 1 protein).

    ``sum`` (Sum_of_all_Contacts) : The total sum of all contacts.

    ``lnr`` (Lipid_Number) : The average number of lipids in contact (the total number of contacts
    normalized with respect to the number of frames).

    ``nlnr`` (Normalized_Lipid_Number) : Lipid_Number normalized with respect to the
    number of different lipids (e.g. number of different cholesterols).

    ``occ`` (Occupancy) : For each frame, we give a value of 0 if no lipid of interest is in
    contact with the residue and 1 otherwise. Occupancy is then: sum_of_all_1s / nr_of_frames.


    Parameters
    ----------

    p : Per residue contacts.
        This is the output of ProLint.retrieve_contacts() call.

    func : callable object.
        Custom function or one of the default ones: FUNC[x] where x is
        one of the following: 'mean', 'sum', 'max', 'lnr', 'nlnr', 'occ'. If you
        want to call FUNCS['occ'] make sure to specify the argument contact='occupancy' to
        the retrieve_contacts call first.

    co : int.
        Discard contacts shorter than co value. Useful when you want to discard contacts
            with a short duration

    t : MDTraj.Trajectory

    norm : bool, default=False,
        Normalize contacts with respect to time

    unit : str
        Time unit to use for the normalization. Either 'us' or 'ns'.

    *args : additional arguments.
        These arguments are passed on to the callable object: func.

    Returns
    -------
    contacts : dict
        Dictionary where the keys are lipids and values are tuples containing
        the mean and standard deviation of the calcualted metric.


    """

    per_replicate = {}
    for lipid, v in p.items():
        r = [func(c, co, t, len(v), *args) if len(c[c>co]) > 0 else 0 for c in v]
        if norm:
            r = [x / 1000. if unit == 'ns' else x / 1000000. for x in r]
        per_replicate[lipid] = (np.mean(r), np.std(r))

    return per_replicate


def residence_time(c, ln, time, delta_t_range=None, range_type='mixed', step=None):
    """Calculate Residence Time. 
    See the following sources for a detailed explanation:
    https://pubs.acs.org/doi/abs/10.1021/ja310577u and https://www.pnas.org/content/117/14/7803

    This metric is currently not calcualted by default. It needs more testing. 
    Its output is also more difficult to control and predict. 

    Parameters
    ----------
    c : contacts
        This is the output of retrieve_contacts_flat, or it can be any flat numpy array.

    ln : int
        The number of the lipid of interest (i.e. the number of residues)

    time : float
        The total trajectory time.

    delta_t_range : list
        The range of the delta_t parameter. If None, it is calcculated.

    range_type : str
        Residence time can be sped up by supplying a default range. You can specify either a geometric,
        linear, or mixed range ('geo', 'mixed', 'linear'). Fastest and default is 'geo'. Only used
        if delta_t_range=None.

    step : int
        Custom step size for delta_t_range. Not implemented.


    """


    from scipy.optimize import curve_fit

    if delta_t_range is None:
        if range_type == 'mixed':
            delta_t_range = list(np.arange(5,100,1)) + list(np.arange(100,1000,5)) + list(np.arange(1000,time,10))
        elif range_type == 'geo':
            delta_t_range = list(np.geomspace(0.01,time,1000))
        elif range_type == 'linear':
            delta_t_range = list(np.arange(5,time,1))

    sigma = {0: 1}
    sigma0 = sum([r for r in c]) / (float(time) * ln)

    for delta_t in delta_t_range:
        denominator = (float(time) - delta_t) * ln * sigma0

        if denominator != 0:
            sigma[delta_t] = sum([r - delta_t for r in c if r >= delta_t]) / denominator
        else:
            sigma[delta_t] = 0

    try:
        popt, pcov = curve_fit(lambda t,a,b,c,d:c*np.exp(-a*t)+d*np.exp(-b*t),
                            list(sigma.keys()),
                            list(sigma.values()),
                            (1, 0.1, 1, 1))
    except RuntimeError:
        return (0, 0)

    return (1 / min([abs(k) for k in popt[:2]]), pcov)
