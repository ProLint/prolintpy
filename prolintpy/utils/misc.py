import mdtraj as md
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial import ConvexHull

def save_pdb(t, protein, filename, bfactors=None):
    """Save a PDB file. It supports adding contact information
    in a beta-factor column.

    Parameters
    ----------

    t : MDTraj.Trajectory

    protein : ProLint.Protein

    filename : str

    bfactors : array

    """

    df = protein.dataframe[0]
    if protein.resolution == 'atomistic':
        indices = df[(df.name == "CA")].index.to_numpy()
    elif protein.resolution == 'martini':
        indices = df[(df.name == "BB")].index.to_numpy()

    t = t[0].atom_slice(indices)

    t.save_pdb(filename, bfactors=bfactors)


def unit_poly_verts(theta, centre ):
    """Return vertices of polygon for subplot axes.
    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [centre ] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

def radar_patch(r, theta, centre ):
    """ Returns the x and y coordinates corresponding to the magnitudes of
    each variable displayed in the radar plot
    """
    # offset from centre of circle
    offset = 0
    yt = (r*centre + offset) * np.sin(theta) + centre
    xt = (r*centre + offset) * np.cos(theta) + centre
    return list(xt), list(yt)



def star_curv(old_x, old_y):
    """ Interpolates every point by a star-shaped curv. It does so by adding
    "fake" data points in-between every two data points, and pushes these "fake"
    points towards the center of the graph (roughly 1/4 of the way).
    """

    def midpoint(p1, p2, sf=1):
        xm = ((p1[0]+p2[0])/2) * sf
        ym = ((p1[1]+p2[1])/2) * sf
        return (xm, ym)

    try:
        points = np.array([old_x, old_y]).reshape(7, 2)
        hull = ConvexHull(points)
        x_mid = np.mean(hull.points[hull.vertices,0])
        y_mid = np.mean(hull.points[hull.vertices,1])
    except:
        x_mid = 0.5
        y_mid = 0.5

    c=1
    x, y = [], []
    for i, j in zip(old_x, old_y):
        x.append(i)
        y.append(j)
        try:
            xm_i, ym_i = midpoint((i, j),
                                midpoint((i, j), (x_mid, y_mid)))

            xm_j, ym_j = midpoint((old_x[c], old_y[c]),
                                midpoint((old_x[c], old_y[c]), (x_mid, y_mid)))

            xm, ym = midpoint((xm_i, ym_i), (xm_j, ym_j))
            x.append(xm)
            y.append(ym)
            c += 1
        except IndexError:
            break


    orig_len = len(x)
    x = x[-3:-1] + x + x[1:3]
    y = y[-3:-1] + y + y[1:3]


    t = np.arange(len(x))
    ti = np.linspace(2, orig_len + 1, 10 * orig_len)

    kind='quadratic'
    xi = interp1d(t, x, kind=kind)(ti)
    yi = interp1d(t, y, kind=kind)(ti)

    return xi, yi





