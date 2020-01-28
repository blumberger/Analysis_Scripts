import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import RDF


def create_coords(length):
    """
    Creates random coordinates
    """
    return np.random.random((length, 3))


def plot_3D_coords(coords):
    """
    Will plot an array of 3D coordinates
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'ko')
    plt.show()
    
    return fig, ax
        

def plot_random_rdf():
    """
    Will plot an rdf distribution against the radius for some random numbers.
    """    
    rand_coords = np.random.random((6000, 3))

    r, rdf = RDF.calc_RDF(rand_coords)      

    plt.plot(r, rdf)
    plt.title("g(r) of random points")
    plt.xlabel("R")
    plt.ylabel("g(R)")
    plt.show()


def plot_ordered_rdf():
    """
    Will plot the rdf of structured data.
    """
    x, y, z = np.arange(0, 10, 1), np.arange(0, 10, 1), np.arange(0, 10, 1)
    crds = np.meshgrid(x, y, z)
    crds = [[crds[xyz][i, j, k] for xyz in range(len(crds))]
            for i in range(len(y))
            for j in range(len(x))
            for k in range(len(z))]
    crds = np.array(crds)

    r, rdf = RDF.calc_RDF(crds) 
    
    plt.plot(r, rdf)
    plt.title("g(r) of points arranged in a grid")
    plt.xlabel("R")
    plt.ylabel("g(R)")
    plt.show()
    
    
def plot_3D_crds(crds):
    """
    Will plot some 3D coordinates, in the normal format (num_atom, 3)...
    """
    f = plt.figure()
    a = f.add_subplot(111, projection="3d")
    a.plot(crds[:, 0], crds[:, 1], crds[:, 2], 'ko')
    plt.show()


plot_ordered_rdf()
plot_random_rdf()
