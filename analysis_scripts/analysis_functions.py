import numpy as np
from scipy import sparse
from synmorph.analysis import geometrical as geo
from synmorph.analysis import topological as top
from synmorph.utils import *
import pickle


def get_edges_from_tri(tri):
    """
    Given a triangulation (nv x 3), get a (3nv x 2) list of cell ids defining the edges in the triangulation/tesselation
    """
    return np.row_stack([np.column_stack((tri[:,i],tri[:,(i+1)%3])) for i in range(3)])


def calculate_t1s_frame_to_frame_by_ctype(edges0,edges1,c_types):
    """
    Find the number of t1s by each cell_type (c_type) between a pair of frames
    edges0 is the list of edges at the first timepoint. Likewise edges1 is that for the second time point.
    c_types: indices of cell types for each cell index.
    """
    n_ctypes = np.max(c_types)+1
    n_t1s = np.zeros((n_ctypes,n_ctypes),dtype=np.int32)
    if np.any(edges0!=edges1):
        adj_mat0 = sparse.coo_matrix((np.ones_like(edges0[:,0],dtype=bool),(edges0[:,0],edges0[:,1])),shape=(np.max(edges0)+1,np.max(edges0)+1))
        adj_mat1 = sparse.coo_matrix((np.ones_like(edges1[:,0],dtype=bool),(edges1[:,0],edges1[:,1])),shape=(np.max(edges0)+1,np.max(edges0)+1))

        t1_mat = (adj_mat0!=adj_mat1)
        for i in range(n_ctypes):
            for j in range(n_ctypes):
                n_t1s[i,j] = int(t1_mat[c_types == i].T[c_types == j].sum()/2) ##NOTE THE DENOMINATOR IS DIFFERENT. This is because it is measuring how many pairs of cells gain or lose contact.
    return n_t1s


def calculate_t1s_across_sim_by_ctype(tri_save,c_types):
    """
    Performs the t1 by cell-type calculation across the simulation
    """
    edges_save = np.array(list(map(get_edges_from_tri,tri_save)))
    edges0s = edges_save[:-1]
    edges1s = edges_save[1:]

    def _calculate_t1s_frame_to_frame_by_ctype(edges0,edges1):
        return calculate_t1s_frame_to_frame_by_ctype(edges0, edges1, c_types)

    t1s_by_frame = np.array([_calculate_t1s_frame_to_frame_by_ctype(edges0,edges1) for (edges0,edges1) in zip(edges0s,edges1s)])
    return t1s_by_frame


def get_AVE_x(x, c_types, c_type=0):
    """
    Find the positions of AVE/DVE cells by subsetting the matrix of cell centres
    """
    mask = c_types == c_type
    AVE_x = x[:, mask]
    return AVE_x

def get_average_centroid(x_sample):
    """
    Calculate the average centroid of a subset of the matrix of cell centres; x_sample.
    """
    return x_sample.mean(axis=1)

def get_scalar_from_vector(v):
    """
    Find the magnitude of a vector, where the last two dimensions are the x and y coordinate of the vector.
    """
    return np.sqrt(v[..., 0] ** 2 + v[..., 1] ** 2)



def run_analysis(h5_file,run_options,L):
    """
    h5_file: path to hdf5 file
    run_options: dictionary detailing run options of the simulation. Can leave as {}
    L: float value defining size of domain.
    """

    sim_dict = load_hdf5_skeleton(h5_file)

    x = np.array(sim_dict["x_save"],dtype=np.float32)
    tri = np.array(sim_dict["tri_save"],dtype=np.int32)
    c_types = np.array(sim_dict["c_types"],dtype=np.int32)
    meshes = geo.mesh_assembler(x, tri, L, run_options)

    AVE_x = get_AVE_x(x, c_types)
    av_AVE_x = get_average_centroid(AVE_x)
    av_AVE_d = get_scalar_from_vector(av_AVE_x - av_AVE_x[0])

    def get_ave_connected_components(mesh):
        """Calculate connected components of the triangulation by cell type"""
        return top.count_connected_components(mesh.tri, c_types, len(mesh.x))[0]

    ave_connected_components = np.array(list(map(get_ave_connected_components, meshes)))

    t1s_across_sim = calculate_t1s_across_sim_by_ctype(tri,c_types)

    return av_AVE_d, ave_connected_components,t1s_across_sim


