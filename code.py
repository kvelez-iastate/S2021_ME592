import sys
sys.path.append(r'C:\Users\ks71275\OneDrive-Deere&Co\OneDrive - Deere & Co\Documents\ME592_HW2\tsne_python\tsne_python')
# import tsne
import numpy as np
import os

#%%

def read_result_file(filepath):
    
    # read the data file
    data = np.genfromtxt(filepath, skip_header = 1)
    
    # split the data file into the 3 surfaces
    # transpose so that rows correspond to x,y,z coordinates
    surf_1 = data[:204,:].T
    surf_2 = data[204:408,:].T
    surf_3 = data[408:,:].T
    
    # create a single list for all 3 surfaces
    surfaces = [surf_1, surf_2, surf_3]
    
    # return the surface x,y,z data
    return surfaces

def read_input_file(path):
    
    # read the data files
    # transpose so that rows correspond to x,y,z coordinates
    surf_1 = np.genfromtxt(path + '\\smesh.1.1.txt', skip_header=5, skip_footer=1, usecols=(0,1,2)).T
    surf_2 = np.genfromtxt(path + '\\smesh.1.2.txt', skip_header=5, skip_footer=1, usecols=(0,1,2)).T
    surf_3 = np.genfromtxt(path + '\\smesh.1.3.txt', skip_header=5, skip_footer=1, usecols=(0,1,2)).T
    
    # create a single list for all 3 surfaces
    surfaces = [surf_1, surf_2, surf_3]
    
    # return the surface x,y,z data
    return surfaces

def read_result_filename(filename,temperatures,pressures):
    
    # split the filename using '_' as a delimiter
    filename_split = filename.split('_')
    
    # extract the desired parameters from the split filename
    timestep = int(filename_split[1])
    temperature_index = int(filename_split[2][1])
    pressure_index = int(filename_split[3][1])
    geom = int(filename_split[4])
    
    # convert temperature index to actual temperature
    # subtract 1 because python is zero-indexed
    temperature = temperatures[temperature_index-1]
    
    # convert pressure index to actual pressure
    # subtract 1 because python is zero-indexed
    pressure = pressures[pressure_index-1]
    
    # return the desired data
    return [timestep, temperature, pressure, geom]

def whiten_surf(a):
    
    # zero-center the x,y,z data
    means = np.mean(a.T,axis=0)
    a_centered = (a.T - means).T
    
    # calculate the covariance matrix
    cov = np.cov(a_centered)
    
    # calculate the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    
    # calculate the inverse sqrt of eigenvalues
    diag_eigenvalues = np.diag(1/(eigenvalues**0.5))
    
    # perform whitening
    wpca = np.dot(np.dot(diag_eigenvalues, eigenvectors.T), a_centered)
    
    # return the whitened data
    return wpca

def list_folder_contents(path):
    
    # initialize an empty list
    entries = []
    
    # loop through the folder to retieve its contents
    for entry in os.listdir(path):
        entries.append(entry)
        
    # return the folder contents
    return entries
    
    

###############################################################################
# Part 1
###############################################################################

# user-defined inputs
input_path = r'C:\Users\ks71275\OneDrive-Deere&Co\OneDrive - Deere & Co\Documents\ME592_HW2\DM\data\Input_geometry'
result_path = r'C:\Users\ks71275\OneDrive-Deere&Co\OneDrive - Deere & Co\Documents\ME592_HW2\DM\data\final_geometry'
temperatures = [300, 350, 400, 450, 500]
pressures = [76, 80, 84]

# create a list of result file names
result_file_names = list_folder_contents(result_path)

# initialize an empty list of ordered pairs for the 80 and 140 timesteps
ordered_pairs_80 = []
ordered_pairs_140 = []

# initialize a counter varible so progress can be printed
n = 0

# create an ordered pair for every results file
for filename in result_file_names:
    
    # extract the required parameters from the file name
    timestep, temperature, pressure, geom = read_result_filename(filename, temperatures, pressures)
    
    # load the x,y,z data for all 3 surfaces in the results file
    result_surfs = read_result_file(result_path + '\\' + filename)

    # load the x,y,z data for all 3 surfaces in the corresponding input file
    input_surfs = read_input_file(input_path + '\\run' + str(geom))
    
    # create an ordered pair per the assignment description
    ordered_pair = [[input_surfs, temperature, pressure], result_surfs]
    
    # append to the appropriate ordered pair list depending on the timestep
    if timestep == 80:
        ordered_pairs_80.append(ordered_pair)
    elif timestep == 140:
        ordered_pairs_140.append(ordered_pair)
    
    # print the number of files that have been read (of 1294 files)
    n=n+1
    print(n)
    
###############################################################################
# Part 2
###############################################################################

#%%
# create a list of the ordered pairs at each time step that can be iterated on
ordered_pairs = [ordered_pairs_80, ordered_pairs_140]

# initialize counter varibles
n = 0
n_surfaces_per_geom = 3

# for each ordered pair list
for i, ordered_pair_list in enumerate(ordered_pairs):
    
    # for each ordered pair
    for j, ordered_pair in enumerate(ordered_pair_list):
    
        # for each of  3 input surfaces
        for k, surf in enumerate(ordered_pair[0][0]):
            
            # whiten the surface
            ordered_pairs[i][j][0][0][k] = whiten_surf(surf)
        
        # print the number of files that have been read (of 1294 files)
        n=n+1
        print(n)
        
###############################################################################
# Part 3
###############################################################################
#%%


def create_vector_row(ordered_pair):
    
    # pull the required data from the ordered pair
    result_surf_1 = ordered_pair[1][0]
    result_surf_2 = ordered_pair[1][1]
    result_surf_3 = ordered_pair[1][2]
    temperature = np.ones((3,1))*ordered_pair[0][1]
    pressure = np.ones((3,1))*ordered_pair[0][1]

    # initialize an empty row array
    row = []

    # create the row array
    row = np.append(result_surf_1, result_surf_2, axis=1)
    row = np.append(row, result_surf_3, axis = 1)
    row = np.append(row, temperature, axis = 1)
    row = np.append(row, pressure, axis = 1)
    
    # return the row array
    return row

# initialize an empty vector
vector = []

for i, ordered_pair in ordered_pairs_80:
    
    # create the row array
    row = create_vector_row(ordered_pair)
    
    # append the row array
    vector.append(row)
