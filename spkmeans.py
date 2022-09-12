import myspk as mkm
import numpy as np
import pandas as pd
import sys
from numpy import linalg as LA


def exit_func (code):
    if code == 1:
        print("Invalid Input!")
    if code == 0:
        print("An Error Has Occurred")
    exit(1)


def is_valid_filename():
    try:
        f = open(sys.argv[3], 'r')
        f.close()
    except:
        exit_func(1)


def is_valid_func_name():
    func_name=sys.argv[2]
    if func_name != "wam" and func_name != "ddg" and func_name != "lnorm" and func_name != "jacobi" and func_name != "spk":
        exit_func(1)


def check_legal():
    if len(sys.argv) != 4:  # there are 5 or 6 arguments depending on max iter
        exit_func(1)
    if not sys.argv[1].isdigit():
        exit_func(1)
    is_valid_filename()
    is_valid_func_name()


def is_cluster_legal(points_list):
    num_of_clusters = int(sys.argv[1])
    if num_of_clusters < 0:
        exit_func(1)


def read_input(file_name1):
    file1 = pd.read_csv(file_name1, sep=",", header=None)
    vectors = file1.to_numpy()
    vectors = vectors.tolist()
    return vectors


def probabilities(dls, sum_dls):
    return dls/sum_dls


def Kmeans_pp(vectors, k):  # k is number of clusters
    #print("clusters: ",k)
    #print("vectors_list: \n",vectors)
    vectors = np.array([np.array(vector) for vector in vectors])
    #print("vectors numpy: \n",vectors)
    n = vectors.shape[0]  # this is number of vectors
    dim = vectors.shape[1]
    rand_cent_list=[]
    i = 1  # number of initialized centroids
    np.random.seed(0)
    centroids = np.zeros((k, dim))
    random_centroid = np.random.choice(n)
    rand_cent_list.append(random_centroid)
    random_centroid = vectors[random_centroid]
    centroids[0] = random_centroid
    dls = np.zeros((n, 1))
    while i < k:
        for l in range(n):
            dl = float('inf')
            for j in range(i):
                dl = min(dl, np.linalg.norm(vectors[l] - centroids[j])**2)
            dls[l] = dl
        probs = probabilities(dls, np.sum(dls))
        probs = probs.flatten()
        random_centroid = np.random.choice(n, p=probs)
        rand_cent_list.append(random_centroid)
        centroids[i] = vectors[random_centroid]
        i += 1
    centroids = centroids.tolist()
    vectors = vectors.tolist()
    centroids = mkm.fit(vectors, centroids, k, dim, n, 3, 0)
    chosen_str = ""
    for i in range(len(rand_cent_list)-1):
        chosen_str = chosen_str + str(rand_cent_list[i]) + ","
    chosen_str = chosen_str + str(rand_cent_list[len(rand_cent_list) - 1])
    print(chosen_str)
    print_matrix(centroids)


def print_matrix(matrix):
    for row in matrix:
        print(','.join([f'{num:.4f}' for num in row]))


def master(vectors):
    k = int(sys.argv[1])
    goal = sys.argv[2]
    if goal == "wam":
        matrix = mkm.mainPy(vectors, len(vectors), len(vectors[0]), 0, k)
        print_matrix(matrix)
    elif goal == "ddg":
        matrix = mkm.mainPy(vectors, len(vectors), len(vectors[0]), 1, k)
        print_matrix(matrix)
    elif goal == "lnorm":
        matrix = mkm.mainPy(vectors, len(vectors), len(vectors[0]), 2, k)
        w, v = LA.eigh(matrix)
        print_matrix(matrix)
    elif goal == "jacobi":
        matrix = mkm.mainPy(vectors, len(vectors), len(vectors[0]), 3, k)
        print_matrix(matrix)
    elif goal == "spk":
        # print('before mkm')
        matrix = mkm.mainPy(vectors, len(vectors), len(vectors[0]), 4, k)
        # print('matrix after jacobi before kmeans\n')
        # print_matrix(matrix)
        # print('end of matrix \n')
        # print('b4 kmeans #########################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%^^^^^^^^^^^^^^^^^^^^^^^^^^&&&&&&&&&&&&&&&&&&&&&&&&&&&&&********************')
        Kmeans_pp(matrix, len(matrix[0]))


if __name__ == '__main__':
    check_legal()
    vectors = read_input(sys.argv[3])
    master(vectors)
