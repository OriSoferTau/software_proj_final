import pandas as pd
import numpy as np
import sys

def read_input(file_name1):
    file1 = pd.read_csv(file_name1,sep=",",header=None)
    vectors = file1.to_numpy()
    print(vectors)
    #vectors=vectors.tolist()
    return vectors





def wamcheck(vectors):
    lst=[]
    for i in range(len(vectors)):
        lst_vec=[]
        for j in range(len(vectors)):
            if i==j:
                lst_vec.append(0)
                continue
            norm=np.sqrt((vectors[i][0]-vectors[j][0])**2+((vectors[i][1]-vectors[j][1])**2))
            norm=norm/2
            norm=np.exp(-norm)
            lst_vec.append(norm)
        lst.append(lst_vec)
    return lst
def print_matrix(matrix):
    cent_str=""
    for i in range(len(matrix)):
        for j in range(len(matrix[i])-1):
            temp=str(float("{:.3f}".format((matrix[i][j]))))
            cent_str=cent_str+temp+","
        temp=str(float("{:.3f}".format((matrix[i][len(matrix[i])-1]))))
        cent_str=cent_str+str(temp)+"\n"
    print(cent_str)


v=read_input(sys.argv[1])
w,v=np.linalg.eigh(v)
print(w)
#print(w)


