import pandas as pd
import numpy as np
import sys

def read_input(file_name1):
    file1 = pd.read_csv(file_name1,sep=",",header=None)
    vectors = file1.to_numpy()
    vectors=vectors.tolist()
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


#v=read_input(sys.argv[1])
v=[[-5.056, 11.011],
   [-6.409, -7.962],
   [5.694, 9.606],
   [6.606, 9.396],
   [-6.772, -5.727],
   [-4.498, 8.399],
   [-4.985, 9.076],
   [4.424, 8.819],
   [-7.595, -7.211]]

x=(wamcheck(v))
print_matrix(x)
