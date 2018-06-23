from ingest_data import *
import pickle
import numpy as np

sequences_input_fn = 'new_filtered_bt_intersect_ML_input.csv'

with open(sequences_input_fn, 'r') as f:

    sequences = f.readlines()

subseq_length = 100
subseq_window_stride = 20

ks = [2,3]

vec_length = [4 ** x for x in ks]
vec_length = sum(vec_length)
print(vec_length)

train_lRNA_collector = []
test_lRNA_collector = []

for s in sequences[:]:

    ind = sequences.index(s)

    if ind % 200 == 0:

        print(ind)

    row = s.split(',')

    for i in range(0,len(row[1]), subseq_window_stride):

        try:
            subseq = row[1][i : i + subseq_length]

            #print(i, k_mer_frequencies(subseq, ks, include_missing=True, vector=True))

            z = list(k_mer_frequencies(subseq, ks, include_missing=True, vector=True))
            z.append(int(row[-1]))

            if len(z) == vec_length + 1:

                if ind % 5 == 0:

                    test_lRNA_collector.append(z)

                else:

                    train_lRNA_collector.append(z)

        except:

            pass

#print(lRNA_collector)
train_a = np.array(train_lRNA_collector)
test_a = np.array(test_lRNA_collector)
print(train_a.shape)
print(test_a.shape)

with open('train_featurized_5.p', 'wb') as f:

    pickle.dump(train_a, f)

with open('test_featurized_5.p', 'wb') as f:

    pickle.dump(test_a, f)
