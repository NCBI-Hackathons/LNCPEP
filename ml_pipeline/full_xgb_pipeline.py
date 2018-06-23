import csv
from Bio import SeqIO
import pickle
import numpy as np

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

from ingest_data import *


'''
Three steps:

1.
Use bedtools intersect coordinates and the chicken reference genome
to create a csv file with the raw sequences corresponding to the RNA transcripts

2.
break the csv rows into subsequences and featurize the subsequences into kmer frequency vectors
seperate these into training/test on a per transcript basis

3.
train an XGBClassifier
'''

# Filenames for step 1
intersect_output_fn = 'data/bed_intersect_1/srr924202_new_screen_intersect.txt'
reference_genome_fn = 'data/Gallus_gallus-5.0_genomic.fasta'
raw_seq_csv_fn = 'new_filtered_bt_intersect_ML_input.csv'

# Params and filenames for step 2
subseq_length = 100
subseq_window_stride = 100
ks = [2,3] # k_mer sizes for featurization
split_denom = 5
# One in every [split_denom] RNA transcripts will be used to create the test set
pickled_featurized_subsequences_fn_suffix = '_featurized_subseqs_1.p'
# train and test will be prepended to this

# Filename for step 3
model_save_fn = 'xgb_pTUPC_1.p'



print('Step 1:')
print('Collecting raw sequences from bed intersect coordinates')
print('And writing to csv:', raw_seq_csv_fn)

with open(intersect_output_fn, 'r') as f:

    intersects = f.readlines()

intersects = [x.split() for x in intersects]

print('Number of lncRNA transcripts:', len(intersects))

hit = 0
no_hit = 0

for i in intersects:

    if i[-1] == '0':

        no_hit += 1

    else:

        hit += 1

print('pTUPCs (>0 intersects with peptide dataset):', hit)
print('non-TUPC lncRNAs (0 intersects with peptide dataset):', no_hit)

linc_collector = []

gallus_reference_genome = list(SeqIO.parse(reference_genome_fn, 'fasta'))

print(len(gallus_reference_genome))
print(gallus_reference_genome[0].id)

ref_genome_dict = {}

for trans in gallus_reference_genome:

    ref_genome_dict[trans.id] = trans

with open(raw_seq_csv_fn, "w") as csv_file:

    writer = csv.writer(csv_file, delimiter=',')

    for i in intersects[:]:

        ind = intersects.index(i)

        if ind % 200 == 0:

            print(ind)



        id = create_id_from_row(i, [0,1,2,3,5])

        transcript = i[0]
        start = int(i[1])
        stop = int(i[2])

        #print(transcript,start,stop)

        full_transcript_seq = ref_genome_dict[transcript].seq

        subseq = full_transcript_seq[start:stop]

        #print(len(subseq))

        #seq = ''

        pTUPC = int(bool(int(i[-1])))

        row = [id, str(subseq).upper(), int(i[-1]), pTUPC]

        #print(row)

        #linc_collector.append(row)
        writer.writerow(row)

print("")
print('Step 2:')
print('Create kmer frequency vectors from subsequences')
print('k values:', ks)
print('split into train and test sets')
print('test fraction:', float(1)/split_denom)

with open(raw_seq_csv_fn, 'r') as f:

    sequences = f.readlines()

vec_length = [4 ** x for x in ks]
vec_length = sum(vec_length)
print('Feature vector length:', vec_length)

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

                if ind % split_denom == 0:

                    test_lRNA_collector.append(z)

                else:

                    train_lRNA_collector.append(z)

        except:

            pass

#print(lRNA_collector)
train_a = np.array(train_lRNA_collector)
test_a = np.array(test_lRNA_collector)
print('Train shape:', train_a.shape)
print('Test shape:', test_a.shape)

with open('train' + pickled_featurized_subsequences_fn_suffix, 'wb') as f:

    pickle.dump(train_a, f)

with open('test' + pickled_featurized_subsequences_fn_suffix, 'wb') as f:

    pickle.dump(test_a, f)

print("")
print('Step 3:')
print('Train an XGBClassifier and apply it to the test set')
print('Pickle the model for future use:', model_save_fn)

X_train = train_a[:,:-1]
y_train = train_a[:,-1:]

X_test = test_a[:,:-1]
y_test = test_a[:,-1:]

print('Fitting model')
model = XGBClassifier()
model.fit(X_train, y_train)

print('making predictions')
y_pred = model.predict(X_test)
predictions = [round(value) for value in y_pred]

accuracy = accuracy_score(y_test, predictions)
print("Accuracy on test set: %.2f%%" % (accuracy * 100.0))

with open(model_save_fn, 'wb') as f:

    pickle.dump(model, f)
