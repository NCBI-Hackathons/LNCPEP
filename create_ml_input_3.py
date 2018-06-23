import csv
from Bio import SeqIO

def create_id_from_row(row, cols_to_concatenate):

    id = ''

    for c in cols_to_concatenate:

        z = row[c]
        z = z.replace('"', '')
        z = z.replace(';', '')

        id += z + '_'

    return id

#intersect_output_fn = 'data/bed_intersect_1/intersect_1.txt'
intersect_output_fn = 'data/bed_intersect_1/srr924202_new_screen_intersect.txt'

with open(intersect_output_fn, 'r') as f:

    intersects = f.readlines()

intersects = [x.split() for x in intersects]

print(len(intersects))

hit = 0
no_hit = 0

for i in intersects:

    if i[-1] == '0':

        no_hit += 1

    else:

        hit += 1

print(hit, no_hit)

linc_collector = []

gallus_reference_genome = list(SeqIO.parse('data/Gallus_gallus-5.0_genomic.fasta', 'fasta'))

print(len(gallus_reference_genome))
print(gallus_reference_genome[0].id)

ref_genome_dict = {}

for trans in gallus_reference_genome:

    ref_genome_dict[trans.id] = trans

with open('new_filtered_bt_intersect_ML_input.csv', "w") as csv_file:

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
