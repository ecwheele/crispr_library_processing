from collections import Counter
import numpy as np

after_guide = 'gttttagagctaGAAAtagcaagttaaaataa'.upper()
before_guide = 'TTCTTGTGGAAAGGACGAAACACCG'



def hamming_distance(s1, s2):
    """
    calculate hamming distance between two strings (must be same length)
    :param s1: first string
    :param s2: second string
    :return: number of mismatches between strings (if not same length, returns 100)
    """
    if len(s1) != len(s2):
        return 100
    else:
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def count_guides(sample_fastq, front_seq=before_guide, back_seq=after_guide):
    """

    :param sample_fastq: demuxed fastq file to count guides
    :param front_seq: sequence before guide
    :param back_seq: sequence after guide
    :return: Counter of everything detected in between
    """

    gRNA_counts = Counter()
    read_count = 0
    good_read_count = 0

    with open(sample_fastq, 'r') as fastq_file:
        while True:
            try:
                fastq_file.next()
                seq = fastq_file.next()
                fastq_file.next()
                fastq_file.next()

                read_count +=1

                start_position = seq.find(front_seq)+len(front_seq)
                end_position = seq.find(back_seq)

                if (start_position != (-1+len(front_seq))) & (end_position != -1):

                    gRNA = seq[start_position: end_position]
                    gRNA_counts[gRNA] += 1
                    good_read_count +=1

                else:

                    gRNA_counts['None'] +=1

            except StopIteration:
                break

    return gRNA_counts, read_count, good_read_count


def assign_guides(gRNA_counts):

    cutoff = np.mean(gRNA_counts.values())*10

    to_keep = dict()

    for item in gRNA_counts.keys():
        count = gRNA_counts[item]
        if count > cutoff:
            to_keep[item] = count

    return to_keep


def make_summary_df(read_count, good_read_count, to_keep, demuxed_fastq):
    name = os.path.basename(demuxed_fastq).split("_")[0]
    index = os.path.basename(demuxed_fastq).split("_")[1].rstrip(".fastq")
    rows = [name, index, 'total_reads', read_count]
    rows.append([name, index, 'reads_with_guide', good_read_count])


    if len(to_keep) > 0:
        for item in to_keep.keys():
            row = [name, index, item, to_keep[item]]
            rows.append(row)

    df = pd.DataFrame(rows, orient='index')

    return df


def process_one_sample(sam)






