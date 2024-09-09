from collections import defaultdict, Counter

import argparse

import file_io


def main():
    # CL arguments
    args   = handle_args()
    in_fp  = args['data']
    out_fp = args['output']
    v_lim  = args['variant_threshold']

    reads = file_io.parse_bowtie(in_fp)

    # Count variants in all sites
    sequence_counters = defaultdict(Counter)
    coverages = defaultdict(int)
    sites_to_reads = defaultdict(list)

    for idx, read in enumerate(reads):
        if read['regions'] > 1:
            continue

        for i, n in enumerate(read['sequence'], start=read['start']):
            sequence_counters[i].update([n] * read['clones'])
            coverages[i] += read['clones']
            sites_to_reads[i].append(idx)

    # Find credible variants
    variant_sites = {}
    for key, val in sequence_counters.items():
        variants  = val.most_common()
        threshold = variants[0][1] * v_lim
        variants  = [x for x in variants if x[1] > threshold]

        if len(variants) == 1:
            continue
        if len(variants) > 2:
            raise RuntimeError

        variant_sites[key] = [x[0] for x in variants]

    # TODO: linkage


def handle_args():
    parser = argparse.ArgumentParser(description='HLA phase deconvolvinator')
    parser.add_argument('-d', '--data',
        help='Bowtie alignments', required=True, type=str)
    parser.add_argument('-o', '--output',
        help='Output filepath', required=True, type=str)
    parser.add_argument('-v', '--variant-threshold',
        help='Fraction of most common variant to use as threshold for a credible variant', default=.1, type=float)

    return vars(parser.parse_args())


if __name__ == "__main__":
    main()
