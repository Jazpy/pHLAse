from collections import defaultdict, Counter

import argparse

import file_io
from allele_dag import ADAG


def main():
    # CL arguments
    args   = handle_args()
    in_fp  = args['data']
    out_fp = args['output']
    v_lim  = args['variant_threshold']
    should_draw = args['draw']

    reads = file_io.parse_bowtie(in_fp)

    # Count variants in all sites
    sequence_counters = defaultdict(Counter)
    coverages = defaultdict(int)
    sites_to_reads = defaultdict(set)

    for idx, read in enumerate(reads):
        if read['regions'] > 1:
            continue

        for i, n in enumerate(read['sequence'], start=read['start']):
            sequence_counters[i].update([n] * read['clones'])
            coverages[i] += read['clones']
            sites_to_reads[i].add(idx)

    # Find credible variants
    variant_sites = []
    variant_site_info = {}
    for key, val in sequence_counters.items():
        variants  = val.most_common()
        threshold = variants[0][1] * v_lim
        variants  = [x for x in variants if x[1] > threshold]

        if len(variants) == 1:
            continue
        if len(variants) > 2:
            raise RuntimeError

        variant_sites.append(key)
        variant_site_info[key] = [x[0] for x in variants]

    # Linkage
    dag = ADAG()
    dag.add_unlinked_level(variant_site_info[variant_sites[0]])
    for i in range(len(variant_sites) - 1):
        s0 = variant_sites[i]
        s1 = variant_sites[i + 1]

        # Reads that cover both sites
        joint_reads = sites_to_reads[s0].intersection(sites_to_reads[s1])

        if not joint_reads:
            dag.add_unlinked_level(variant_site_info[s1])
            continue

        allele_linkages = []
        for read in joint_reads:
            r  = reads[read]
            a0 = r['sequence'][s0 - r['start']]
            a1 = r['sequence'][s1 - r['start']]

            allele_linkages.extend([(a0, a1)] * r['clones'])

        top_two = Counter(allele_linkages).most_common(2)
        dag.add_level(top_two)

    # Data output
    paths = dag.get_paths()

    seq_len = max(sequence_counters.keys())
    ref_seq = []
    for site in range(seq_len):
        if site in sequence_counters and site not in set(variant_sites):
            ref_seq.append(sequence_counters[site].most_common(1)[0][0])
        elif site in sequence_counters and site in set(variant_sites):
            ref_seq.append('V')
        else:
            ref_seq.append('X')

    with open(out_fp, 'w') as out_f:
        for path in paths:
            variant_seq = ref_seq.copy()
            for site, a in zip(variant_sites, path):
                variant_seq[site] = a

            out_f.write(f'{"".join(variant_seq)}\n')

    # Optional visualization
    if should_draw:
        dag.draw()


def handle_args():
    parser = argparse.ArgumentParser(description='HLA phase deconvolvinator')
    parser.add_argument('-d', '--data',
        help='Bowtie alignments', required=True, type=str)
    parser.add_argument('-o', '--output',
        help='Output filepath', required=True, type=str)
    parser.add_argument('-v', '--variant-threshold',
        help='Fraction of most common variant to use as threshold for a credible variant', default=.1, type=float)
    parser.add_argument('--draw',
        help='Draw allele DAG', action=argparse.BooleanOptionalAction)

    return vars(parser.parse_args())


if __name__ == "__main__":
    main()
