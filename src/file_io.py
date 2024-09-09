def parse_bowtie(fp):
    '''
    Returns a list object where each element is a dictionary of:

    FASTA sequence.
    Number of identical reads.
    Regions it mapped to.
    List of scores for each mapped region.
    Start position (0-indexed).
    '''

    ret = []
    with open(fp) as in_f:
        for info_line in in_f:
            if not info_line.startswith('>'):
                continue

            toks     = info_line[1:].split(':')
            scores   = [int(x.split('x')[1]) for x in toks[2].split('_')[1:]]
            seq_line = next(in_f)

            ret.append({'sequence': seq_line.strip().strip('~'), 'clones': int(toks[0].split('_')[1]),
                'regions': int(toks[1].split('_')[1]), 'scores': scores, 'start': int(toks[3].split('_')[1]) - 1})

    return ret
