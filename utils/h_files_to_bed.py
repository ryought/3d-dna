#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

def parse_cprops(cprops_fname):
    names = {}
    lengths = {}
    with open(cprops_fname) as f:
        for line in f:
            record = line.rstrip('\n').split(' ')
            names[int(record[1])] = record[0]
            lengths[int(record[1])] = int(record[2])
    return names, lengths

def parse_asm(asm_fname):
    scafs = {}
    with open(asm_fname) as f:
        for i, line in enumerate(f):
            record = line.rstrip('\n').split(' ')
            scaf = [(
                int(r[1:]) if r[0] == '-' else int(r),
                True if r[0] == '-' else False,
                ) for r in record]
            scafs[i+1] = scaf
    return scafs

def calc_start_pos(lengths, scafs):
    start_pos = {}
    current_pos = 0
    for sid in scafs.keys():
        for cid, cori in scafs[sid]:
            start_pos[cid] = current_pos
            current_pos += lengths[cid]
    return start_pos

if __name__ == '__main__':
    basedir = '/work/ryought/hi-c-assembly/human/3ddna-prob-100000-new-more-iteration/'
    names, lengths = parse_cprops(basedir + 'GSE95797_Hs1.3.cprops')
    scafs_final    = parse_asm(basedir + 'GSE95797_Hs1.3.asm')
    start_poss = calc_start_pos(lengths, scafs_final)
    scale = 2

    iterations = [
            (0, '255,255,0'),
            (1, '0,0,255'),
            (2, '0,255,0'),
            (3, '0,255,255'),
            (4, '255,0,255'),
            (5, '0,0,100'),
            (6, '0,100,0'),
            (7, '0,100,100'),
            (8, '100,0,100'),
            (9, '100,100,0'),
            # (10, '100,100,100'),
            ]
    N = 10
    iterations = [
            (i, '0,{},{}'.format(int(255 * (i/(N-1))), int(255 * (1 - i/(N-1)))))
            for i in range(N)
            ]
    print(iterations)

    for iteration, color in iterations:
        bed = 'chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment\n'
        scafs = parse_asm(basedir + 'hfiles_3/h.scaffolds.original.notation.step.{}.txt'.format(iteration))
        for i, sid in zip(range(1000000), scafs.keys()):
            # begin pos
            head_id, _ = scafs[sid][0]
            tail_id, _ = scafs[sid][-1]
            head_pos = start_poss[head_id]
            tail_pos = start_poss[tail_id]
            if head_pos < tail_pos:
                start_pos = head_pos
                end_pos = tail_pos + lengths[tail_id]
            else:
                start_pos = tail_pos
                end_pos = head_pos + lengths[head_id]

            bed += 'assembly\t{start}\t{end}\tassembly\t{start}\t{end}\t{color}\t{comment}\n'.format(
                    start=start_pos // scale,
                    end=end_pos // scale,
                    color=color,
                    comment=str(iteration)+':'+str(sid)+':'+','.join([str(cid) for cid, _ in scafs[sid]]) if iteration == 0 else str(iteration)+':'+str(sid)
                    )
        with open(basedir + 'notation.step.{}.bed'.format(iteration), mode='w') as f:
            f.write(bed)

