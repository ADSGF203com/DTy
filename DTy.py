import numpy as np, pandas as pd, re, click, os #, filelock
import hashlib
from _collections import OrderedDict

try :
    from .configure import dbs
    from .uberBlast import uberBlast, readFastq, rc
except :
    from configure import dbs
    from uberBlast import uberBlast, readFastq, rc


#logging.basicConfig(format='%(asctime)s - %(message)s', level=20)

refseqs = None
core_genes = None
repr = None
hiercc_repr = None
species = None
key_level = None
scheme_info = None



def get_md5(value, dtype=str):
    m = hashlib.md5(str(value).encode()).hexdigest()
    if dtype == str:
        return m
    else:
        return int(m, 16)


@click.command()
@click.option('-q', '--query', help='input assembly. fasta or fastq format, can be gziped.')
@click.option('-d', '--dbname', help='scheme that defined in configure.py. default: Neisseria_dcgMLST', default='Neisseria_dcgMLST')
@click.option('-o', '--outfile', help='alleles output file')
@click.option('-n', '--n_thread', help='n_threads [default:8]', default=8, type=int)
def main(query, dbname, outfile, n_thread) :
    db = dbs[dbname]
    
    global refseqs, core_genes, repr, hiercc_repr, species, key_level, scheme_info
    
    refseqs = db['refseqs']
    core_genes = db['core_genes']
    repr = db['repr']
    hiercc_repr = db['hiercc_repr']
    species = db['species']
    key_level = db['levels']
    scheme_info = db['scheme_info']

    
    alleles, hiercc = cgmlst(query, key_level, n_thread=n_thread)
    with open(outfile, 'wt') as fout :
        for gene, allele in sorted(alleles.items()) :
            if 'sequence' in allele :
                fout.write('>{gene_name} value_md5={value_md5} CIGAR={CIGAR} accepted={flag} reference={reference} identity={identity:.3f} coordinates={coordinates}\n{sequence}\n'.format(**allele))
    print('\t'.join([query] + ['{0}={1}'.format(k, hiercc[k]) for k in [key_level, 'species', 'reference', 'allelic_distance'] ]))
    return outfile


def cgmlst(query, key_level, n_thread=8) :
    levels = [int(hc) for hc in re.findall('\d+', key_level)]
    alleles = nomenclature(query, n_thread)

    # lock = filelock.FileLock(repr + '.lock')
    # with lock :
    repr_profile = pd.read_parquet(repr)
    repr_profile = repr_profile.set_index(repr_profile.columns[0])
    repr_hiercc = pd.read_csv(hiercc_repr, sep='\t')
    repr_hiercc = repr_hiercc.set_index(repr_hiercc.columns[0])
    hc_species = pd.read_csv(species, sep='\t', dtype=str)
    
    profile = np.array([('-' if v.startswith('-') else v) for v in [v.get('value_md5', '-') for k, v in sorted(alleles.items())]])
    relshare = (np.sum((repr_profile.values == profile) & (profile != ''), 1)*profile.size+0.1)/(np.max(( np.sum((repr_profile.values != '') & (profile != '-'), 1), \
                                                                                np.sum((repr_profile.values != ''), 1)*0.97 ), 0)+0.1)
    max_idx = np.argmax(relshare)
    min_dist = int(profile.size-relshare[max_idx]+0.5)
    ref_repr = repr_profile.index[max_idx]
    clusters = repr_hiercc.loc[ref_repr].values.astype(str).tolist()
        
    clusters[:min_dist] = ['N'] * min_dist
    hiercc = {'reference':ref_repr, 'allelic_distance':min_dist}
    hiercc[key_level] = '.'.join([clusters[d] for d in levels])

    sp = 'unknown species'
    for hc, tax in hc_species.values :
        if hiercc[key_level].startswith(hc+'.') :
            sp = tax
            break

    hiercc['species'] = sp
    return alleles, hiercc
    
    
def nomenclature(query, n_thread=8) :
    core = {}
    with open(core_genes, 'rt') as fin :
        for line in fin :
            core[line.strip().split()[0]] = 1
            
    blastab = uberBlast(
        '-r {0} -q {1} -f --blastn --diamond --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 21,21 -m --merge_gap 300'.format(
            query, refseqs, scheme_info['min_iden']-0.05, 0.05, n_thread).split())
    merges = {}
    for b in blastab.T[16] :
        if len(b) > 4 :
            key = tuple(b[3:])
            merges[key] = b[:3]
    bsn = { b[15]:b[:15] for b in blastab }
    for ids, (score, iden, size) in merges.items() :
        bs = np.array([ bsn[i] for i in ids ], dtype=object)
        if np.unique(bs.T[1]).size == 1 :
            bs[0][2], bs[0][3], bs[0][7], bs[0][9], bs[0][11], bs[0][14] = iden, size, bs[-1][7], bs[-1][9], score, 'COMPLEX'
        for i in ids[1:] :
            bsn.pop(i)
        bsn[ids[0]] = bs[0]
    blastab = np.array(list(bsn.values()))

    blastab = blastab[(blastab.T[2] >= scheme_info['min_iden']) & ((blastab.T[7]-blastab.T[6]+1) >= scheme_info['min_frag'] * blastab.T[12])]
    blastab.T[11] = blastab.T[11]*(blastab.T[7] - blastab.T[6]+1)/blastab.T[12]
    blastab = blastab[np.lexsort((blastab[:, 8], blastab[:, 1]))]
    for i, b0 in enumerate(blastab[:-1]) :
        if b0[0] == '' :
            continue
        s0, e0 = sorted(b0[8:10])
        todel = []
        for b1 in blastab[i+1:] :
            s1, e1 = sorted(b1[8:10])
            if b0[1] != b1[1] or e0 < s1 :
                break
            ovl = min(e0, e1) - max(s0, s1) + 1
            if ovl >= 0.5 * (e0-s0+1) or ovl >= 0.5 * (e1-s1+1) :
                sc0, sc1 = abs(b0[11]), abs(b1[11])
                g0, g1 = b0[0].rsplit('_', 1)[0], b1[0].rsplit('_', 1)[0]
                if b0[2] < b1[2]*scheme_info['max_iden'] or (b1[2] >= b0[2]*scheme_info['max_iden']
                    and (sc0 < sc1 or (sc0 == sc1 and b0[0] > b1[0]))) :
                    b0[11] = -sc0
                    if g0 == g1 or sc0 < sc1 * scheme_info['max_iden'] or b0[2] < b1[2]*scheme_info['max_iden'] :
                        b0[0] = ''
                        break
                else :
                    b1[11] = -sc1
                    if g0 == g1 or sc1 < sc0 * scheme_info['max_iden'] or b1[2] < b0[2]*scheme_info['max_iden'] :
                        todel.append(b1)
        if b0[0] and len(todel) :
            for b1 in todel :
                b1[0] = ''
    blastab = blastab[blastab.T[0] != '']
    blastab = blastab[np.lexsort([-blastab.T[11], [b.rsplit('_', 1)[0] for b in blastab.T[0]]])]
    alleles = OrderedDict()
    for bsn in blastab :
        gene = bsn[0].rsplit('_', 1)[0]
        if gene in alleles :
            if alleles[gene]['score']*scheme_info['max_iden'] > bsn[11] :
                continue
            alleles[gene]['coordinates'].append((bsn[1], bsn[8], bsn[9]))
            alleles[gene]['flag'] |= 32
            continue
        flag = 0
        if bsn[6] > 1 or bsn[7] < bsn[12] :
            flag = 64
            if bsn[14] != 'COMPLEX' :
                if bsn[6] > 1 :
                    bsn[14] = '{0}D{1}'.format(bsn[6]-1, bsn[14])
                if bsn[7] < bsn[12] :
                    bsn[14] = '{0}{1}D'.format(bsn[14], bsn[12]-bsn[7])
        alleles[gene] = {'gene_name': gene, 'CIGAR':bsn[0]+':'+bsn[14], 'reference':os.path.basename(query),
                         'identity': bsn[2], 'coordinates':[(bsn[1], bsn[8], bsn[9])], 'flag':flag, 'score':bsn[11]}

    seq, qual = readFastq(query)
    for gene, allele in sorted(alleles.items()) :
        if allele['flag'] & 96 == 96 :
            alleles.pop(gene)
            continue
        if allele['flag'] & 32 > 0 :
            allele['sequence'] = 'DUPLICATED'
        else :
            c, s, e = allele['coordinates'][0]
            ss = seq[c][s - 1:e] if s < e else rc(seq[c][e - 1:s])
            qs = (min(qual[c][s - 1:e] if s < e else qual[c][e-1:s])) if len(qual) else 0
            
            allele['sequence'] = ss
            if qs < 10 :
                allele['flag'] |= 2
            if allele['flag'] == 0 :
                allele['flag'] = 1
        allele['coordinates'] = ','.join(['{0}:{1}..{2}'.format(*c) for c in  allele['coordinates']])
        allele['value_md5'] = ('' if allele['flag'] < 16 else '-') + get_md5(allele['sequence'])
    return {g:alleles.get(g, {"gene_name":g, "value_md5":"-"}) for g in core }


if __name__ == '__main__':
    main()
