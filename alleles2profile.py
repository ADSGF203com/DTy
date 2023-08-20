import sys, pandas as pd, click
try :
    from .configure import dbs
except :
    from configure import dbs




@click.command()
@click.option('-c', '--scheme', help='scheme name. default: Neisseria_dcgMLST', default='Neisseria_dcgMLST')
@click.option('-p', '--prefix', help='prefix for parquet/tsv output of the allelic profiles', required=True)
@click.argument('fns', nargs=-1)
def main(scheme, prefix, fns) :
    db = dbs[scheme]
    core_genes = db['core_genes']
    genes = {}
    with open(core_genes, 'rt') as fin :
        for line in fin :
            genes[line.strip().split()[0]] = 1

    fields = ['ID'] + sorted(genes.keys())
    data = []
    for fid, fn in enumerate(fns) :
        print(fid, fn)
        tag = fn.rsplit('.', 1)[0]
        profiles = {g:'' for g in genes}
        with open(fn, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    p = line.strip()[1:].split()
                    gene = p[0].split(':', 1)[-1]
                    if gene in genes :
                        id = [ pp.split('=',1)[1] for pp in p if pp.startswith('value_md5=') ][0]
                        if not id.startswith('-') :
                            profiles[gene] = id.replace('-', '').lower()
        data.append([tag] + [profiles[g] for g in fields[1:]])
    
    data = pd.DataFrame(data, columns=fields)
    data.to_parquet(prefix+'.parq', compression='gzip')                    
    data.to_csv(prefix+'.tsv', sep='\t')



if __name__ == '__main__' :
    main(sys.argv[1:])
