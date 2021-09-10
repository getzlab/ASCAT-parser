import pandas
import numpy
import scipy
import scipy.stats
import sys


def is_between(s, lower, upper):
    return (lower <= s) & (upper >= s)

def parser(caveman_path, copynumber_path, copynumber_normal_path):
    caveman = pandas.read_csv(caveman_path,
                              names=['chromosome', 'start_bp', 'end_bp', 'cp_major_normal', 'cp_minor_normal',
                                     'cp_major_tumor', 'cp_minor_tumor'])
    
    cp = pandas.read_csv(copynumber_path, sep='\t')
    cp_normal = pandas.read_csv(copynumber_normal_path, sep='\t')
    
    # remove rows with nan in segmented baf
    cp = cp[~(cp['segmented LogR'].isna() | cp['segmented BAF'].isna())]
    
    def hets(tumor, normal):
        def foo(x):
            thres = 0.1
            index = list()
            for i, row in x[[c for c in x.columns if 'Count' in c]].iterrows():
                temp = sorted(row, reverse=True)
                
                if sum(temp) == 0 or temp[2]/sum(temp) > 0.01:
                    index.append(False)
                    continue
                    
                temp = temp[:2]
                phet = scipy.stats.beta.cdf(0.5, temp[0] + 1, temp[1] + 1)
                index.append( ((phet > thres) or ((1-phet) > thres)) and (sum(temp) >= 20))
                
            return x[index]
                

        tt = tumor['Chromosome'].astype(str).str.cat(tumor['Position'].astype(str), sep='-')
        
        # select hets
        tn = foo(normal)
        
        # format
        tn = tn['#CHR'].astype(str).str.cat(tn['POS'].astype(str), sep='-')
        
        return tt.isin(tn)
    
    cp = cp[hets(cp, cp_normal)]
    cp = cp[~cp['Chromosome'].astype(str).str.contains("X|Y|M")]
    
    parsed = pandas.DataFrame(
        columns=['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'length', 'n_hets', 'f',
                 'tau', 'sigma.tau', 'mu.minor', 'sigma.minor', 'mu.major', 'sigma.major','cp_major_tumor', 'cp_minor_tumor'])

    for i, row in caveman.iterrows():
        temp = pandas.DataFrame(
            columns=['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'length', 'n_hets', 'f',
                     'tau', 'sigma.tau', 'mu.minor', 'sigma.minor', 'mu.major', 'sigma.major', 'cp_major_tumor',
                     'cp_minor_tumor'])

        if row['chromosome'] in ['X', 'Y', 'MT']:  # redundant legacy check
            continue

        selected = cp[(cp['Position'] >= row['start_bp']) & (cp['Position'] <= row['end_bp']) & 
                      (cp['Chromosome'] == int(row['chromosome']))
                     & (cp['BAF'] < 0.95) & (cp['BAF'] > 0.05)] # remove extreme values

        selected.loc[selected['BAF'] > 0.5, 'BAF'] = 1 - selected[selected['BAF'] > 0.5]['BAF'].to_numpy()
        
        aux = {'Chromosome': int(row['chromosome']),
               'Start.bp': row['start_bp'],
               'End.bp': row['end_bp'],
               'n_probes': cp[(cp['Chromosome'].astype('int') == int(row['chromosome'])) & 
                                   (cp['Position'] >= row['start_bp']) & 
                                   (cp['Position'] <= row['end_bp']) & (cp['BAF'] <= 1) & (cp['BAF'] >= 0) ].shape[0],
               'n_hets': selected.shape[0],
               'f': selected['segmented BAF'].mean(), # mean if < 0.5 + mean 1-baf if >0.5
               'tau': 2 * 2 ** selected['segmented LogR'].mean(),
               'cp_major_tumor': row['cp_major_tumor'],
               'cp_minor_tumor': row['cp_minor_tumor']
        }

        parsed = parsed.append(aux.copy(), ignore_index=True, verify_integrity=True)

    parsed = parsed[parsed['n_hets'] > 1]  # cleanup missed sites

    parsed['mu.minor'] = parsed['tau'] * parsed['f']  
    parsed['mu.major'] = parsed['tau'] * (1 - parsed['f'])  
    parsed['length'] = parsed['End.bp'] - parsed['Start.bp']

    means = parsed[(is_between(parsed['tau'], 1.8, 2.2))]['tau']
    print(parsed)
    denominator = parsed['n_probes'].sum() * ((parsed['n_probes'] != 0).sum() - 1) / (parsed['n_probes'] != 0).sum()

    t = numpy.sqrt(
                    (parsed[(is_between(parsed['tau'], 1.8, 2.2))]['tau'] - means.mean()
                   ).pow(2).multiply(parsed['n_probes'], axis=0).sum() / denominator)

    t = t * numpy.sqrt(means.shape[0] - 1)

    denominator = parsed['n_probes']
    denominator[denominator == 0] = 1
    parsed['sigma.tau'] = t / denominator
    parsed['sigma.minor'] = parsed['sigma.tau'].mul(parsed['f'])
    parsed['sigma.major'] = parsed['sigma.tau'].mul(1-parsed['f'])

    parsed['sample'] = 'test'

    del parsed['cp_major_tumor']
    del parsed['cp_minor_tumor']
    
    parsed[["Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets"]] = parsed[["Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets"]].astype('int')

    return parsed
    
if __name__ == '__main__':
    parser(sys.argv[1], sys.argv[2], sys.argv[3]).to_csv(sys.argv[4], sep='\t')
