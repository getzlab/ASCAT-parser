import pandas
import numpy
import scipy
import scipy.stats
import sys
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import minimize
#import statistics

def parse_args(argv):
    output = argparse.ArgumentParser("ASCAT to Allelic Capseg/GISTIC format parser.")

    # file with CNV segments
    output.add_argument('caveman_segments', type=str, help="ASCAT tumor copynumber.caveman.csv output.")
    # file with target-het-level data
    output.add_argument('copynumber', type=str, help="ASCAT tumor copynumber.txt.gz output.")
    # file with normal het sites
    output.add_argument('normal_hets', type=str, help="ASCAT normal het site output.")
    # approximate read depth
    output.add_argument('depth', type=int, help="approximate read depth.",default=100)

    # output acs file for ABSOLUTE
    output.add_argument('allelic_output', type=str, help="Path for allelic capseg output.")
    # output seg file for IGV and GISTIC
    output.add_argument('gistic_output', type=str, help="Path for seg file output.")
    
    # ignore het site BAFs and use segmentedBAF
    # **** WARNING: this option seriously undercalls het site allele shifts ****
    output.add_argument("--use-ascat-normalization", action="store_true", help="Use segmented BAF for f calculation.")
    
    return output.parse_args(argv)
    
def helper(args):
    return parser(args.caveman_segments, args.copynumber, args.normal_hets, args.depth, args.use_ascat_normalization)

def is_between(s, lower, upper):
    return (lower <= s) & (upper >= s)

def negLL_beta2(parameters,v):
  maf,dep = parameters
  tiny = 1e-30
  a = dep*maf+1
  b = dep*(1-maf)+1
  p1 = numpy.exp(scipy.stats.beta.logpdf(v, a, b))/2+numpy.exp(scipy.stats.beta.logpdf(v, b, a))/2
  neg_LL=-1*numpy.sum(numpy.log(p1+tiny))
  return neg_LL

def beta2fit(baf,DEP):
    NSIGMA = 2
    v1 = baf
    v2 = 1-baf
    v = numpy.concatenate((v1, v2), axis=0)
    par0 = [0.5, DEP/1.5]
    #print(par0)
    mle_model = minimize(negLL_beta2, par0, args=v, method='Nelder-Mead')
    par0 = mle_model.x
    if (par0[0] > 0.5):
        par0[0] = 1 - par0[0]

    # Inverse Hessian -> Covariance matrix doesn't seem to work...
    mle_model2 = minimize(negLL_beta2, par0, args=v, method='BFGS')
    maf=mle_model2.x[0]
    dep=mle_model2.x[1]
    mafse = numpy.sqrt(mle_model2.hess_inv[0, 0])
    depse = numpy.sqrt(mle_model2.hess_inv[1, 1])

    # estimate maf CI - fit maf loglike
    ll0 = negLL_beta2([maf,dep], v)

    # maf high CI 2sigma (95%) - log like up by 1 unit
    loghigh=numpy.log10(0.5-maf)
    maf1 = numpy.logspace(-4 ,loghigh, num=100) + maf
    ll1 = maf1*0
    i = 0
    for m1 in maf1:
        ll1[i] = negLL_beta2([m1, dep], v) - ll0 - NSIGMA/2.0
        i = i+1

    khigh = numpy.where(numpy.diff(numpy.sign(ll1)))[0]
    mafH=0.5
    if len(khigh) == 1:
        mafH=float(maf1[khigh])

    # maf low CI
    loglow=numpy.log10(maf)
    maf2 = maf - numpy.logspace(-4 ,loglow, num=100)
    ll2 = maf2*0
    i = 0
    for m1 in maf2:
        ll2[i] = negLL_beta2([m1, dep], v) - ll0 - NSIGMA/2.0
        i=i+1

    klow = numpy.where(numpy.diff(numpy.sign(ll2)))[0]
    mafL=0
    if len(klow) == 1:
        mafL=float(maf2[klow])

    mafCI=[mafL, mafH]

    return [maf,mafCI,dep,depse]

def parser(caveman_path, copynumber_path, copynumber_normal_path, depth, use_ascat_normalization=False):
    # load segments from caveman
    HET_MAF_MIN   = 0.3  # het sites with MAF below not used a het
    HET_ADEP_MIN  = 7    # minimum normal allelic depth for normal het sites
    HET_OTHER_MAX = 0    # max dep for bases not in two two for hets
    TAU_SIGMA0    = 0.01 # added to sigma.tau in quadrature
    F_SIGMA0      = 0.002 # added to het MAF std err in quadrature
    HET_GAP_MAX   = 0.5 # maximum allowed gap in segment of hets (~anti centromere filter)

    caveman_segs = pandas.read_csv(caveman_path,
                              names=['chromosome', 'start_bp', 'end_bp', 'cp_major_normal', 'cp_minor_normal',
                                     'cp_major_tumor', 'cp_minor_tumor'])

    # load tumor het site LogCR and BAF data
    cp = pandas.read_csv(copynumber_path, sep='\t')
    cp = cp[~cp['Chromosome'].astype(str).str.contains("X|Y|M")]
    cp = cp[cp['Chromosome'].astype(int)<=22]
    cp['site'] = cp['Chromosome'].astype(str).str.cat(cp['Position'].astype(str), sep='-')

    # load matched normal het data
    het_normal = pandas.read_csv(copynumber_normal_path, sep='\t')
    het_normal = het_normal[~het_normal['#CHR'].astype(str).str.contains("X|Y|M")]
    het_normal['site'] = het_normal['#CHR'].astype(str).str.cat(het_normal['POS'].astype(str), sep='-')
    # subset normal_hets to valid_hets (depth threshold , MAF threshold, base noise threshold)
    het_normal = het_normal[het_normal.Good_depth>=(2*HET_ADEP_MIN)]
    het = numpy.transpose(numpy.array([het_normal.Count_A, het_normal.Count_C,het_normal.Count_G, het_normal.Count_T ]))
    het = numpy.sort(-1*het,axis=1)
    het = -1*het
    AB = het[:,:2]
    ABD=numpy.amin(AB,axis=1)
    MAF = numpy.divide(AB[:,1],numpy.sum(AB,1))
    OTHER = numpy.sum(het[:,2:],axis=1)
    ok = [MAF[i1]>=HET_MAF_MIN and OTHER[i1]<=HET_OTHER_MAX and ABD[i1]>=HET_ADEP_MIN for i1 in range(len(MAF))]
    #sum(ok)
    het_normal = het_normal[ok]
    cp['het'] = cp['site'].isin(het_normal['site'])
    cp['targ'] = cp['Log R'] > -1e10
    del(het_normal)
    # remove rows with nan in segmented baf
    #cp = cp[~(cp['segmented LogR'].isna() | cp['segmented BAF'].isna())]

    # create acs dataframe
    acs_df = pandas.DataFrame(
        columns=['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'length', 'n_hets', 'f',
                 'tau', 'sigma.tau', 'mu.minor', 'sigma.minor', 'mu.major', 'sigma.major'] ) #,'cp_major_tumor', 'cp_minor_tumor'])

    longest_segment=0
    longest_segment_stdev_tau=0

    for i, row in caveman_segs.iterrows():
        # acs_df1 = pandas.DataFrame(
        #     columns=['Chromosome', 'Start.bp', 'End.bp', 'n_probes', 'length', 'n_hets', 'f',
        #              'tau', 'sigma.tau', 'mu.minor', 'sigma.minor', 'mu.major', 'sigma.major']) #, 'cp_major_tumor','cp_minor_tumor'])

        if row['chromosome'] in ['X', 'Y', 'MT']:  # redundant legacy check
            continue

        cpseg = cp[(cp['Position'] >= row['start_bp']) & (cp['Position'] <= row['end_bp']) &
                      (cp['Chromosome'] == int(row['chromosome']))]

        # check for gap fraction of segment.
        gap=numpy.amax(numpy.diff(cpseg[cpseg['het']]['Position'] ) )
        xgap = float(gap)/(row['end_bp']-row['start_bp'])
        if xgap>HET_GAP_MAX:
            continue
        NTARG = cpseg['targ'].shape[0]
        if NTARG>longest_segment:
            longest_segment=NTARG
            longest_segment_stdev_tau = numpy.std(2 * 2 ** cpseg[cpseg['targ']]['Log R'])

        DEP = depth
        # baf = cpseg['het']['BAF']
        stdev=numpy.std(cpseg[cpseg['het']]['BAF'])
        fse=stdev/numpy.sqrt(cpseg[cpseg['het']].shape[0])
        #
        [maf1,mafCI,dep1,depse1]=beta2fit(cpseg[cpseg['het']]['BAF'], DEP)
        #
        if (mafCI[1]) >= 0.4999:
            maf1 = 0.5
        # [fse, numpy.diff(mafCI)]
        fse2=numpy.diff(mafCI)[0]/2.0
        fse2=numpy.sqrt(numpy.power(fse2,2)+numpy.power(F_SIGMA0,2))
        tau1=2 * 2 ** cpseg['segmented LogR'].median()
        mu_minor1 = tau1 * maf1
        mu_major1 = tau1 * (1-maf1)

        acs_df1 = {'Chromosome': int(row['chromosome']),
               'Start.bp': row['start_bp'],
               'End.bp': row['end_bp'],
               'n_probes': cpseg[cpseg['targ']].shape[0],
               'n_hets': cpseg[cpseg['het']].shape[0],
               'f': maf1, #'f': selected['segmented BAF'].mean() if use_ascat_normalization else selected['BAF'].mean(), # mean if < 0.5 + mean 1-baf if >0.5
               'tau':tau1,
               'mu.major': mu_major1,
               'mu.minor': mu_minor1,
               'sigma.tau': 0,
               'sigma.major': fse2,
               'sigma.minor': fse2
               }
        print(i,maf1,fse,mafCI,float(numpy.diff(mafCI)),dep1,depse1)
        acs_df = acs_df.append(acs_df1.copy(), ignore_index=True, verify_integrity=True)
        # if i>2:
        #     break

    acs_df = acs_df[acs_df['n_hets'] > 1]  # cleanup missed sites

    acs_df['length'] = acs_df['End.bp'] - acs_df['Start.bp']

    acs_df['sigma.tau'] = longest_segment_stdev_tau / numpy.sqrt(acs_df['n_probes'].astype('float'))
    acs_df['sigma.tau'] = numpy.sqrt( numpy.power(acs_df['sigma.tau'],2) + numpy.power(TAU_SIGMA0,2) )

    tau_relerr=acs_df['sigma.tau'].div(acs_df['tau'])
    f_relerr=acs_df['sigma.major'].div(acs_df['f'])
    acs_df['sigma.minor'] = acs_df['mu.minor'].mul(numpy.sqrt(numpy.power(tau_relerr,2) + numpy.power(f_relerr,2)))
    acs_df['sigma.major'] = acs_df['mu.major'].mul(numpy.sqrt(numpy.power(tau_relerr,2) + numpy.power(f_relerr,2)))

    acs_df['sample'] = caveman_path.split('/')[-1].split('.')[0]

    acs_df[["Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets"]] = acs_df[["Chromosome", "Start.bp", "End.bp", "n_probes", "length", "n_hets"]].astype('int')

    seg_df = acs_df.copy()

    seg_df['mean_log2_copy_ratio'] = list(map(lambda x: numpy.log( x / 2 ) / numpy.log(2), list(seg_df['tau'])))

    seg_df = seg_df[['sample','Chromosome','Start.bp','End.bp','n_probes','mean_log2_copy_ratio']]

    return acs_df, seg_df, cp


def xhg38(a,p):
    L = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616,
    64444167, 46709983, 50818468, 156040895, 57227415, 16569]
    a1 = list(a)
    p1 = list(p)
    C = [1]+list(numpy.cumsum(L))
    #x = 0*a
    x = [C[a1[i]-1]+p1[i] for i in range(len(a1))]
    return x

def create_plots(args, acs_df, cp):
    aL = [i+1 for i in range(23)]
    pL = [0 for i in range(23)]
    xL = xhg38(aL, pL)
    xM = [(xL[i] + xL[i + 1]) / 2 for i in range(22)]
    aL=aL[:22]
    pL=pL[:22]

    plt.figure(figsize=(18, 10), dpi=100)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)


    x = xhg38(cp[cp['targ']]['Chromosome'].astype(int), cp[cp['targ']]['Position'])
    y = list(cp[cp['targ']]['Log R'])
    ax1.scatter(x, y, s = 0.5, alpha=0.05 )
    ax1.set_title(args.caveman_segments.split('/')[-1].split('.')[0] ) #+ " Log R")
    ax1.set_xlabel('genomic position')
    ax1.set_ylabel('Log R')
    ax1.set_xlim(0, xL[-1])
    ax1.set_ylim(-4, 4)
    ax1.set_xticks(xM)
    ax1.set_xticklabels(aL)
    for i in range(0,22,2):
        p = plt.Rectangle((xL[i], -4), xL[i+1]-xL[i], 8, alpha=0.2, edgecolor=None, facecolor=(0.7, 0.7, 0.7))
        ax1.add_patch(p)

    x1 = xhg38(acs_df['Chromosome'].astype(int), acs_df['Start.bp'])
    x2 = xhg38(acs_df['Chromosome'].astype(int), acs_df['End.bp'])
    y1 = numpy.log2(list(acs_df['tau']/2))
    y2 = y1
    #hold(True)
    for i in range(acs_df['Chromosome'].shape[0]):
        ax1.plot([x1[i],x2[i]],[y1[i],y2[i]], ls='-', color = 'k', linewidth = 3)

    ax1.grid(True)
    #hold(False)

    
    x = xhg38(cp[cp['het']]['Chromosome'].astype(int), cp[cp['het']]['Position'])
    y = list(cp[cp['het']]['BAF'])
    ax2.scatter(x, y, s = 0.5, alpha=0.05 )  # TODO plot segmented or normal BAF?
    #plt.title(args.caveman_segments.split('/')[-1].split('.')[0] + " BAF (B-Allele Fraction)")
    ax2.set_xlabel('genomic position')
    ax2.set_ylabel('BAF')
    ax2.set_xlim(0, xL[-1])
    ax2.set_ylim(0, 1)
    ax2.set_xticks(xM)
    ax2.set_xticklabels(aL)
    for i in range(0,22,2):
        p = plt.Rectangle((xL[i], 0), xL[i+1]-xL[i], 1, alpha=0.2, edgecolor=None, facecolor=(0.7, 0.7, 0.7))
        ax2.add_patch(p)
    x1 = xhg38(acs_df['Chromosome'].astype(int), acs_df['Start.bp'])
    x2 = xhg38(acs_df['Chromosome'].astype(int),acs_df['End.bp'])
    y1 = list(acs_df['f'])
    y2 = y1
    #hold(True)
    for i in range(acs_df['Chromosome'].shape[0]):
        ax2.plot([x1[i],x2[i]],[y1[i],y2[i]], ls='-', color = 'k', linewidth = 3)

    ax2.grid(True)
    #hold(False)
    #plt.show()

    plt.savefig(args.caveman_segments.split('/')[-1].split('.')[0] + ".ascat.acs.png")
    
    
if __name__ == '__main__':
    #expect caveman_tumor, copynumber_tumor, copynumber_normal, allelic_output_path, gistic_output_path
    args = parse_args(sys.argv[1:])
    
    print("args:")
    print(args)

    allelic_capseg, gistic, cp  = helper(args)
    
    allelic_capseg.to_csv(args.allelic_output, sep='\t', index=False, float_format="%.5f")
    
    gistic.to_csv(args.gistic_output, sep='\t', index=False, float_format="%.5f")
    
    create_plots(args, allelic_capseg, cp)

