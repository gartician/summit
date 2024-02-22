from os import makedirs
from pandas import DataFrame, concat, read_table
from os.path import basename, dirname, exists
from numpy import zeros, vectorize, where
from sys import argv, exit
from json import dump
from datetime import date
from pybedtools import BedTool
from scipy.stats import binomtest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from pysam import AlignmentFile
from bisect import bisect_left, bisect_right
from timeit import default_timer as timer

VERSION="1.0.0"

def parse_args():

    # Define program arguments
    parser = ArgumentParser(
        formatter_class=RawDescriptionHelpFormatter,
        description=
        """
GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing data. 
GoPeaks by default works best with narrow peaks such as H3K4me3 and
transcription factors. GoPeaks can be used with the "--broad" flag to
call broad peaks like H3K27Ac/H3K4me1. We encourage users to explore
the parameters of GoPeaks to analyze their data.
        """
    )

    # introduce argument groups
    required = parser.add_argument_group('Required Arguments')
    params = parser.add_argument_group('Summit Parameters')
    optional = parser.add_argument_group('Optional Arguments')

    # input/output file options
    required.add_argument("-b", "--bam", help = "Input BAM file with histone mark or protein-of-interest", required = True, type = str)
    required.add_argument("-o", "--output", help = "Output prefix e.g. data/peaks/sample-1", required = True, type = str)

    # essential program parameters
    params.add_argument("-c", "--control", help = "Control BAM file", required = False, type = str)
    params.add_argument("-m", "--minreads", help = "Analyze bins greater than mindreads. Default: 15", type = int, default = 15)
    params.add_argument("-d", "--mindist", help = "Merge bins within <mindist> base pairs. Default: 1000", type = int, default = 1000)
    params.add_argument("-w", "--minwidth", help = "Remove peaks less than <mindwidth> base pairs wide. Default: 150", type = int, default = 150)
    params.add_argument("-t", "--step", help = "Bin size for genome bins. Default: 100", type = int, default = 100)
    params.add_argument("-l", "--slide", help = "Slide size for genome bins. Default: 50", type = int, default = 50)
    params.add_argument("-p", "--pval", help = "Filter away bins above <pval> significance level. Default: 0.05", type = float, default = 0.05)
    params.add_argument("-s", "--chromsize", help = "Tab-delimited file containing chromosome name and sizes.", type = str)

    # optional program parameters
    optional.add_argument("--broad", help = "Call broad peaks with step=5000, slide=1000, mindist=3000.", action = "store_true", default = False)
    optional.add_argument("--chromosomes", help = "Use a regex string to subset chromosomes.", type = str, required=False)
    optional.add_argument("--canonical", help = "Subset canonically-named chromosomes. Default: 'chr[0-9XY]+$'", nargs="?", type = str, required=False, const = "chr[0-9XY]+$")
    optional.add_argument("-v", "--verbose", help = "Run the program in verbose mode", action = "store_true", default=False)
    optional.add_argument("--version", help = "Display Summit version", action = "store_true", default=False)
    # optional.add_argument("--counts", help = "Output prefix for normalized bin counts", action = "store_true", default = False)

    args = parser.parse_args()

    return(args)

def generate_bins(bam: AlignmentFile, step: int):

    """

    Given a SAM/BAM header, generate a dictionary of regions.
    Genome tracks are often represented as python dictionaries where dicts 
    can represent bins like this...

    FORMAL DEFINITION
    {
        'chrN': [0, step, step + slide, step + 2slide, step + 3slide, etc.]
    }

    Note that the list contains only the start coordinate of a bin. This reduces the size of the list
    and one can easily calculate the end coordinate with an addition function.

    EXAMPLE
    
    with step=150 and slide=100...

    {
        'chr1': [0, 150, 250, 350, 450, etc.],
        'chr2': [0, 150, 250, 350, 450, etc.]
    }

    This python dictionary where key is a chromosome and the value as a list
    is a common structure used throughout this program.

    """

    # filter chromosomes if needed
    if args.chromosomes or args.canonical:

        if args.chromosomes != None:
            pattern = args.chromosomes
        else:
            pattern = args.canonical

        import re
        reg = re.compile(pattern)
        valid = [ i for i, v in enumerate(bam.references) if reg.match(v) ]
        references = [ bam.references[i] for i in valid  ]
        lengths = [ bam.lengths[i] for i in valid  ]

    # use provided tab-delimited file for chromosome names and sizes
    elif args.chromsize:
        cs = read_table(args.chromsize, sep="\t")
        references = cs.loc[:, 1]
        lengths = cs.loc[:, 2]

    # use provided chromosome names and sizes provided in the bam file (recommended)
    else:
        references = bam.references
        lengths = bam.lengths

    # Generate dictionary of regions
    bins = dict()

    # loop across chromosomes
    for chr, last_bp in zip(references, lengths):

        # loop across regions within a chromosome
        regions = []

        for j in range(0, last_bp // step):

            start = j * step
            regions.append(start)

        bins[chr] = regions

    # 8 --> 4 seconds to build the bins for the human genome!
    return(bins)

def normalize_tracks(ttracks: dict, ctracks: dict, bins: dict):

    norm_track = init_bin_values(bins)
    
    # vectorize normalization function
    vecnorm = vectorize(norm)

    for chr in bins.keys():
        norm_track[chr] = vecnorm(ttracks[chr], ctracks[chr])
    
    return(norm_track)

def norm(treatment, control):
    if treatment > control:
        count = round(treatment * (1 - cpm(control)/cpm(treatment)))
        return(count)
    else:
        return(0.0)

def cpm(count: float):
    return(count / 1_000_000)

def count_tracks(bam: AlignmentFile, bins: dict, nreads: bool):

    """
    If a read overlaps a bin, that bin is incremented by one. If a read overlaps
    multiple bins (often the case), then each bin gets an increment of one.

    We cannot expect to exactly replicate GoPeaks counting scheme because
    replicating the BAM-processing library pbenner/gonetics is impossible.
    """

    # initialize zero-filled pileups for all bins
    pileups = init_bin_values(bins)

    # loop over read pairs
    i=0
    for r1, r2 in read_pair_generator(bam):
        
        if r1.is_forward & r2.is_reverse:
            r_start, r_end = r1.reference_start, r2.reference_end
        elif r1.is_reverse & r2.is_forward:
            r_start, r_end = r2.reference_start, r1.reference_end
        else:
            exit(f"Read pair R1 {r1} and R2 {r2} are not complementary! Either both are forward or reverse reads!")

        # define chromosome
        chr = r1.reference_name
        
        # for every read pair, find the left and rightmost bins
        leftmost_bindex = bisect_left(bins[chr], r_start)
        rightmost_bindex = bisect_right(bins[chr], r_end)

        # increment the affected bins by one
        pileups[chr][leftmost_bindex: rightmost_bindex] += 1
        i += 1

        # dev options
        if args.verbose:
            print(bam.filename.decode(),
                r1.reference_name,
                r1.reference_start,
                r1.reference_end,
                r1.is_forward,
                r2.reference_name,
                r2.reference_start,
                r2.reference_end,
                r2.is_reverse,
                i,
                sep="\t"
            )

    if nreads:
        return(pileups, i)
    else:
        return(pileups)

def read_pair_generator(bam, region_string=None):

    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """

    if args.chromosomes or args.canonical:

        if args.chromosomes != None:
            pattern = args.chromosomes
        else:
            pattern = args.canonical

        import re
        reg = re.compile(pattern)
        valid = [ i for i, v in enumerate(bam.references) if reg.match(v) ]
        references = [ bam.references[i] for i in valid  ]
        check_chrs = True

    else:
        check_chrs = False

    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):

        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue

        qname = read.query_name
        
        # skip a read if chromosome filter is set
        if check_chrs:
            if not read.reference_name in references:
                continue

        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def init_bin_values(bins: dict, zero = True):

    # retrieve chromosome names
    chrs = bins.keys()

    # create empty list OR np.array of zeroes per bin for every chromosome
    bin_values = dict()
    for chr in chrs:
        if zero:
            bin_values[chr] = zeros(len(bins[chr]))
        else:
            bin_values[chr] = []
    
    return(bin_values)

def filterBins(tracks: dict, minreads: int):

    """
    Get the indices of all the bins that have read signals
    Get the number of bins that have read signal
    """

    filteredBins = init_bin_values(bins, zero=False)
    numNZBins = []
    
    # determine all non-zero signal bins
    for chr in tracks.keys():
        filteredBins[chr] = where(tracks[chr] > minreads)[0]
        # numNZBins += len(filteredBins[chr])
        numNZBins.append( len(where(tracks[chr] > 0)[0]) )

    numNZBins = sum(numNZBins)
    # it's kinda weird to divide by the number of nz bins when we only analyze bins > minreads no?
    # maybe that's just an artifact of the binomial test, and we can find other distributions

    return(filteredBins, numNZBins)

def nzSignal(filteredBins: dict, tracks: dict):

    """
    Calculate average signal in non-zero bins
    """

    # calculate global signal at non-zero bins
    nzSignal = 0
    numfilteredBins = 0

    # count the number of qualifying bins for FDR adjustment later...
    for chr in tracks.keys():
        for bin in filteredBins[chr]:
            nzSignal += tracks[chr][bin]
            numfilteredBins += 1

    meanNZSignal = nzSignal / numfilteredBins

    return(meanNZSignal)

def callpeaks(tracks: dict, bins: dict, nreads: int):

    # pre peak calling
    filteredBins, minBins = filterBins(tracks, args.minreads)
    meanSignal = nzSignal(filteredBins, tracks)

    if args.verbose:
        print("mu:", meanSignal)
        print("# candidate bins:", minBins)

    # define distribution parameters
    p = float(meanSignal)/float(minBins)/float(nreads)

    # loop through bins and test for significance
    # aggregate results into a DataFrame for each bin
    chrs = []
    starts = []
    ends = []
    pvals = []
    signals = []
    for chr in filteredBins.keys():
        for bin in filteredBins[chr]:

            start = bins[chr][bin]
            signal = tracks[chr][bin]
            pval = binomtest(k = int(signal), n = int(nreads), p = float(p), alternative='greater').pvalue

            chrs.append(chr)
            starts.append(start)
            ends.append(start + args.slide)
            pvals.append(pval)
            signals.append(signal)
    
    df = DataFrame({
        'chr': chrs,
        'start': starts,
        'end': ends,
        'pval': pvals,
        'signal': signals
    })

    # adjust pvals of bins
    df = fdr(df, nreads)

    # stitch peaks together
    intervals = BedTool.from_dataframe(df).merge(d=args.mindist).to_dataframe()

    # filter peaks by minimum width
    intervals = intervals[intervals["end"] - intervals["start"] > args.minwidth]

    return(intervals)

def fdr(df: DataFrame, nreads: int):

    df = df.sort_values('pval')

    # assign pval ranks by UNIQUE pvals
    pval_map = dict()
    for i, v in enumerate(df['pval'].unique()):
        pval_map[v] = i + 1
    df["rank"] = df["pval"].apply(lambda x: pval_map[x])

    # calculate fdr-adjusted pvals
    df["pval"] = (df["pval"] * nreads) / df["rank"]
    df["rank"] = None

    # filter away insignificant bins
    df = df[df["pval"] < args.pval]

    # re-arrange df into bed-like format
    df = df.sort_values(['chr', 'start', 'end'])

    return(df)

def check_tracks(ttrack, ctrack):
    if not ttrack.keys() == ctrack.keys():
        exit("ERROR: non-matching chromosomes in the treatment and control tracks!")

def export_bins(df: DataFrame, output: str):
    
    # define output files
    outpeaks = f"{output}-peaks.bed"

    # create folders if needed
    dname = dirname(output)
    if dname:
        if not exists(dname):
            makedirs(dname, 770, exist_ok=True)

    # export bed file
    df.to_csv(outpeaks, sep = "\t", index = False, header = False)    

def export_counts(tracks, bins):

    # define output file
    outfile = f"{args.output}-counts.txt.gz"

    # put each chr into a DF into a list
    chrs = [ i for i in bins.keys() ]

    dfs = []
    for chr in chrs:
        starts = bins[chr]
        counts = tracks[chr]
        d = DataFrame({
            "chr": chr,
            "start": starts,
            "count": counts
        })
        d["end"] = d["start"].apply(lambda x: x + args.slide)
        d = d[["chr", "start", "end", "count"]]
        dfs.append(d)

    df = concat(dfs)
    df.to_csv(outfile, sep = "\t", index=False)

def export_log(peaks):

    end = timer()

    log = {
        "Version": VERSION,
        "Date": date.today().strftime("%b-%d-%Y"),
        "Elapsed": end - start,
        "Sample": basename(args.bam).split(".")[0],
        "Command": " ".join(argv),
        "Peaks": peaks.shape[0]
    }

    outjson = f"{args.output}-summit.json"
    with open(outjson, "w") as f:
        dump(log, f)

def preflight():

    # override default step and slide with broad flag
    if args.broad:
        args.step = 5000
        args.slide = 1000
        args.mindist = 3000

    # display package version
    if args.version:
        print("Summit Version: {}".format(VERSION))
        exit(0)

    # return starting time
    start = timer()

    return(start)

if __name__ == "__main__":

    args = parse_args()

    start = preflight()

    # import treatment sample
    bam = AlignmentFile(args.bam)
    
    # Generate bins for the genome with step and slide
    bins = generate_bins(bam, args.step)

    # Get tracks for treatment
    ttracks, nreads = count_tracks(bam, bins, nreads = True)

    if args.control != None:

        # import ctrl sample
        ctrl = AlignmentFile(args.control)

        # Get track for control if provided
        ctracks = count_tracks(ctrl, bins, nreads = False)

        # Normalize treatment to control
        norm_tracks = normalize_tracks(ttracks, ctracks, bins)

        # pre-flight check
        check_tracks(ttracks, ctracks)

        # if args.counts:
        #     export_counts(norm_tracks, bins)

        # call peaks on normalized tracks
        peaks = callpeaks(tracks = norm_tracks, bins = bins, nreads = nreads)
    
    else:

        # else call peaks on single treatment sample track
        peaks = callpeaks(tracks = ttracks, bins = bins, nreads = nreads)

    export_bins(peaks, args.output)
    export_log(peaks)