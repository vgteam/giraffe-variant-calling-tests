import functools
import tempfile
import os

configfile: "lr-config.yaml"

# Where are the input graphs?
#
# For each reference (here "chm13"), this directory must contain:
#
# hprc-v1.1-mc-chm13.d9.gbz
# hprc-v1.1-mc-chm13.d9.dist
#
# Also, it must either be writable, or already contain zipcode and minimizer
# indexes for each set of minimizer indexing parameters (here "k31.w50.W"),
# named like:
#
# hprc-v1.1-mc-chm13.d9.k31.w50.W.withzip.min
# hprc-v1.1-mc-chm13.d9.k31.w50.W.zipcodes
#
# And for haplotype sampling, the full graph (without any .dXX) must be available:
#
# hprc-v1.1-mc-chm13.gbz
#
# A top-level chains only "distance" index will be made if not present:
#
# hprc-v1.1-mc-chm13.tcdist
#
# As will the haplotype sampling indexes:
#
# hprc-v1.1-mc-chm13.ri
# hprc-v1.1-mc-chm13.fragment.hapl
#
GRAPHS_DIR = config.get("graphs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs"

#
# For SV calling, we need a truth vcf and a BED file with confident regions
#
# giab6_chm13.vcf.gz
# giab6_chm13.confreg.bed
#
# And a tandem repeat bed for sniffles
# human_chm13v2.0_maskedY_rCRS.trf.bed (from https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed)
#
# TODO: For HG001 we might need https://platinum-pedigree-data.s3.amazonaws.com/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.sort.oa.vcf.gz

SV_DATA_DIR = config.get("sv_data_dir", None) or "/private/home/jmonlong/workspace/lreval/data"

# Where are the reads to use?
#
# This directory must have "real" and "sim" subdirectories. Within each, there
# must be a subdirectory for the sequencing technology, and within each of
# those, a subdirectory for the sample.
#
# For real reads, each sample directory must have a ".fq.gz" or ".fastq.gz" file.
# The name of the file must contain the sample name. If the directory is not
# writable, and you want to trim adapters off nanopore reads, there must also
# be a ".trimmed.fq.gz" or ".trimmed.fastq.gz" version of this file, with the
# first 100 and last 10 bases trimmed off. The workflow will generate
# "{basename}.{subset}.fq.gz" files for each subset size in reads ("1k", "1m:,
# etc.) that you want to work with.
#
# For simulated reads, each sample directory must have files
# "{sample}-sim-{tech}-{subset}.gam" for each subset size as a number (100, 1000,
# 1000000, etc.) that you want to work with. If the directory is not writable,
# it must already have abbreviated versions ("1k" or "1m" instead of the full
# number) of the GAM files, and the corresponding extracted ".fq" files.
#
# There can also be a "{sample}-sim-{tech}.{category}.txt" file with the names
# of reads in a gategory (like "centromeric") for analysis of different subsets
# of reads.
#
# Simulated reads should be made with the "make_pbsim_reads.sh" script in this
# repository, or, for paired-end read simulation for Illumina reads:
#
# vg sim -r -n 2500000 -a -s 12345 -p 570 -v 165 -i 0.00029 --multi-position
#
# Using an XG and a GBWT for the target sample-and-reference graph,
# --sample-name with the target sample's name, and the full real read set for
# the training FASTQ.
#
# A fully filled out reads directory might look like:
#.
#├── real
#│   ├── hifi
#│   │   └── HG002
#│   │       ├── HiFi_DC_v1.2_HG002_combined_unshuffled.1k.fq
#│   │       └── HiFi_DC_v1.2_HG002_combined_unshuffled.fq.gz
#│   └── r10
#│       └── HG002
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.fastq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.fastq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.10k.fq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1k.fq.gz
#│           └── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1m.fq.gz
#└── sim
#    ├── hifi
#    │   └── HG002
#    │       ├── HG002-sim-hifi.centromeric.txt
#    │       ├── HG002-sim-hifi-1000.gam
#    │       ├── HG002-sim-hifi-10000.gam
#    │       ├── HG002-sim-hifi-1000000.gam
#    │       ├── HG002-sim-hifi-10k.fq
#    │       ├── HG002-sim-hifi-10k.gam
#    │       ├── HG002-sim-hifi-1k.fq
#    │       ├── HG002-sim-hifi-1k.gam
#    │       ├── HG002-sim-hifi-1m.fq
#    │       └── HG002-sim-hifi-1m.gam
#    └── r10
#        └── HG002
#            ├── HG002-sim-r10.centromeric.txt
#            ├── HG002-sim-r10-1000.gam
#            ├── HG002-sim-r10-10000.gam
#            ├── HG002-sim-r10-1000000.gam
#            ├── HG002-sim-r10-10k.fq
#            ├── HG002-sim-r10-10k.gam
#            ├── HG002-sim-r10-1k.fq
#            ├── HG002-sim-r10-1k.gam
#            ├── HG002-sim-r10-1m.fq
#            └── HG002-sim-r10-1m.gam
#
READS_DIR = config.get("reads_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads"

# Where are the linear reference files?
#
# For each reference name (here "chm13") this directory must contain:
#
# A FASTA file with PanSN-style (CHM13#0#chr1) contig names: 
# chm13-pansn-newY.fa
#
# (For grch38 we don't use the -newY part.)
#
# For the calling references (chm13v2.0 and grch38) we also need a plain .fa
# and .fa.fai without pansn names, and _PAR.bed files with the pseudo-autosomal
# regions.
#
# We also use, but can generate:
#
# Index files for Minimap2 for each preset (here "hifi", can also be "ont" or "sr", and can be generated from the FASTA):
# chm13-pansn-newY.hifi.mmi
# 
# A Winnowmap repetitive kmers file:
# chm13-pansn-newY.repetitive_k15.txt
#
# Minimap2 and BWA-MEM indexes
#
REFS_DIR = config.get("refs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references"

# Where are variant call truth set files kept (for when there isn't a handy hosted URL somewhere).
# These are organized by reference and then sample
TRUTH_DIR = config.get("truth_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/truth-sets"

# Where are custom trained DeepVariant models for Giraffe?
# Models should be in directories by tech and then by date in YYYY-MM-DD format.
MODELS_DIR = config.get("models_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/"

# What stages does the Giraffe mapper report times for?
STAGES = ["minimizer", "seed", "tree", "fragment", "chain", "align", "winner"]

# What stages does the Giraffe mapper report times for on the non-chainign codepath?
NON_CHAINING_STAGES = ["minimizer", "seed", "cluster", "extend", "align", "pairing", "winner"]

# To allow for splitting and variable numbers of output files, we need to know
# the available subset values to generate rules.
KNOWN_SUBSETS = ["100", "1k", "10k", "100k", "1m"]
CHUNK_SIZE = 10000

# For each Slurm partition name, what is its max wall time in minutes?
# TODO: Put this in the config
SLURM_PARTITIONS = [
    ("short", 60),
    ("medium", 12 * 60),
    ("long", 7 * 24 * 60)
]

# Where is a large temp directory?
LARGE_TEMP_DIR = config.get("large_temp_dir", "/data/tmp")

#Different phoenix nodes seem to run at different speeds, so we can specify which node to run
#This gets added as a slurm_extra for all the real read runs
REAL_SLURM_EXTRA = config.get("real_slurm_extra", None) or ""

# How many threads do we want mapping to use?
MAPPER_THREADS = 32

# We may not want to populate the MiniWDL task cache because it makes us take
# more shared disk space.
FILL_WDL_CACHE = "true" if config.get("fill_wdl_cache", True) else "false"

wildcard_constraints:
    expname="[^/]+",
    refgraphbase="[^/]+?",
    refgraph="(?!model)(?!legacy)(?!olddv)(?!newdv)[^/_]+?",
    reference="chm13|grch38|chm13v1",
    # We can have multiple versions of graphs with different modifications and clipping regimes
    modifications="(?!-mc)(-[^.-]+(\\.trimmed)?(\\.clip\\.[0-9]*\\.[0-9]*)?(\\.ec[0-9kmKM]+)?)*",
    clipping="\\.d[0-9]+|",
    full="\\.full|",
    chopping="\\.unchopped|",
    # After a full graph, we might have sampling stuff in a graph name
    sampling="(-[^-]+)*",
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)",
    basename=".+(?<!\\.trimmed)",
    subset="([0-9]+[km]?|full)",
    category="((not_)?(centromeric))?|",
    # We can restrict calling to a small region for testing
    region="(|chr21)",
    # We use this for an optional separating dot, so we can leave it out if we also leave the field empty
    dot="\\.?",
    callparams="(|(\\.model[0-9-]+[a-zA-Z]*|\\.nomodel)?(\\.legacy|\\.olddv|\\.newdv)*)",
    tech="[a-zA-Z0-9]+",
    statname="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    statnamex="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    statnamey="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    realness="(real|sim)",
    realnessx="(real|sim)",
    realnessy="(real|sim)",

def auto_mapping_threads(wildcards):
    """
    Choose the number of threads to use map reads, from subset.
    """
    number = subset_to_number(wildcards["subset"])
    mapping_threads = 0
    if number > 100000:
        mapping_threads = MAPPER_THREADS
    elif number > 10000:
        mapping_threads = 16
    else:
        mapping_threads = 8

    if wildcards.get("mapper", "").startswith("graphaligner") and wildcards.get("realness", "") == "sim":
        #Graphaligner is really slow so for simulated reads where we don't care about time
        #double the number of threads
        #At most 128 because it errors with too many threads sometimes
        return min(mapping_threads * 2, 128) 
    else:
        return mapping_threads

def auto_mapping_slurm_extra(wildcards):
    """
    Determine Slurm extra arguments for a timed, real-read mapping job from subset.
    """
    if exclusive_timing(wildcards):
        return "--exclusive " + REAL_SLURM_EXTRA
    else:
        return REAL_SLURM_EXTRA

def auto_mapping_memory(wildcards):
    """
    Determine the memory to use for Giraffe mapping, in MB, from subset and tech.
    """
    thread_count = auto_mapping_threads(wildcards)

    base_mb = 60000

    if wildcards["tech"] == "illumina" or wildcards["tech"] == "element":
        scale_mb = 200000
    elif wildcards["tech"] == "hifi":
        scale_mb = 240000
    elif wildcards["tech"] == "r10":
        scale_mb = 600000
    else:
        scale_mb = 210000

    # Scale down memory with threads
    return scale_mb / 64 * thread_count + base_mb

def choose_partition(minutes):
    """
    Get a Slurm partition that can fit a job running for the given number of
    minutes, or raise an error.
    """
    for name, limit in SLURM_PARTITIONS:
        if minutes <= limit:
            return name
    raise ValueError(f"No Slurm partition accepts jobs that run for {minutes} minutes")

def remote_or_local(url):
    """
    Wrap a URL as a Snakemake "remote file", but pass a local path through.
    """

    if url.startswith("https://"):
        # Looks like a remote.
        return storage.http(url)
    else:
        return url

def to_local(possibly_remote_file):
    """
    Given a result of remote_or_local, turn it into a local filesystem path.

    Snakemake must have already downloaded it for us if needed.
    """

    # TODO: Do we still have the list problem in Python code with storage
    # plugins?

    if isinstance(possibly_remote_file, list):
        return possibly_remote_file[0]
    else:
        return possibly_remote_file

def subset_to_number(subset):
    """
    Take a subset like 1m or full and turn it into a number.
    """
    if subset == "full":
        return float("inf")
    elif subset.endswith("m"):
        multiplier = 1000000
        subset = subset[:-1]
    elif subset.endswith("k"):
        multiplier = 1000
        subset = subset[:-1]
    else:
        multiplier = 1

    return int(subset) * multiplier

def reference_basename(wildcards):
    """
    Find the linear reference base name without extension from a reference.

    This reference is used for mapping. Calling may be against a different
    calling reference (with possibly a different Y) or renamed contigs.

    This reference will use PanSN contig names.
    """
    parts = [wildcards["reference"], "pansn"]
    if wildcards["reference"] == "chm13":
        # We want to use a version of the reference FASTA with the "new"
        # non-HG002, non-GRCh38 Y contig.
        parts.append("newY")
    elif wildcards["reference"] == "chm13v1":
        # Use this for the older version, so just chm13 without the v1 and without the newY
        parts[0] = "chm13"
    return os.path.join(REFS_DIR, "-".join(parts))

def reference_fasta(wildcards):
    """
    Find the linear reference FASTA from a reference.
    """
    return reference_basename(wildcards) + ".fa"

def reference_dict(wildcards):
    """
    Find the linear reference FASTA dictionary from a reference.
    """
    return reference_fasta(wildcards) + ".dict"

def reference_path_list_callable(wildcards):
    """
    Find the path list file for a linear reference that we can actually call on, from reference and region.
   
    We "can't" call on chrY for CHM13 because the one we use in the graphs
    isn't the same as the one in CHM13v2.0 where the calling happens.
    """
    return reference_fasta(wildcards) + ".paths" + wildcards.get("region", "") + ".callable.txt"

def reference_prefix(wildcards):
    """
    Find the PanSN prefix we need to remove to convert form PanSN names to
    non-PanSN names, from reference.
    """
    return {
        "chm13": "CHM13#0#",
        "chm13v1": "CHM13#0#",
        "grch38": "GRCh38#0#"
    }[wildcards["reference"]]

def calling_reference_fasta(wildcards):
    """
    Find the linear reference FASTA with non-PanSN names from a reference (for
    interpreting VCFs).

    For CHM13, we always use CHM13v2.0 as the calling reference since that's
    the one we can get a truth on.
    """
    match wildcards["reference"]:
        case "chm13":
            return os.path.join(REFS_DIR, "chm13v2.0.fa")
        case "chm13v1":
            return os.path.join(REFS_DIR, "chm13v2.0.fa")
        case reference:
            return os.path.join(REFS_DIR, reference + ".fa")

def calling_reference_fasta_index(wildcards):
    """
    Find the index for the linear calling (non-PanSN) reference, from reference.
    """
    return calling_reference_fasta(wildcards) + ".fai"

def calling_reference_restrict_bed(wildcards):
    """
    Find the BED for the linear calling (non-PanSN) reference region we think we can call on, from reference.
    """
    return calling_reference_fasta(wildcards) + ".callable.from." + wildcards["reference"] + ".bed"

def reference_restrict_bed(wildcards):
    """
    Find the BED for the linear PanSN reference region we think we can call on, from reference.

    (We need this for SV calling since it is implemened in PanSN space and not
    on the "calling" references. TODO: refactor to change that.)
    """
    return reference_fasta(wildcards) + ".callable.from." + wildcards["reference"] + ".bed"

def calling_reference_par_bed(wildcards):
    """
    Find the BED for the psuedo-autosomal regions of a reference, from reference.
    """
    return os.path.splitext(calling_reference_fasta(wildcards))[0] + "_PAR.bed"

def haploid_contigs(wildcards):
    """
    Get a list of all haploid contigs in a sample, with no prefixes, from sample.

    Will be either [] or ["chrX", "chrY"].
    """
    XX_SAMPLES = {"HG001"}
    XY_SAMPLES = {"HG002"}
    if wildcards["sample"] in XX_SAMPLES:
        return []
    elif wildcards["sample"] in XY_SAMPLES:
        return ["chrX", "chrY"]
    else:
        raise RuntimeError(f"Unknown karyotype for {wildcards['sample']}")

def wdl_cache(wildcards):
    """
    Get a WDL workflow step cache directory path from root.
    """
    return os.path.join(os.path.abspath("."), wildcards["root"], "miniwdl-cache")

def model_files(wildcards):
    """
    Get the DV model files to use, or None, from tech, mapper, and callparams.
    """

    if wildcards.tech in ("hifi", "r10y2025") and "giraffe" in wildcards.mapper:
        if wildcards.callparams.startswith(".model"):
            # Grab the model name
            model_name = wildcards.callparams.split(".")[1][len("model"):]
        elif wildcards.callparams.startswith(".nomodel"):
            # Just use the model as shipped in DV
            return None
        else:
            if wildcards.tech == "hifi":
                # Use a particular trained model as the default for HiFi.
                if ".olddv" in wildcards.callparams:
                    # Old DV needs the old version of the model
                    model_name = "2025-03-26noinfo"
                elif ".newdv" in wildcards.callparams:
                    # New DV needs the new version with the new model example_info JSON
                    model_name = "2025-03-26"
                else:
                    # We know we are defaulting to old DV until https://ucsc-gi.slack.com/archives/C01D0M09G5D/p1753482919257519?thread_ts=1753201495.785309&cid=C01D0M09G5D is fixed. So use the old style model by default.
                    model_name = "2025-03-26noinfo"
            else:
                # For r10, keep defaulting to no trained model for now.
                return None
        model_base_path = os.path.join(MODELS_DIR, wildcards.tech, model_name)
        # Get all the files as absolute paths in sorted order, minus any README.
        # We want to be the same as the manual order so we don't get WDL cache misses.
        model_files = sorted([os.path.join(model_base_path, f) for f in os.listdir(model_base_path) if f != "README.txt"])
        return model_files
    return None

def truth_vcf_url(wildcards):
    """
    Find the URL or local file for the variant calling truth VCF, from reference and sample.
    """

    if wildcards["sample"] == "HG002":
        # These are available online directly
        return  {
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz",
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        }[wildcards["reference"]]
    elif wildcards["sample"] == "HG001":
        return  {
            # On CHM13 we don't have a real benchmark set, so we have to use the raw Platinum Pedigree dipcall calls.
            "chm13": os.path.join(TRUTH_DIR, wildcards["reference"], wildcards["sample"], wildcards["sample"] + ".dip.vcf.gz"),
            # On CHM13 we don't have a real benchmark set, so we have to use the raw Platinum Pedigree dipcall calls.
            "chm13v1": os.path.join(TRUTH_DIR, wildcards["reference"], wildcards["sample"], wildcards["sample"] + ".dip.vcf.gz"),
            # On GRCh38 there's a GIAB truth set
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        }[wildcards["reference"]]
    else:
        raise RuntimeError("Unsupported sample: " + wildcards["sample"])

def truth_vcf_index_url(wildcards):
    """
    Find the URL or local file for the variant calling truth VCF index, from reference and sample.
    """
    return truth_vcf_url(wildcards) + ".tbi"

def truth_bed_url(wildcards):
    """
    Find the URL or local file for the variant calling truth high confidence BED, from reference.

    If compressed, must end in ".gz".
    """

    if wildcards["sample"] == "HG002":
        # For HG002, these are available online directly.
        return {
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
        }[wildcards["reference"]]
    elif wildcards["sample"] == "HG001":
        return  {
            # TODO: On CHM13 we don't have Platinum Pedigree high-confidence regions, so
            # we need to just use the HG002 ones for other samples and hope they're close enough.
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            # TODO: On CHM13 we don't have Platinum Pedigree high-confidence regions, so
            # we need to just use the HG002 ones for other samples and hope they're close enough.
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            # On GRCh38 there's a GIAB truth set
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
        }[wildcards["reference"]]
    else:
        raise RuntimeError("Unsupported sample: " + wildcards["sample"])

def surjectable_gam(wildcards):
    """
    Find a GAM mapped to a graph built on a reference that we can use to
    surject to the reference we're interested in.
    """
    # TODO: Make similar redirects for the benchmarks!

    format_data = dict(wildcards)

    return "{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam".format(**format_data)

def graph_base(wildcards, force_all_refs=False):
    """
    Find the base name for a collection of graph files from reference and either refgraph or all of refgraphbase, modifications, and clipping.

    For graphs ending in "-sampled" and sampling parameters, autodetects the right haplotype-sampled full graph name to use from tech and sample.

    For GraphAligner, selects an unchopped version of the graph based on mapper.

    If force_all_refs is set, will drop 'o' from the sampling flags, if present,
    to get the version of the graph that has all the references present in the
    base graph.
    """

    # For membership testing, we need a set of wildcard keys
    wc_keys = set(wildcards.keys())
    reference = wildcards["reference"]
    # TODO: This sucks
    if reference == "chm13v1":
        reference = "chm13"
    modifications = []

    if "refgraphbase" in wc_keys and "modifications" in wc_keys:
        # We already have the reference graph base (hprc-v1.1-mc) and clipping
        # (.d9) and chopping (.unchopped) and full (.full) if allowed cut
        # apart. Sampling info may also be there (we assume it already
        # describes what to sample for specifically).
        refgraphbase = wildcards["refgraphbase"]
        modifications.append(wildcards["modifications"])
        if "clipping" in wc_keys:
            modifications.append(wildcards["clipping"])
        if "full" in wc_keys:
            modifications.append(wildcards["full"])
        if "chopping" in wc_keys:
            modifications.append(wildcards["chopping"])
        if "sampling" in wc_keys:
            modifications.append(wildcards["sampling"])
    else:
        assert "refgraph" in wc_keys, f"No refgraph wildcard in: {wc_keys}"
        # We need to handle hprc-v1.1-mc and hprc-v1.1-mc.full and hprc-v1.1-mc-d9 and hprc-v2.prereease-mc-R2-d32 and hprc-v2.prereease-mc-R2-sampled10d.
        # Also probaby primary.
        # They need the reference inserted after the -mc. but before the other stuff, and -d32 needs to become .d32.
        # And -sampled10d needs to be expanded to say what sample and read set we are haplotype sampling from.

        # TODO: If it's not -mc, where would the reference go?
        
        refgraph = wildcards["refgraph"]

        parts = refgraph.split("-")
        last = parts[-1]
        if re.fullmatch("unchopped", last):
            # We have a choppedness modifier, which gets a dot.
            modifications.append("." + last)
            parts.pop()

        last = parts[-1]
        if re.fullmatch("sampled[0-9]+d?o?", last):
            # We have a generic haplotype sampling flag.
            # Autodetect the right haplotype-sampled graph to use.

            # We need to convert this into the graph sampled for the right sample.
            # Which means we need some info about the sample we are working on.
            assert "tech" in wc_keys, "No tech known for haplotype sampling graph " + wildcards["refgraph"]
            assert "sample" in wc_keys, "No sample known for haplotype sampling graph " + wildcards["refgraph"]
           
            # TODO: We always sample for the full real trimmed-if-R10 version
            # of whatever reads we're going to map, so we can consistently use
            # one graph. Note r10y2025 doesn't get trimmed.
            sampling_trimmedness = ".trimmed" if wildcards["tech"] == "r10" else ""

            if force_all_refs and last.endswith("o"):
                # We need to get the version of the graph sampled with all
                # references instead.
                last = last[:-1]

            modifications.append(f"-{last}_fragmentlinked-for-real-{wildcards['tech']}-{wildcards['sample']}{sampling_trimmedness}-full")
            parts.pop()
        elif re.fullmatch("d[0-9]+", last):
            # We have a clipping modifier, which gets a dot.
            modifications.append("." + last)
            parts.pop()
        
        while len(parts) > 3:
            # We have more than just the 3-tuple of name, version, algorithm. Take the last thing into modifications.
            modifications.append("-" + parts[-1])
            parts.pop()
        
        last = parts[-1]
        if last.endswith(".full"):
            # Move a .full over the reference
            modifications.append(".full")
            last = last[:-5]
            parts[-1] = last
        

        # Now we have all the modifications. Flip them around the right way.
        modifications.reverse()
        
        # The first 3 or fewer parts are the graph base name.
        refgraphbase = "-".join(parts)

    result = os.path.join(GRAPHS_DIR, refgraphbase + "-" + reference + "".join(modifications))
    return result

def gbz(wildcards):
    """
    Find a graph GBZ file from reference.
    """
    return graph_base(wildcards) + ".gbz"

def all_refs_gbz(wildcards):
    """
    Find a graph GBZ file with all linear references from reference.
    """
    return graph_base(wildcards, force_all_refs=True) + ".gbz"

def minimizer_k(wildcards):
    """
    Find the minimizer kmer size from mapper.
    """
    if wildcards["mapper"].startswith("giraffe"):
        # Looks like "giraffe-k31.w50.W-lr-default-noflags".
        # So get second part on - and first part of that on . and number-ify it after the k.
        return int(wildcards["mapper"].split("-")[1].split(".")[0][1:])
    else:
        mode = minimap_derivative_mode(wildcards)
        match mode:
            # See minimap2 man page
            case "map-ont":
                return 15
            case "map-pb":
                return 19
            case "sr":
                return 21
        raise RuntimeError("Unimplemented mode: " + mode)

def dist_indexed_graph(wildcards):
    """
    Find a GBZ and its dist index from reference.
    """
    base = graph_base(wildcards)
    return {
        "gbz": gbz(wildcards),
        "dist": base + ".dist"
    }

def indexed_graph(wildcards):
    """
    Find an indexed graph and all its indexes from reference and minparams.

    Also checks vgversion to see if we need a no-zipcodes index for old vg.
    """
    base = graph_base(wildcards)
    indexes = {
        "minfile": base + "." + wildcards["minparams"] + ".path.min",
        "zipfile": base + "." + wildcards["minparams"] + ".path.zipcodes"
    }
    indexes.update(dist_indexed_graph(wildcards))
    return indexes

def base_fastq_gz(wildcards):
    """
    Find a full compressed FASTQ for a real sample, based on realness, sample,
    tech, and trimmedness.

    If an untrimmed version exists and the trimmed version does not, returns
    the name of the trimmed version to make.
    """
    import glob
    full_gz_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*{trimmedness}.f*q.gz".format(**wildcards))
    results = glob.glob(full_gz_pattern)
    if wildcards["trimmedness"] != ".trimmed":
        # Don't match trimmed files when not trimmed.
        results = [r for r in results if ".trimmed" not in r]
    # Skip any subset files
    subset_regex = re.compile("([0-9]+[km]?|full)")
    results = [r for r in results if not subset_regex.fullmatch(r.split(".")[-3])]
    # Skip any chunk files
    chunk_regex = re.compile("chunk[0-9]*")
    results = [r for r in results if not chunk_regex.fullmatch(r.split(".")[-3])]
    if len(results) == 0:
        # Can't find it
        if wildcards["trimmedness"] == ".trimmed":
            # Look for an untrimmed version
            untrimmed_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*.f*q.gz".format(**wildcards))
            results = glob.glob(untrimmed_pattern)
            # Skip any subset files
            results = [r for r in results if not subset_regex.fullmatch(r.split(".")[-3])]
            if len(results) == 1:
                # We can use this untrimmed one to make a trimmed one
                without_gz = os.path.splitext(results[0])[0]
                without_fq, fq_ext = os.path.splitext(without_gz)
                trimmed_base = without_fq + ".trimmed" + fq_ext + ".gz"
                return trimmed_base
        raise FileNotFoundError(f"No  non-chunk, non-subset files found matching {full_gz_pattern}")
    elif len(results) > 1:
        raise RuntimeError("Multiple non-chunk, non-subset files matched " + full_gz_pattern)
    return results[0]

def fastq_finder(wildcards, compressed = False):
    """
    Find a FASTQ or compressed FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """

    extra_ext = ".gz" if compressed else ""

    #If we chunked the reads, we need to add it. Otherwise there isn't a chunk
    readchunk = wildcards.get("readchunk", "")

    import glob
    fastq_by_sample_pattern = os.path.join(READS_DIR, ("{realness}/{tech}/{sample}/*{sample}*{trimmedness}[._-]{subset}" + readchunk+ ".f*q" + extra_ext).format(**wildcards))
    results = glob.glob(fastq_by_sample_pattern)
    if wildcards["trimmedness"] != ".trimmed":
        # Don't match trimmed files when not trimmed.
        results = [r for r in results if ".trimmed" not in r]
    if len(results) == 0:
        if wildcards["realness"] == "real":
            # Make sure there's a full .fq.gz to extract from (i.e. this doesn't raise)
            full_file = base_fastq_gz(wildcards)
            # And compute the subset name
            without_gz = os.path.splitext(full_file)[0]
            without_fq = os.path.splitext(without_gz)[0]
            return without_fq + (".{subset}" + readchunk+ ".fq" + extra_ext).format(**wildcards)
        elif wildcards["realness"] == "sim":
            # Assume we can get this FASTQ.
            # For simulated reads we assume the right subset GAM is there. We
            # don't want to deal with the 1k/1000 difference here.
            return os.path.join(READS_DIR, ("{realness}/{tech}/{sample}/{sample}-{realness}-{tech}{trimmedness}-{subset}" + readchunk+ ".fq" + extra_ext).format(**wildcards))
        else:
            raise FileNotFoundError(f"No files found matching {fastq_by_sample_pattern}")
    elif len(results) > 1:
        raise AmbiguousRuleException("Multiple files matched " + fastq_by_sample_pattern)
    return results[0]

def fastq(wildcards):
    """
    Find a FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    return fastq_finder(wildcards, compressed=False)

def fastq_gz(wildcards):
    """
    Find a compressed FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    return fastq_finder(wildcards, compressed=True)

def mapper_stages(wildcards):
    """
    Find the list of mapping stages from mapper.
    """

    if wildcards["mapper"].startswith("giraffe"):
        parts = wildcards["mapper"].split("-")
        assert len(parts) >= 3
        # Should be giraffe, then the minimizer parameters, then the parameter preset.
        if parts[3] == "default":
            # Default mapping preset is non-chaining
            return NON_CHAINING_STAGES
        else:
            return STAGES
    else:
        return []


def all_experiment_conditions(expname, filter_function=None, debug=False):
    """
    Yield dictionaries of all conditions for the given experiment.
    
    The config file should have a dict in "experiments", of which the given
    expname should be a key. The value is the experiment dict.

    The experiment dict should have a "control" dict, listing names and values
    of variables to keep constant.

    The experiment dict should have a "vary" dict, listing names and values
    lists of variables to vary. All combinations will be generated.

    The experiment dict should have a "constrain" list. Each item in the list
    is a "pass", which is a list of constraints. Each item in the pass is a
    dict of variable names and values (or lists of values). A condition must
    match *at least* one of these dicts on *all* values in the dict in order to
    survive the pass. And it must survive all passes in order to be run.

    If optional variables (defined internally) don't appear in the experiment
    "vary" or "control" sections, they will be given an empty-string value.

    If filter_function is provided, only yields conditions that the filter
    function is true for.

    Yields variable name to value dicts for all passing conditions for the
    given experiment.
    """

    if "experiments" not in config:
        raise RuntimeError(f"No experiments section in configuration; cannot run experiment {expname}")
    all_experiments = config["experiments"]
    
    if expname not in all_experiments:
        raise RuntimeError(f"Experiment {expname} not in configuration")
    exp_dict = all_experiments[expname]

    # Make a base dict of all controlled variables.
    base_condition = exp_dict.get("control", {})

    to_vary = exp_dict.get("vary", {})

    constraint_passes = exp_dict.get("constrain", [])

    total_conditions = 0
    for condition in augmented_with_all(base_condition, to_vary):
        # For each combination of independent variables on top of the base condition

        # Fill in optional experiment variables with empty strings
        for optional_key in ["trimmedness", "callparams"]:
            if optional_key not in condition:
                condition[optional_key] = ""

        # We need to see if this is a combination we want to do
        if matches_all_constraint_passes(condition, constraint_passes):
            if not filter_function or filter_function(condition):
                total_conditions += 1
                yield condition
            else:
                if debug:
                    print(f"Condition {condition} does not match requested filter function")
        else:
            if debug:
                print(f"Condition {condition} does not match a constraint in some pass")
    print(f"Experiment {expname} has {total_conditions} eligible conditions")
    

def augmented_with_each(base_dict, new_key, possible_values):
    """
    Yield copies of base_dict with each value from possible_values under new_key.
    """

    for value in sorted(possible_values):
        clone = dict(base_dict)
        clone[new_key] = value
        yield clone

def augmented_with_all(base_dict, keys_and_values):
    """
    Yield copies of base_dict augmented with all combinations of values from
    keys_and_values, under the corresponding keys.
    """

    if len(keys_and_values) == 0:
        # Base case: nothing to add
        yield base_dict
    else:
        # Break off one facet
        first_key = next(iter(keys_and_values.keys()))
        first_values = keys_and_values[first_key]
        rest = dict(keys_and_values)
        del rest[first_key]

        for with_first in augmented_with_each(base_dict, first_key, first_values):
            # Augment with this key
            for with_rest in augmented_with_all(with_first, rest):
                # And augment with the rest
                yield with_rest

def matches_constraint_value(query, value):
    """
    Returns True if query equals value, except if value is a list, query has to
    be in the list instead.
    """

    if isinstance(value, list):
        return query in value
    else:
        return query == value

def matches_constraint(condition, constraint, debug=False):
    """
    Returns True if all keys in constraint are in condition with the same
    values, or with values in the list in constraint.
    """
    for key, match in constraint.items():
        if key not in condition or not matches_constraint_value(condition[key], match):
            if debug:
                print(f"Condition {condition} mismatched constraint {constraint} on {key}")
            return False
    return True

def matches_any_constraint(condition, constraints):
    """
    Return True if, for some constraint dict, the condition dict matches all
    values in the constraint dict.
    """

    for constraint in constraints:
        if matches_constraint(condition, constraint):
            return True
    return False

def matches_all_constraint_passes(condition, passes):
    """
    Return True if the condfition matches some constraint in each pass in passes.
    """
    
    if len(passes) > 0 and not isinstance(passes[0], list) and isinstance(passes[0], dict):
        # Old style config where there's just one pass of constraints. Fix it up.
        passes = [passes]

    for constraints in passes:
        if not matches_any_constraint(condition, constraints):
            return False
    return True

def wildcards_to_condition(all_wildcards):
    """
    Filter down wildcards to just the condition parameters for the experiment
    in expname.
    
    Raises an error if any varied variable in the experiment cannot be
    determined (unless it's a variable that only matters for calling).
    """

    exp_dict = config.get("experiments", {}).get(all_wildcards["expname"], {})
    base_condition = exp_dict.get("control", {})
    to_vary = exp_dict.get("vary", {})
    
    # For membership testing, we need a set of wildcard keys
    wc_keys = set(all_wildcards.keys())

    condition = {}

    for var in base_condition.keys():
        if var in wc_keys:
            # Constant variables across the whole experiment don't need to
            # actually be in the condition name unless we actually have them
            # available.
            condition[var] = all_wildcards[var]

    for var in to_vary.keys():
        if len(all_wildcards.get(var, "")) != 0:
            condition[var] = all_wildcards[var]
        elif var == "mapper":
            #If we haven't defined a mapper, but we have the components of a giraffe mapper
            minimizer_params = all_wildcards.get("minparams", "k"+all_wildcards["k"]+".w"+all_wildcards["w"]+all_wildcards["weightedness"])
            condition[var] = "giraffe-" + minimizer_params + "-" + all_wildcards["preset"] + "-" + all_wildcards["vgversion"] + "-" + all_wildcards["vgflag"]
        elif var == "refgraph" and len(all_wildcards.get("refgraphbase", "")) != 0:
            #If we haven't defined refgraph, but we have the components
            #TODO: This isn't great but it works in the one case I need it to
            graph = all_wildcards.get("refgraphbase")
            if len(all_wildcards.get("sampling", "")) != 0:
                graph += "-"
                graph += all_wildcards.get("sampling")
            graph += all_wildcards.get("clipping","")
            graph += all_wildcards.get("full","")
            graph += all_wildcards.get("chopping","")

            condition[var] = graph
        elif var == "callparams" and len(all_wildcards.get("callparams", "")) == 0:
            # There's no callparams, but that's a calling-stage variable and we
            # might be asking mapping-stage questions.
            pass
        else:
            #Catch any case where it fails
            condition[var] = all_wildcards[var]



    return condition

def condition_name(wildcards):
    """
    Determine a human-readable condition name from expname and the experiment's variable values.
    """
 
    # Get what changes in the experiment
    exp_dict = config.get("experiments", {}).get(wildcards["expname"], {})
    to_vary = exp_dict.get("vary", {})

    # Get the condition dict in use here
    condition = wildcards_to_condition(wildcards)

    name_parts = []

    for varied_key in to_vary:

        #Don't include realness
        if varied_key == "realness":
            continue

        if varied_key == "callparams" and varied_key not in condition:
            # Skip calling params if we're not asking about a calling-stage thing
            continue

        # Look at the value we have for this varied variable
        condition_value = condition[varied_key]
        # And the other possible values that are used
        alternatives = set(to_vary.get(varied_key, []))

        if varied_key == "mapper" and "giraffe" in condition_value:
            # If we're working on a Giraffe mapper name, only compare against other Giraffe mapper names
            alternatives = {a for a in alternatives if "giraffe" in a}

        # Find all the name parts used in alternatives.
        parts_in_alternatives = [set(alternative.split("-")) for alternative in alternatives]
        # Find those used in all alternatives
        universal_alternative_parts = functools.reduce(lambda a, b: a & b, parts_in_alternatives)

        # Find all the name parts used in us
        condition_value_parts = condition_value.split("-")
       
        # Drop parts that are universal among applicable alternatives, keeping
        # only the parts that represent differences.
        interesting_parts = [p for p in condition_value_parts if p not in universal_alternative_parts]
        
        if varied_key == "mapper" and "giraffe" in condition_value:
            # Make sure to mark Giraffe conditions even when all of them have "giraffe" in them.
            interesting_parts = ["giraffe"] + interesting_parts

        # And add that to the name
        name_parts.append("-".join(interesting_parts))

    return ",".join(name_parts)

def all_experiment(wildcard_values, pattern, filter_function=None, empty_ok=False, debug=False):
    """
    Produce all values of pattern substituted with the wildcards and the experiment conditions' values, from expname.

    If provided, restricts to conditions passing the filter function.

    Throws an error if nothing is produced and empty_ok is not set.

    Needs to be used like:
        lambda w: all_experiment(w, "your pattern")
    """

    empty = True
    for condition in all_experiment_conditions(wildcard_values["expname"], filter_function=filter_function):
        merged = dict(wildcard_values)
        merged.update(condition)

        # TODO: Hackily fill in the optional {dot} which should be set if {category} is set.
        # There's no way in a Snakemake expandion to tack on a leader sequence.
        if "category" in merged and len(merged["category"]) > 0:
            merged["dot"] = "."
        else:
            merged["dot"] = ""
        if "category" not in merged:
            # Category might appear in templates but in experiments it is meant to be optional.
            merged["category"] = ""

        if debug:
            print(f"Evaluate {pattern} in {merged} from {wildcard_values} and {condition}")
        filename = pattern.format(**merged)
        yield filename
        empty = False
    if empty:
        if debug:
            print("Produced no values for " + pattern + " in experiment!")
        if not empty_ok:
            if debug:
                print("Failing!")
            raise RuntimeError("Produced no values for " + pattern + " in experiment!")

def has_stat_filter(stat_name):
    """
    Produce a filter function for conditions that might have the stat stat_name.

    Applies to stat files like:
    {root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv

    Use with all_experiment() when aggregating stats to compare, to avoid
    trying to agregate from conditions for which the stat cannot be measured.
    """

    def filter_function(condition):
        """
        Return True if the given condition dict should have the stat named stat_name.
        """

        if stat_name in {"correct", "accuracy", "wrong"}:
            # These stats only exist for conditions with a truth set (i.e. simulated ones)
            if condition["realness"] != "sim":
                return False

        if stat_name.startswith("indel_") or stat_name.startswith("snp_"):
            # These are calling stats, and we shoudl only do calling for real reads.
            if condition["realness"] != "real":
                return False

        if stat_name.startswith("time_used") or stat_name in ("mapping_speed", "chain_coverage"):
            # This is a Giraffe time used stat or mean thereof. We need to be a
            # Giraffe condition.
            if not condition["mapper"].startswith("giraffe"):
                return False

        return True

    return filter_function

def get_vg_flags(wildcard_flag):
    match wildcard_flag:
        case "gapExt":
            return "--do-gapless-extension"
        case "mqCap":
            return "--explored-cap"
        case downsample_number if downsample_number[0:10] == "downsample":
            return "--downsample-min " + downsample_number[10:]
        case "mapqscale":
            return "--mapq-score-scale 0.01"
        case scorescale_number if scorescale_number[0:10] == "scorescale":
            return "--mapq-score-scale " + scorescale_number[10:]
        case "moreseeds":
            return "--downsample-window-length 400"
        case "mqWindow":
            return "--mapq-score-scale 1 --mapq-score-window 150"
        case mcspb_number if mcspb_number[0:5] == "mcspb":
            return "--min-chain-score-per-base " + mcspb_number[5:]
        case mcspem_number if mcspem_number[0:6] == "mcspem":
            return "--min-chain-score-per-base 0.0 --min-chain-score-per-explored-minimizer " + mcspem_number[6:]
        case sp_number if sp_number[0:2] == "sp":
            return "--softclip-penalty " + sp_number[2:]
        case mmcs_number if mmcs_number[0:4] == "mmcs":
            return "--max-min-chain-score " + mmcs_number[4:]
        case "candidate1":
            return "--max-min-chain-score 60 --chain-score-threshold 150"
        case "noflags":
            return ""
        case unknown:
            #otherwise this is a hash and we get the flags from ParameterSearch
            return PARAM_SEARCH.hash_to_parameter_string(wildcard_flag)

def get_vg_version(wildcard_vgversion):
    if wildcard_vgversion == "default":
        return "vg"
    else:
        return "./vg_"+wildcard_vgversion

rule giraffe_real_reads:
    input:
        unpack(indexed_graph),
        fastq_gz=fastq_gz,
    output:
        # Giraffe can dump out pre-annotated reads at annotation range -1.
        gam="{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    log:"{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    benchmark: "{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=REAL_SLURM_EXTRA,
        full_cluster_nodes=0
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)
        pairing_flag="-i" if wildcards.preset == "default" else ""
        zipcodes_flag=f"-z {input.zipfile}" if "zipfile" in dict(input).keys() else ""

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress -Z {input.gbz} -d {input.dist} -m {input.minfile} -f {input.fastq_gz} " + zipcodes_flag + " " + flags + " " + pairing_flag + " >{output.gam} 2>{log}")

rule surject_gam:
    input:
        gbz=gbz,
        # We surject onto all target reference paths, even those we can't call
        # on, like Y in CHM13 (due to different Ys being used in different
        # graphs), and the _random contigs in GRCh38. Otherwise mappers mapping
        # against GRCh38 with its _random contigs have an advantage when
        # calling there, because Giraffe will stick all the _random reads in
        # the main genome.
        reference_dict=reference_dict,
        gam=surjectable_gam
    output:
        # Don't keep this around because we're going to keep the sorted version and use that for most things.
        bam=temp("{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam")
    wildcard_constraints:
        mapper="(giraffe.*|graphaligner-.*)"
    params:
        # The default preset is paired
        paired_flag=lambda w: "-i" if re.match("giraffe-[^-]*-default-[^-]*-[^-]*", w["mapper"]) else ""
    threads: 40
    resources:
        mem_mb=lambda w: (600000 / 64 * 40) if w["tech"] in ("r10", "r10y2025") else (150000 / 64 * 40),
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg surject -F {input.reference_dict} -x {input.gbz} -t {threads} --bam-output --sample {wildcards.sample} --read-group \"ID:1 LB:lib1 SM:{wildcards.sample} PL:{wildcards.tech} PU:unit1\" --prune-low-cplx {params.paired_flag} {input.gam} > {output.bam}"

rule sort_bam:
    input:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam",
        bai="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam.bai"
    threads: 16
    resources:
        mem_mb=16000,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        with tempfile.TemporaryDirectory() as sort_scratch:
            shell("samtools sort -T " + os.path.join(sort_scratch, "scratch") + " --threads {threads} {input.bam} -O BAM > {output.bam} && samtools index -b {output.bam} {output.bai}")

rule call_variants_dv:
    input:
        sorted_bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam",
        sorted_bam_index="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam.bai",
        reference_path_list_callable=reference_path_list_callable,
        calling_reference_fasta=calling_reference_fasta,
        calling_reference_fasta_index=calling_reference_fasta_index,
        calling_reference_restrict_bed=calling_reference_restrict_bed,
        calling_reference_par_bed=calling_reference_par_bed,
        truth_vcf=lambda w: remote_or_local(truth_vcf_url(w)),
        truth_vcf_index=lambda w: remote_or_local(truth_vcf_index_url(w)),
        truth_bed=lambda w: remote_or_local(truth_bed_url(w))
    output:
        wdl_output_file="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.json",
        # TODO: make this temp so we can delete it?
        wdl_output_directory=directory("{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.wdlrun"),
        vcf="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.vcf.gz",
        vcf_index="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.vcf.gz.tbi",
        happy_evaluation_archive="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.happy_results.tar.gz"
    log:
        logfile="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.log"
    params:
        wdl_input_file="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.input.json",
        job_store="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.jobstore",
        batch_logs="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.batchlogs",
        reference_prefix=reference_prefix,
        haploid_contigs=haploid_contigs,
        wdl_cache=wdl_cache,
        model_files=model_files
    threads: 2
    resources:
        mem_mb=4000,
        # Since Toil tasks need to schedule through the cluster again, this
        # might need to wait a long time.
        # TODO: Run single machine so Slurm can't plan to run our jobs *after*
        # the leader finishes?
        runtime=8640,
        slurm_partition=choose_partition(8640)
    run:
        import json
        if ".olddv" in wildcards.callparams:
            dv_docker = "gcr.io/deepvariant-docker/deepvariant:head756846963"
        elif ".newdv" in wildcards.callparams:
            dv_docker = "google/deepvariant:CL782981885"
        else:
            # Default to old DeepVariant until https://ucsc-gi.slack.com/archives/C01D0M09G5D/p1753482919257519?thread_ts=1753201495.785309&cid=C01D0M09G5D is fixed.
            dv_docker = "gcr.io/deepvariant-docker/deepvariant:head756846963"
        wf_source = "#workflow/github.com/vgteam/vg_wdl/DeepVariant:lr-giraffe-paper"
        wf_inputs = {
            "DeepVariant.MERGED_BAM_FILE": input.sorted_bam,
            "DeepVariant.MERGED_BAM_FILE_INDEX": input.sorted_bam_index,
            "DeepVariant.SAMPLE_NAME": wildcards.sample,
            "DeepVariant.PATH_LIST_FILE": input.reference_path_list_callable,
            "DeepVariant.REFERENCE_PREFIX": params.reference_prefix,
            "DeepVariant.REFERENCE_PREFIX_ON_BAM": True,
            "DeepVariant.REFERENCE_FILE": input.calling_reference_fasta,
            # TODO: Should we not left align for minimap2 or GraphAligner?
            "DeepVariant.LEFTALIGN_BAM": True,
            # Indel realignment tools aren't available for long reads.
            # Also, we can't realign indels on our CHM13 non-2.0 mapping
            # reference when using our CHM13v2.0 calling reference.
            # TODO: Turn indel realignment back on if we ever we manage to use
            # CHM13v2.0. Which we probably won't because it uses HG002 Y and
            # we hold that out of the eval graphs.
            "DeepVariant.REALIGN_INDELS": False,
            "DeepVariant.DV_MODEL_TYPE": {"hifi": "PACBIO", "r10": "ONT_R104", "r10y2025": "ONT_R104", "illumina": "WGS"}[wildcards.tech],
            # We don't send MIN_MAPQ and let the model example info define it.
            # We also let the model control DV_KEEP_LEGACY_AC, which we made default to not passed in the workflow.
            # Read normalization should work for long reads as of gcr.io/deepvariant-docker/deepvariant:head756846963 (which is post-1.9)
            # But we also don't send DV_NORM_READS because the model defines it as of google/deepvariant:CL782981885
            # Work around <https://github.com/google/deepvariant/issues/989>
            "DeepVariant.OTHER_MAKEEXAMPLES_ARG": "--small_model_call_multiallelics=false" if wildcards.tech == "r10y2025" else None,
            "DeepVariant.DV_GPU_DOCKER": dv_docker + "-gpu",
            "DeepVariant.DV_NO_GPU_DOCKER": dv_docker,
            "DeepVariant.TRUTH_VCF": to_local(input.truth_vcf),
            "DeepVariant.TRUTH_VCF_INDEX": to_local(input.truth_vcf_index),
            "DeepVariant.EVALUATION_REGIONS_BED": to_local(input.truth_bed),
            "DeepVariant.RESTRICT_REGIONS_BED": input.calling_reference_restrict_bed,
            "DeepVariant.PAR_REGIONS_BED_FILE": input.calling_reference_par_bed,
            "DeepVariant.HAPLOID_CONTIGS": params.haploid_contigs,
            # We just need hap.py; we don't need a stand-alone vcfeval run.
            "DeepVariant.RUN_STANDALONE_VCFEVAL": False,
            # We also don't need BAMs.
            "DeepVariant.OUTPUT_SINGLE_BAM": False,
            "DeepVariant.OUTPUT_CALLING_BAMS": False,
            "DeepVariant.CALL_CORES": 8 * 4,
            "DeepVariant.CALL_MEM": 50 * 4,
            "DeepVariant.MAKE_EXAMPLES_MEM": 50
        }
        if dv_docker == "gcr.io/deepvariant-docker/deepvariant:head756846963":
            # We can't use a model example_info file, so we still need to set some flags.
            # TODO: Should this be different for the non-Giraffe controls????
            wf_inputs["DeepVariant.MIN_MAPQ"] = 0
            wf_inputs["DeepVariant.DV_NORM_READS"] = True
            wf_inputs["DeepVariant.DV_KEEP_LEGACY_AC"] = False
        if params.model_files:
            # Use a model that's not built in.
            wf_inputs["DeepVariant.DV_MODEL_FILES"] = params.model_files
        json.dump(wf_inputs, open(params["wdl_input_file"], "w"))
        shell("rm -Rf {params.job_store}")
        # Run and keep the first manageable amount of logs not sent to the log
        # file in case we can't start. Don't stop when we hit the log limit.
        # See https://superuser.com/a/1531706
        toil_command = "MINIWDL__CALL_CACHE__GET=true MINIWDL__CALL_CACHE__PUT={FILL_WDL_CACHE} MINIWDL__CALL_CACHE__DIR={params.wdl_cache} toil-wdl-runner '" + wf_source + "' {params.wdl_input_file} --clean=onSuccess --jobStore file:{params.job_store} --batchLogsDir {params.batch_logs} --wdlOutputDirectory {output.wdl_output_directory} --wdlOutputFile {output.wdl_output_file} --batchSystem slurm --slurmTime 11:59:59 --disableProgress --caching=False --logFile={log.logfile} 2>&1 | (head -c1000000; cat >/dev/null)"
        print("Running Toil: " + toil_command)
        shell(toil_command)

        wdl_result=json.load(open(output.wdl_output_file))
        shell("cp " + wdl_result["DeepVariant.output_vcf"] + " {output.vcf}")
        shell("cp " + wdl_result["DeepVariant.output_vcf_index"] + " {output.vcf_index}")
        shell("cp " + wdl_result["DeepVariant.output_happy_evaluation_archive"] + " {output.happy_evaluation_archive}")
        
rule extract_happy_summary:
    input:
        happy_evaluation_archive="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.happy_results.tar.gz"
    output:
        happy_evaluation_summary="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.eval.summary.csv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "tar -xOf {input.happy_evaluation_archive} happy_results/eval.summary.csv >{output.happy_evaluation_summary}"

rule stat_from_happy_summary:
    input:
        happy_evaluation_summary="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.eval.summary.csv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_{colname}.tsv"
    wildcard_constraints:
        vartype="(snp|indel)",
        colname="(f1|precision|recall|fn|fp|tp)"
    params:
        colnum=lambda w: {"f1": 14, "precision": 12, "recall": 11, "fn": 5, "fp": 7, "tp":4}[w["colname"]]
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        shell("cat {input.happy_evaluation_summary} | grep '^" + wildcards["vartype"].upper() + ",PASS' | cut -f{params.colnum} -d',' >{output.tsv}")

rule total_errors_from_fp_and_fn:
    input:
        snp_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.snp_fn.tsv",
        snp_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.snp_fp.tsv",
        indel_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.indel_fn.tsv",
        indel_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.indel_fp.tsv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.total_errors.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.snp_fn} {input.snp_fp} {input.indel_fn} {input.indel_fp} | awk '{{sum += $1 }} END {{ print sum }}' >{output.tsv}"

rule vartype_errors_from_fp_and_fn:
    input:
        vartype_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_fn.tsv",
        vartype_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_fp.tsv",
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_errors.tsv"
    wildcard_constraints:
        vartype="(snp|indel)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.vartype_fn} {input.vartype_fp} | awk '{{sum += $1 }} END {{ print sum }}' >{output.tsv}"

# Some experiment stats can come straight from stats for the individual conditions
rule condition_experiment_stat:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{conditionstat}.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{conditionstat}.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        conditionstat="((overall_fraction_)?(wrong|correct|eligible|(positive_)?(unmapped|mismapped))|accuracy|(snp|indel)_(f1|precision|recall|fn|fp)|(snp|indel|total)_errors|[a-zA-Z0-9_]*.total)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_calling_stat_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{statname}.tsv", filter_function=has_stat_filter(w["statname"]))
    output:
        table="{root}/experiments/{expname}/results/{statname}.tsv"
    wildcard_constraints:
        statname="((snp|indel)_(f1|precision|recall|fn|fp)|(snp|indel|total)_errors)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "cat {input} >{output.table}"

rule experiment_total_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/total_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/total_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Total Point Variant Errors' --y_label 'Total Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_snp_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/snp_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/snp_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} SNP Errors' --y_label 'Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_indel_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/indel_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/indel_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Indel Errors' --y_label 'Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"
