# Variant calling comparisons across `vg giraffe` versions

This repository is a fork of [`long-read-giraffe-experiments`](https://github.com/vgteam/long-read-giraffe-experiments).
To make this one easier to use, it has strictly less functionality and will only let you compare 
variant calling results between `vg giraffe` versions. If you only need read alignment comparisons,
please check out the [`giraffe-parameter-search`](https://github.com/vgteam/giraffe-parameter-search) repo.

**Set up this repository in `/private/groups`**. The pipeline will create very large files
(e.g. read alignments) which will not fit within `/private/home` directories.

## Usage

The general flow for using this repository is:

1. Get static `vg` binaries for the versions you want to compare. For each, navigate to inside your `vg` installation folder.
    1. `git checkout` the branch (e.g. `git checkout mybranch`) or commit hash (e.g. `git checkout e2547e979d86a446e8151c5489850b6176a82c0c`)
    2. Note the first six characters of the commit hash. If you don't know them, run `git rev-parse HEAD`.
    3. Make a static binary with:
        ````bash
        make -j32 static
        strip -d bin/vg
        ```
    4. Copy the binary file to inside this directory, where `COMMIT` is replaced by the six-character commit hash:
        ```
        # For example; replace with where this directory is
        VC_TEST_REPO=/private/groups/patenlab/$USER/giraffe-variant-calling-tests
        cp bin/vg $VC_TEST_REPO/vg_COMMIT
        ```
2. Update the `compare_variant_calling` experiment by replacing the commit hashes
  (`hifi-COMMIT-noflags` & `r10-COMMIT-noflags`) with the six-character hashes of your binaries.
  There should be one `hifi` and one `r10` version for each of the commits you're testing.
  Make sure the `constrain` section has the `hifi` versions for `hifi` reads, and visa versa.

3. Run the Snakefile (invoke from right inside this folder)
    ```
    (umask 002 && snakemake -j128 --rerun-incomplete --use-singularity --singularity-args "-B /private" --latency-wait 120 --executor slurm --keep-going ./output/experiments/compare_variant_calling/plots/{snp,indel,total}_errors.png)
    ```

## Dependencies

The following conda environment was sufficient to run vg version comparison variant calling tests on the Phoenix cluster:
```
conda create -n long-read-exp -c conda-forge -c bioconda snakemake=9.13.7 \
    toil snakemake-executor-plugin-slurm snakemake-storage-plugin-http singularity \
    bidict matplotlib meryl 'minimap2>=2.28' seqkit
```