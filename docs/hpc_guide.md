# Run TOFU-MAaPO on HPC systems

Check these points first:
1. Create a configuration file for TOFU-MAaPO for your computing environment as [as explained here](./installation.md#configuration).
2. Make sure your compute nodes can access the same shared filesystem as the launch node.
3. Make sure your container cache is set and writable. If needed, pre-download the containers before submitting compute jobs by performing a [test run](#perform-a-test-run).
4. Make sure the download jobs are run on nodes with internet access. [See our guide for HPCs with restricted internet access](/docs/hpc_internet_access.md#running-on-hpc-systems-with-restricted-internet-access).
5. Start TOFU-MAaPO from a node that has internet access.
   This is usually a login node, head node, or dedicated submission node. If no node has access, perform a [test run](#perform-a-test-run) on a system with internet access and migrate all downloaded data.
6. Keep the Nextflow main process alive in `tmux` or `screen` as [explained here](#recommendation-to-start-tofu-maapo).


>***Note***: If the container cache is only stored on the launch node's local disk, compute nodes will not be able to reuse those containers. The default example in [`conf/custom.config`](../conf/custom.config) stores the cache under `${launchDir}/singularity_cache`. On HPC systems you may want to change this to a shared filesystem path.

# Perform a test run
We recommend to perform one preparation run in an environment that has internet access, such as a login node, submission node, or another local system. 
To do this, use the same TOFU-MAaPO command you plan to run later, but add `-stub-run`.

The `-stub-run` parameter is useful to check whether the pipeline launches correctly on your system without calculating data, but it will also:
- pull required software containers into the shared container cache
- download required databases to their configured persistent locations
- fetch and download SRA input data to be reused in later runs

>**Note**: `-stub-run` will produce dummy results, we therefore recommend to use a distinct `--outdir` for this run.

## Example preparation run

Use the same flags you would use for the real analysis, including any database update flags and SRA input parameters, and add `-stub-run`. Such a preparation command would look like this:

```bash
nextflow run ikmb/tofu-maapo \
  -profile custom \
  -c /path/to/tofu.config \
  --sra 'SRR12345678,SRR12345679' \
  --apikey **YOUR_NCBI_API_KEY** \
  --publish_rawreads \
  --humann \
  --metaphlan \
  --assembly \
  --updatehumann \
  --updatemetaphlan \
  --updategtdbtk \
  --metaphlan_db /shared/db/metaphlan \
  --humann_db /shared/db/humann \
  --gtdbtk_reference /shared/db/gtdbtk \
  --outdir results_stub \
  -stub-run
```



# Recommendation to start TOFU-MAaPO

The best-practice setup is:

1. Start the pipeline from an internet-enabled login, head, or submission node.
2. Run it inside `tmux` or `screen`.
3. Let Nextflow submit the compute jobs to SLURM from there.

Example with `tmux`:

```bash
tmux new -s tofu
nextflow run main.nf -profile custom -c /path/to/tofu.config
```

Detach from `tmux` with `Ctrl-b` then `d`, and reconnect later with:

```bash
tmux attach -t tofu
```

Example with `screen`:

```bash
screen -S tofu
nextflow run main.nf -profile custom -c /path/to/tofu.config
```

Detach from `screen` with `Ctrl-a` then `d`, and reconnect later with:

```bash
screen -r tofu
```

This keeps the pipeline alive even if your SSH connection drops.

If your site does not allow long-running sessions on the login node, use the same approach on a dedicated submission node that has internet access.