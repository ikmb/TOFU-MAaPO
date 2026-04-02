# Running on HPC Systems Without Internet Access on Compute Nodes

Some HPC systems do not allow internet access from regular compute nodes. This matters for TOFU-MAaPO because:

- the Nextflow main process needs internet access to download required containers if they are not cached yet
- download steps in the pipeline need internet access while they run

>**Note**: This page covers runtime behavior on HPC systems (with restricted internet access). For general software installation, database setup, and the baseline structure of `custom.config`, see [Installation](./installation.md).

## Before you start

Check these points first:

1. Start TOFU-MAaPO from a node that has internet access.
   This is usually a login node, head node, or dedicated submission node.
2. Keep the Nextflow main process alive in `tmux` or `screen`.
3. Make sure your compute nodes can access the same shared filesystem as the launch node.
4. Make sure your container cache is on shared storage, or pre-download the containers before submitting compute jobs.

If the container cache is only stored on the launch node's local disk, compute nodes will not be able to reuse those containers. The default example in [`conf/custom.config`](../conf/custom.config) stores the cache under `${launchDir}/singularity_cache`, so on HPC systems you may want to change this to a shared filesystem path.






## Choose your setup
>**Note**: This documentation uses SLURM as the example workload manager, but other solutions such as [SGE and others](https://docs.seqera.io/nextflow/executor) are supported and can be used analogously.


Use the scenario that matches your cluster:

1. The login or submission node has internet access, and your cluster provides internet-enabled SLURM nodes or a dedicated partition.
   Route `local_download` jobs to that partition or node group.
2. Only the login or submission node has internet access.
   Run `local_download` jobs with `executor = 'local'` on the same node as the Nextflow main process.
3. You are unsure which nodes have internet access.
   Check your site documentation or ask your HPC admins which SLURM partition, constraint, or nodes should be used for internet-enabled jobs.



## Recommended setup

The best-practice setup is:

1. Start the pipeline from an internet-enabled login, head, or submission node.
2. Run it inside `tmux` or `screen`.
3. Let Nextflow submit the heavy compute jobs to SLURM from there.

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

## Configuration for HPCs with limited internet access

TOFU-MAaPO defines the `local_download` label for download processes such as:

- `download_sra`
- `download_files`

This means you do not need to modify the pipeline code itself. You only need to map that label to the correct SLURM partition, constraint, or local executor in your configuration.

Add or extend the `process` block in your custom Nextflow configuration as recommended in the following examples:

### Scenario 1: Dedicated internet-enabled SLURM partition or nodes

Some clusters provide a dedicated SLURM partition, queue, or node group for jobs that require internet access. In that case, route only the internet-dependent TOFU-MAaPO jobs there.

In Nextflow, the setting is called `queue`, even when your cluster documentation calls it a SLURM partition.

### Example: dedicated partition

```groovy
process {
  executor = 'slurm'
  queue = 'compute'

  withLabel: 'local_download' {
    executor = 'slurm'
    queue = 'internet'
  }
}
```

In this example:

- most jobs run in the default `compute` partition
- jobs labeled `local_download` run in the `internet` partition

### Example: dedicated nodes via `clusterOptions`

If your site uses a specific node list instead of a dedicated partition:

```groovy
process {
  executor = 'slurm'
  queue = 'compute'

  withLabel: 'local_download' {
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--nodelist=inet-node01'
  }
}
```

### Example: dedicated nodes via SLURM constraint

If your site uses features or constraints:

```groovy
process {
  executor = 'slurm'
  queue = 'compute'

  withLabel: 'local_download' {
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--constraint=internet'
  }
}
```

Use whichever SLURM option matches your local infrastructure.

### Scenario 2: Only the login or submission node has internet access

If no compute nodes can access the internet, but the node running the Nextflow main process can, then run internet-dependent jobs locally on that same node:

```groovy
process {
  executor = 'slurm'
  queue = 'compute'

  withLabel: 'local_download' {
    executor = 'local'
    maxForks = 10
  }
}
```

`maxForks` limits how many download jobs run in parallel on that node. If your admins only allow light work on the login or submission node, lower this value. If your site does not allow such jobs there at all, ask for a dedicated internet-enabled submission node or partition.
