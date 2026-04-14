# Running on HPC Systems with restricted internet access
>**Note**: This page covers runtime behavior on HPC systems (with restricted internet access). For general software installation, database setup, and the baseline structure of `custom.config`, see [Installation](/docs/installation.md).

Some HPC systems do not allow internet access from regular compute nodes. This matters for TOFU-MAaPO because:

- the Nextflow main process needs internet access to download required containers if they are not cached yet
- download steps in the pipeline need internet access while they run

Here, we want to give some guide lines how to configure and test TOFU-MAaPO in these environments.

Table of content:
- [Running on HPC Systems with restricted internet access](#running-on-hpc-systems-with-restricted-internet-access)
- [How to configure TOFU-MAaPO for your HPC with network restrictions](#how-to-configure-tofu-maapo-for-your-hpc-with-network-restrictions)
  - [Configuration for HPCs with limited internet access](#configuration-for-hpcs-with-limited-internet-access)
  - [Scenario 1: Dedicated internet-enabled SLURM partition or nodes](#scenario-1-dedicated-internet-enabled-slurm-partition-or-nodes)
    - [Example: dedicated partition](#example-dedicated-partition)
    - [Example: dedicated nodes via `clusterOptions`](#example-dedicated-nodes-via-clusteroptions)
    - [Example: dedicated nodes via SLURM constraint](#example-dedicated-nodes-via-slurm-constraint)
  - [Scenario 2: Only the login or submission node has internet access](#scenario-2-only-the-login-or-submission-node-has-internet-access)

# How to configure TOFU-MAaPO for your HPC with network restrictions
>**Note**: This documentation uses SLURM as the example workload manager, but other solutions such as [SGE and others](https://docs.seqera.io/nextflow/executor) are supported and can be used analogously.


Use the scenario that matches your cluster:

1. The login or submission node has internet access, and your cluster provides internet-enabled SLURM nodes or a dedicated partition.
   Route `local_download` jobs to that partition or node group [as explained here](#scenario-1-dedicated-internet-enabled-slurm-partition-or-nodes).
2. Only the login or submission node has internet access.
   Run `local_download` jobs with `executor = 'local'` on the same node as the Nextflow main process [as explained here](#scenario-2-only-the-login-or-submission-node-has-internet-access).
3. You are unsure which nodes have internet access.
   Check your site documentation or ask your HPC admins which SLURM partition, constraint, or nodes should be used for internet-enabled jobs.

>**Note**: You can use [`-stub-run`](/docs/hpc_guide.md#perform-a-test-run) to test your pipeline configuration and/or download software containers and databases.

## Configuration for HPCs with limited internet access

TOFU-MAaPO defines the `local_download` label for download processes such as:

- `download_sra`
- `download_files`

This means you do not need to modify the pipeline code itself. You only need to map that label to the correct SLURM partition, constraint, or local executor in your configuration file.

Add or extend the `process` block in your custom Nextflow configuration as recommended in the following examples:

## Scenario 1: Dedicated internet-enabled SLURM partition or nodes

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

## Scenario 2: Only the login or submission node has internet access

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
