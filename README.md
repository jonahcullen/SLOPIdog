# SLOPIdog

**Download the container**

Due to the size of the included canine reference genomes and index files (and depending on your internet speed) this may take some time.

```
singularity pull docker://hexive/slopi:2.0
```

To run the pipeline, update:
config.yaml: 
    all relevant file locations, as well as the desired contigs in their desired order.

submit_slopi.slurm: 
    the supercomputer partitions for your cluster and user email address, as well as what folders to bind to singularity

slurm.go_wags/slurm-submit.py: 
    partitions for your supercomputer cluster and account name