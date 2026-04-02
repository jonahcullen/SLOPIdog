# SLOPIdog

**Download the container**

Due to the size of the included canine reference genomes and index files (and depending on your internet speed) this may take some time.

```
singularity pull docker://hexive/slopi:2.0
```

**Required inputs**
A .txt list of the absolute paths of the lowpass bams to impute, one file per line

A linkage map for your reference genome (see linkage map found at /home/refgen/dog/canfam4/canFam4.linkage.map in the container for a working example)

A reference genome (.fa) with a .fai index

A .vcf.gz phased reference panel with a tabix (.tbi) index

If using canFam4, the linkage map and .fa are optional

**To run the pipeline, update the fields in:**

config.yaml: 
    all relevant file locations, as well as the desired contigs in their desired order.

submit_slopi.slurm: 
    the supercomputer partitions for your cluster and user email address, as well as what folders to bind to singularity

slurm.go_wags/slurm-submit.py: 
    partitions for your supercomputer cluster and account name