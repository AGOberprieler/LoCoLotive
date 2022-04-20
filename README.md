# LoCoLotive
In silico identification of low-copy nuclear loci based on published target capture probe sets and custom reference genomes.

## Requirements
- Linux OS
- [Docker](https://docs.docker.com/) or ([BEDOPS](https://bedops.readthedocs.io/en/latest/), [bedtools](https://bedtools.readthedocs.io/en/latest/), [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), [GenomeTools](http://genometools.org/), [MAFFT](https://mafft.cbrc.jp/alignment/software/), [BLAST](https://blast.ncbi.nlm.nih.gov/), [python3](https://www.python.org/), [R](https://www.r-project.org/)).

Depending on your Linux distribution, the latter tools can possibly be installed with

```raw
apt update && apt install bedops bedtools fastx-toolkit genometools mafft ncbi-blast+ python3 r-base
```

Alternatively, you can use the supplied Docker image, which already contains all dependencies required:

1. Install [Docker](https://docs.docker.com/engine/install/).
2. Run `docker pull ulbio/probdist` to download and import the docker image.

Since compatibility with future versions of the external tools cannot be guaranteed, it is safer to use the Docker image rather than local installations.
Note that this image does NOT support the Windows Subsystem for Linux (WSL).
