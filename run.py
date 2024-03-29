#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import hashlib


def md5sum(fname: str) -> str:
    """Calculate MD5 checksum.

    Parameters
    ----------
    fname : str
        file name

    Returns
    -------
    str
        checksum
    """
    md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5.update(chunk)
    return md5.hexdigest()


def confirmation(question: str) -> bool:
    """Ask question expecting "yes" or "no".

    Parameters
    ----------
    question : str
        Question to be printed

    Returns
    -------
    bool
        True or False for "yes" or "no", respectively
    """
    valid = {"yes": True, "y": True, "no": False, "n": False}

    while True:
        choice = input(question + " [y/n]\n").strip().lower()
        if choice in valid.keys():
            return valid[choice]
        print("Please type 'yes' or 'no'\n")


def main():

    parser = argparse.ArgumentParser()

    # obligate parameters:
    parser.add_argument(
        "targets",
        help = "Target sequences used for probe design (FASTA file)"
    )

    parser.add_argument(
        "genome",
        help = "Reference genome (FASTA file)"
    )

    # optional parameters:
    parser.add_argument(
        "-a", "--annotation",
        required = False,
        default = "",
        help = "Genome annotation (GFF3 file, optional)"
    )

    parser.add_argument(
        "-e", "--evalue",
        required = False,
        default = "1e-5",
        help = "E-value threshold for saving BLAST hits (default: 1e-5)"
    )

    parser.add_argument(
        "-i", "--relax_intron_def",
        required = False,
        default = False,
        action = "store_true",
        help = "Relax definition of intronic regions to include even exon-overlapping parts. (legacy option)"
    )

    parser.add_argument(
        "-m", "--mc_length",
        required = False,
        type = int,
        default = 15,
        help = "Targets with multi-copy regions of at least MC_LENGTH bp will be discarded. (default: 15)"
    )

    parser.add_argument(
        "-r", "--run_all",
        required = False,
        default = False,
        action = "store_true",
        help = "(Re)run all intermediate steps. This may be useful, if previous runs have been interrupted."
    )

    cfg = parser.parse_args()

    fasta_extensions = (
        ".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa", ".fas",
        ".FASTA", ".FNA", ".FFN", ".FAA", ".FRN", ".FA", ".FAS"
    )

    for seq in ["genome", "targets"]:
        _, ext = os.path.splitext(getattr(cfg, seq))

        if ext not in fasta_extensions:
            raise ValueError(f'\"{getattr(cfg, seq)}\" must be FASTA file and have an appropriate file name extension (eg. ".fasta")')

    genome_fname = os.path.split(cfg.genome)[1]
    genome_name = os.path.splitext(os.path.split(cfg.genome)[1])[0]
    targets_name = os.path.splitext(os.path.split(cfg.targets)[1])[0]
    os.makedirs(
        os.path.join(
            targets_name,
            genome_name,
            "e_thresh_" + cfg.evalue,
            "mc_thresh_" + str(cfg.mc_length)
        ),
        exist_ok=True
    )


    # consistency checks:

    print("check input files...")

    md5file = os.path.join(targets_name, f"{targets_name}.md5")
    md5 = md5sum(cfg.targets)

    if not os.path.isfile(md5file):
        with open(md5file, "w") as f:
            f.write(md5)
    else:
        with open(md5file, "r") as f:
            md5_old = f.read()

        if md5 != md5_old:
            print(f"\nWARNING: {cfg.targets} has changed!\nThis may lead to inconsistent results/subdirectories in {targets_name}/.")
            if not confirmation("Do you want to proceed? (existing checksums will be overwritten)"):
                sys.exit()
            else:
                with open(md5file, "w") as f:
                    f.write(md5)


    md5file = os.path.join(targets_name, genome_name, f"{genome_name}.md5")
    md5 = md5sum(cfg.genome)

    if not os.path.isfile(md5file):
        with open(md5file, "w") as f:
            f.write(md5)
    else:
        with open(md5file, "r") as f:
            md5_old = f.read()

        if md5 != md5_old:
            print(f"\nWARNING: {cfg.genome} has changed!\nThis may lead to inconsistent results/subdirectories in {os.path.join(targets_name, genome_name)}/.")
            if not confirmation("Do you want to proceed? (existing checksums will be overwritten)"):
                sys.exit()
            else:
                with open(md5file, "w") as f:
                    f.write(md5)


    if cfg.annotation:
        md5file = os.path.join(targets_name, genome_name, f"{genome_name}_ann.md5")
        md5 = md5sum(cfg.annotation)

        if not os.path.isfile(md5file):
            with open(md5file, "w") as f:
                f.write(md5)
        else:
            with open(md5file, "r") as f:
                md5_old = f.read()

            if md5 != md5_old:
                print(f"\nWARNING: {cfg.annotation} has changed!\nThis may lead to inconsistent results/subdirectories in {os.path.join(targets_name, genome_name)}/.")
                if not confirmation("Do you want to proceed? (existing checksums will be overwritten)"):
                    sys.exit()
                else:
                    with open(md5file, "w") as f:
                        f.write(md5)


    # run pipeline:

    if cfg.annotation:
        if (
            cfg.run_all
            or not os.path.isfile(os.path.join(targets_name, genome_name, "intronic.BED"))
            or not os.path.isfile(os.path.join(targets_name, genome_name, "intronic_strict.BED"))
        ):
            print("extract intronic regions...")

            subprocess.run(
                [
                    "./extract_introns.sh",
                    cfg.annotation,
                    os.path.join(targets_name, genome_name)
                ],
                stderr=subprocess.STDOUT,
                stdout=sys.stdout,
                check=True
            )
        else:
            print("existing intronic.BED and intronic_strict.BED found")


    db_files = [
        os.path.join(targets_name, genome_name, genome_fname + ext)
        for ext in [".nhr", ".nin", ".nsq"]
    ]

    if cfg.run_all or not all(os.path.isfile(f) for f in db_files):
        print("\ncreate BLAST database...")

        subprocess.run(
            [
                "makeblastdb",
                "-in", cfg.genome,
                "-input_type", "fasta",
                "-out", os.path.join(
                    targets_name,
                    genome_name,
                    genome_name + os.path.splitext(cfg.genome)[1]
                ),
                "-max_file_sz", "2GB",
                "-dbtype", "nucl"
            ],
            stderr=subprocess.STDOUT,
            stdout=sys.stdout,
            check=True
        )

    else:
        print("existing BLAST database found")

    if cfg.run_all or not os.path.isfile(
        os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "blast_hits.txt")
    ):
        print("\nBLAST target sequences against reference genome (this may take some time)...")

        with open(os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "blast_hits.txt"), "wb") as f:
            subprocess.run(
                [
                    "blastn",
                    "-task", "dc-megablast",
                    "-query", cfg.targets,
                    "-db", os.path.join(targets_name, genome_name, genome_fname),
                    "-evalue", cfg.evalue,
                    "-outfmt", "7 qseqid qstart qend sseqid sstart send sstrand bitscore evalue pident length qcovhsp mismatch gapopen gaps btop"
                ],
                stderr=subprocess.STDOUT,
                stdout=f,
                check=True
            )
    else:
        print("existing BLAST results found")

    print("\nfilter target sequences:\n")
    subprocess.run(
        [
            "./filter.sh",
            os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "blast_hits.txt"),
            os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "mc_thresh_" + str(cfg.mc_length)),
            str(cfg.mc_length)
        ],
        stderr=subprocess.STDOUT,
        stdout=sys.stdout,
        check=True
    )

    if 0 == len(os.listdir(os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "mc_thresh_" + str(cfg.mc_length), "hits_filtered"))):
        print("\n\nWARNING: No loci left, abort.")
        sys.exit()

    print("\ncreate alignments...")
    if cfg.annotation:
        subprocess.run(
            [
                "./alignments_multi.sh",
                cfg.targets,
                cfg.genome,
                os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "mc_thresh_" + str(cfg.mc_length)),
                os.path.join(
                    targets_name,
                    genome_name,
                    "intronic.BED" if cfg.relax_intron_def else "intronic_strict.BED"
                )
            ],
            stderr=subprocess.STDOUT,
            stdout=sys.stdout,
            check=True
        )
    else:
        subprocess.run(
            [
                "./alignments_multi.sh",
                cfg.targets,
                cfg.genome,
                os.path.join(targets_name, genome_name, "e_thresh_" + cfg.evalue, "mc_thresh_" + str(cfg.mc_length))
            ],
            stderr=subprocess.STDOUT,
            stdout=sys.stdout,
            check=True
        )


if __name__ == "__main__":
    main()
