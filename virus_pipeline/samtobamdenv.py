import sys
import argparse
import glob
import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_bam(bam_file):
    logging.info(f"Validating BAM file: {bam_file}")
    run_command(f"samtools quickcheck {bam_file}")
    stdout, stderr = run_command(f"samtools view -c {bam_file}")
    alignments = int(stdout.strip())
    logging.info(f"BAM file {bam_file} validated with {alignments} alignments")
    return alignments > 0

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing SAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for BAM files.')
    args = parser.parse_args(argv)

    input_dir = os.path.abspath(args.input_dir)
    reference_fasta = os.path.abspath(args.reference_fasta)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sam_files = glob.glob(os.path.join(input_dir, "*_aln.sam"))
    logging.info(f"Found {len(sam_files)} SAM files in {input_dir}: {sam_files}")
    if not sam_files:
        logging.error(f"No SAM files found in {input_dir}")
        raise FileNotFoundError(f"No SAM files found in {input_dir}")

    for sam_file in sam_files:
        logging.info(f"Processing: {sam_file}")
        sample_name = os.path.splitext(os.path.basename(sam_file))[0].replace('_aln', '')

        # Convert SAM to BAM with MAPQ and flag filtering (Fix 5)
        # -q 20: minimum mapping quality 20
        # -F 0x904: exclude unmapped (0x4), secondary (0x100), supplementary (0x800)
        bam_file = os.path.join(output_dir, f"{sample_name}.bam")
        sam_to_bam_command = (
            f"samtools view -S -b -T {reference_fasta} "
            f"-q 20 -F 0x904 "
            f"{sam_file} > {bam_file}"
        )
        run_command(sam_to_bam_command)
        validate_bam(bam_file)

        # Add mate information for duplicate marking (Fix 2)
        fixmate_bam = os.path.join(output_dir, f"{sample_name}.fixmate.bam")
        fixmate_command = f"samtools fixmate -m {bam_file} {fixmate_bam}"
        run_command(fixmate_command)

        # Sort BAM file
        sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
        sort_command = f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {sorted_bam_file} {fixmate_bam}"
        run_command(sort_command)

        # Mark and remove duplicates (Fix 2)
        dedup_bam = os.path.join(output_dir, f"{sample_name}.dedup.bam")
        markdup_command = f"samtools markdup -r -s {sorted_bam_file} {dedup_bam}"
        run_command(markdup_command)
        os.rename(dedup_bam, sorted_bam_file)
        validate_bam(sorted_bam_file)

        # Index BAM file
        index_command = f"samtools index {sorted_bam_file}"
        run_command(index_command)

        # Clean up intermediate files
        for tmp in [bam_file, fixmate_bam]:
            if os.path.exists(tmp):
                os.remove(tmp)

        logging.info(f"Sorted, deduplicated BAM file generated and indexed: {sorted_bam_file}")

if __name__ == '__main__':
    main()
