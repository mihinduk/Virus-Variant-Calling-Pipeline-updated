import sys
import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_fasta(fasta_file):
    valid_bases = set("ACGTNacgtnRYKMSWBDHVrykmswbvh")
    with open(fasta_file, "r") as f:
        sequence = ""
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(f"Non-IUPAC bases found in {fasta_file}: {invalid_bases}")
    logging.info(f"FASTA file {fasta_file} validated: contains only IUPAC bases.")

def validate_bam(bam_file):
    """Validate BAM file using samtools quickcheck and check alignment count."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file {bam_file} does not exist")
    validate_command = f"samtools quickcheck {bam_file}"
    stdout, stderr = run_command(validate_command)
    if stderr:
        raise Exception(f"BAM validation failed for {bam_file}: {stderr}")
    # Check alignment count
    count_command = f"samtools view -c {bam_file}"
    stdout, stderr = run_command(count_command)
    alignment_count = int(stdout.strip())
    if alignment_count == 0:
        raise ValueError(f"BAM file {bam_file} contains no alignments")
    logging.info(f"BAM file {bam_file} validated successfully with {alignment_count} alignments")

def plot_coverage(coverage_file, output_dir):
    positions = []
    depths = []
    with open(coverage_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            position = int(columns[1])
            depth = int(columns[2])
            positions.append(position)
            depths.append(depth)

    plt.plot(positions, depths)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title('Read Coverage')
    plt.grid(True)

    coverage_plot_file = os.path.join(output_dir, os.path.splitext(os.path.basename(coverage_file))[0] + '.png')
    plt.savefig(coverage_plot_file)
    plt.close()

    logging.info(f"Coverage plot generated: {coverage_plot_file}")

def prepare_reference(reference_fasta, output_dir):
    validate_fasta(reference_fasta)
    if not os.path.exists(f"{reference_fasta}.fai"):
        logging.info("Indexing reference FASTA...")
        run_command(f"samtools faidx {reference_fasta}")

    dict_file = os.path.splitext(reference_fasta)[0] + ".dict"
    if not os.path.exists(dict_file):
        logging.info("Creating sequence dictionary for reference...")
        run_command(f"gatk CreateSequenceDictionary -R {reference_fasta} -O {dict_file}")

def add_read_groups(bam_file, sample_name, output_dir):
    rg_bam = os.path.join(output_dir, f"{sample_name}_rg.bam")
    rg_command = (
        f"samtools addreplacerg -r 'ID:{sample_name}\tSM:{sample_name}\tLB:lib1\tPL:ILLUMINA' "
        f"-o {rg_bam} {bam_file}"
    )
    stdout, stderr = run_command(rg_command)
    print("Add read groups output:", stdout)
    print("Add read groups error:", stderr)

    # Validate and sort the BAM file
    run_command(f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {rg_bam}.sorted {rg_bam}")
    os.rename(f"{rg_bam}.sorted", rg_bam)
    run_command(f"samtools index {rg_bam}")
    return rg_bam

def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir):
    validate_bam(bam_file)
    raw_vcf = os.path.join(output_dir, f"{sample_name}.vcf")
    # Fix 3: Added -ploidy 1 for haploid virus genome
    gatk_command = (
        f"gatk HaplotypeCaller "
        f"-R {reference_fasta} "
        f"-I {bam_file} "
        f"-O {raw_vcf} "
        f"-ploidy 1 "
        f"--standard-min-confidence-threshold-for-calling 30 "
        f"--min-base-quality-score 20"
    )
    stdout, stderr = run_command(gatk_command)
    logging.info(f"GATK HaplotypeCaller output: {stdout}")
    logging.info(f"GATK HaplotypeCaller error: {stderr}")
    return raw_vcf

def filter_vcf(raw_vcf, reference_fasta, sample_name, output_dir):
    """Fix 4: Filter VCF for quality, depth, and strand bias."""
    filtered_vcf = os.path.join(output_dir, f"{sample_name}_filtered.vcf")
    filter_command = (
        f"gatk VariantFiltration "
        f"-R {reference_fasta} "
        f"-V {raw_vcf} "
        f"--filter-expression 'QD < 2.0' --filter-name 'LowQD' "
        f"--filter-expression 'FS > 60.0' --filter-name 'StrandBias' "
        f"--filter-expression 'MQ < 40.0' --filter-name 'LowMQ' "
        f"--filter-expression 'DP < 10' --filter-name 'LowDepth' "
        f"-O {filtered_vcf}"
    )
    run_command(filter_command)

    # Select only PASS variants
    pass_vcf = os.path.join(output_dir, f"{sample_name}_pass.vcf")
    select_command = (
        f"gatk SelectVariants "
        f"-R {reference_fasta} "
        f"-V {filtered_vcf} "
        f"--exclude-filtered "
        f"-O {pass_vcf}"
    )
    run_command(select_command)
    logging.info(f"VCF filtering complete: {pass_vcf}")
    return pass_vcf

def trim_primers(bam_file, primer_bed, sample_name, output_dir):
    """Fix 7: Trim primer sequences from aligned reads using ivar trim."""
    trimmed_prefix = os.path.join(output_dir, f"{sample_name}_trimmed")
    trimmed_bam = f"{trimmed_prefix}.bam"
    trimmed_sorted = os.path.join(output_dir, f"{sample_name}_trimmed.sorted.bam")

    trim_command = (
        f"ivar trim -b {primer_bed} "
        f"-p {trimmed_prefix} "
        f"-i {bam_file} "
        f"-q 20 -m 30 -s 4 -e"
    )
    run_command(trim_command)

    # Sort and index the trimmed BAM
    run_command(f"samtools sort -o {trimmed_sorted} {trimmed_bam}")
    run_command(f"samtools index {trimmed_sorted}")

    # Clean up unsorted trimmed BAM
    if os.path.exists(trimmed_bam):
        os.remove(trimmed_bam)

    logging.info(f"Primer trimming complete: {trimmed_sorted}")
    return trimmed_sorted

def run_snpeff_annotation(raw_vcf, sample_name, output_dir, reference_name="denv1"):
    annotated_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    summary_html = os.path.join(output_dir, f"{sample_name}_snpEff_summary.html")
    summary_csv = os.path.join(output_dir, f"{sample_name}_snpEff_summary.csv")
    snpeff_command = (
        f"snpEff -Xmx4g "
        f"-c {os.path.join(output_dir, 'snpEff.config')} "
        f"-v {reference_name} "
        f"-s {summary_html} "
        f"-csvStats {summary_csv} "
        f"{raw_vcf} > {annotated_vcf}"
    )
    stdout, stderr = run_command(snpeff_command)
    logging.info(f"SnpEff annotation output: {stdout}")
    logging.info(f"SnpEff annotation error: {stderr}")
    return annotated_vcf, summary_html, summary_csv, os.path.join(output_dir, f"{sample_name}_snpEff_genes.txt")

def run_snpsift_extract(annotated_vcf, sample_name, output_dir):
    snpsift_output = os.path.join(output_dir, f"{sample_name}_snpSift.txt")
    snpsift_command = (
        f"SnpSift extractFields {annotated_vcf} "
        f"CHROM POS REF ALT "
        f"\"ANN[*].EFFECT\" \"ANN[*].IMPACT\" \"ANN[*].GENE\" \"ANN[*].GENEID\" "
        f"\"ANN[*].FEATURE\" \"ANN[*].HGVS_C\" \"ANN[*].HGVS_P\" \"ANN[*].AA_POS\" "
        f"\"EFF[*].CODON\" \"EFF[*].AA\" \"EFF[*].GENE\" > {snpsift_output}"
    )
    stdout, stderr = run_command(snpsift_command)
    logging.info(f"SnpSift extract output: {stdout}")
    logging.info(f"SnpSift extract error: {stderr}")
    return snpsift_output

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA and VCF files.')
    parser.add_argument('--database_name', type=str, default="denv1", help='Name of the SnpEff database (default: denv1).')
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file with primer coordinates for ivar trim (Fix 7). '
                             'If not provided, primer trimming is skipped.')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    database_name = args.database_name
    primer_bed = args.primer_bed

    if primer_bed and not os.path.exists(primer_bed):
        logging.error(f"Primer BED file not found: {primer_bed}")
        sys.exit(1)
    if primer_bed:
        logging.info(f"Primer trimming enabled with BED file: {primer_bed}")
    else:
        logging.warning("No primer BED file provided -- skipping primer trimming (Fix 7)")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files = glob.glob(os.path.join(input_dir, "*.sorted.bam"))
    if not bam_files:
        logging.error(f"No sorted BAM files found in {input_dir}")
        sys.exit(1)

    prepare_reference(reference_fasta, output_dir)
    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))
    for bam_file in bam_files:
        print("Processing:", bam_file)
        sample_name = os.path.splitext(os.path.basename(bam_file))[0].replace('.sorted', '')
        try:
            # Add read groups before GATK
            rg_bam = add_read_groups(bam_file, sample_name, output_dir)

            # Fix 7: Primer trimming (if BED file provided)
            analysis_bam = rg_bam
            if primer_bed:
                analysis_bam = trim_primers(rg_bam, primer_bed, sample_name, output_dir)

            # Coverage calculation with quality filters (Fix 5 continued)
            coverage_file = os.path.join(output_dir, f"{sample_name}_coverage.txt")
            coverage_command = (
                f"samtools depth -q 20 -Q 20 {analysis_bam} > {coverage_file}"
            )
            stdout, stderr = run_command(coverage_command)
            logging.info(f"Coverage command output: {stdout}")
            logging.info(f"Coverage command error: {stderr}")
            plot_coverage(coverage_file, output_dir)

            # Fix 6: ivar consensus with majority-rule threshold and minimum depth
            pileup_command = (
                f"samtools mpileup -d 10000 -Q 20 -q 20 "
                f"{analysis_bam} | "
                f"ivar consensus -p {sample_name} -q 20 -t 0.5 -m 10 -n N"
            )
            stdout, stderr = run_command(pileup_command)
            logging.info(f"Consensus command output: {stdout}")
            logging.info(f"Consensus command error: {stderr}")
            consensus_fasta = f"{sample_name}.fa"
            os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))

            # Variant calling (Fix 3: ploidy=1 applied inside function)
            raw_vcf = run_variant_calling(analysis_bam, reference_fasta, sample_name, output_dir)

            # Fix 4: VCF filtering before annotation
            filtered_vcf = filter_vcf(raw_vcf, reference_fasta, sample_name, output_dir)

            # Annotate filtered variants
            annotated_vcf, summary_html, summary_csv, summary_txt = run_snpeff_annotation(
                filtered_vcf, sample_name, output_dir, database_name)
            snpsift_output = run_snpsift_extract(annotated_vcf, sample_name, output_dir)
            logging.info(
                f"Processing complete for {sample_name}: "
                f"consensus={os.path.join(output_dir, consensus_fasta)}, "
                f"coverage={coverage_file}, filtered_vcf={filtered_vcf}, "
                f"annotated_vcf={annotated_vcf}"
            )
        except Exception as e:
            print(f"Error occurred during processing {sample_name}: {str(e)}")
            continue

if __name__ == '__main__':
    main()
