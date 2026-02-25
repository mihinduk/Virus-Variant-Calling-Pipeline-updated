#!/usr/bin/env python

import argparse
import os
import sys
import logging
import shutil

from virus_pipeline import (
    create_samplesheet,
    map_reads,
    samtobamdenv,
    create_snpeff_database,
    sam2consensus_test2_ivar,
    variant_calling_consensus,
    summarize_result,
    summarize_snpEff,
)
from virus_pipeline.provenance import ProvenanceTracker

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_file_exists(file_path, description):
    if not os.path.exists(file_path):
        logging.error(f"{description} not found: {file_path}")
        sys.exit(1)

def check_write_permission(directory):
    try:
        test_file = os.path.join(directory, '.test_write')
        with open(test_file, 'w') as f:
            f.write('test')
        os.remove(test_file)
    except (PermissionError, OSError) as e:
        logging.error(f"No write permission for {directory}: {e}")
        sys.exit(1)

def check_tools():
    tools = ['bwa-mem2', 'samtools', 'fastp', 'fastqc', 'gatk', 'ivar', 'bcftools']
    # Check for both snpeff and snpEff
    snpeff_found = False
    for tool in ['snpeff', 'snpEff']:
        if shutil.which(tool):
            snpeff_found = True
            break
    if not snpeff_found:
        logging.error("Required tool not found: snpeff/snpEff")
        sys.exit(1)
    for tool in tools:
        if not shutil.which(tool):
            logging.error(f"Required tool not found: {tool}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Automate variant calling pipeline.')
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--reference_fasta', type=str, required=True)
    parser.add_argument('--genbank_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--database_name', type=str, default='denv1')
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file with primer coordinates for ivar trim. '
                             'If not provided, primer trimming is skipped.')
    args = parser.parse_args()

    # Check inputs and tools
    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
    check_file_exists(args.genbank_file, "GenBank file")
    check_tools()

    # Create and check output directories
    os.makedirs(args.output_dir, exist_ok=True)
    check_write_permission(args.output_dir)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)
    check_write_permission(sam_files_dir)

    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Initialize provenance tracker
    tracker = ProvenanceTracker(args.output_dir)
    tracker.set_pipeline_args(vars(args))
    tracker.detect_all_tool_versions()

    # Run pipeline steps
    try:
        logging.info("Starting create_samplesheet")
        create_samplesheet([args.input_dir, sample_sheet])
        check_file_exists(sample_sheet, "Sample sheet")
        tracker.record_step("Create sample sheet", "create_samplesheet",
                           {'input_dir': args.input_dir, 'output': sample_sheet})
        logging.info("Completed create_samplesheet")
    except Exception as e:
        tracker.record_step("Create sample sheet", "create_samplesheet",
                           {'input_dir': args.input_dir}, status="failed",
                           notes=str(e))
        logging.error(f"Failed in create_samplesheet: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting map_reads")
        map_reads(['--samplesheet', sample_sheet, '--reference', args.reference_fasta])
        tracker.record_step("Read trimming (fastp) + QC (FastQC) + Mapping (bwa-mem2)", "fastp, fastqc, bwa-mem2", {
            'fastp_qualified_quality_phred': 20,
            'fastp_length_required': 50,
            'fastp_cut_front': True,
            'fastp_cut_tail': True,
            'fastp_cut_window_size': 4,
            'fastp_cut_mean_quality': 20,
            'fastp_detect_adapter_for_pe': True,
            'fastp_correction': True,
            'fastp_overlap_len_require': 30,
            'bwa_mem2_reference': args.reference_fasta,
        })
        logging.info("Completed map_reads")
    except Exception as e:
        tracker.record_step("Read trimming + QC + Mapping", "fastp, fastqc, bwa-mem2",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in map_reads: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting samtobamdenv")
        samtobamdenv(['--input_dir', sam_files_dir, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir])
        tracker.record_step("SAM to BAM conversion + filtering + dedup", "samtools", {
            'samtools_view_min_mapq': 20,
            'samtools_view_exclude_flags': '0x904 (unmapped + secondary + supplementary)',
            'samtools_fixmate': True,
            'samtools_markdup_remove': True,
        })
        logging.info("Completed samtobamdenv")
    except Exception as e:
        tracker.record_step("SAM to BAM conversion", "samtools",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in samtobamdenv: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting create_snpeff_database")
        create_snpeff_database(['--genbank_file', args.genbank_file, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir, '--database_name', args.database_name])
        tracker.record_step("Build SnpEff database", "snpEff", {
            'genbank_file': args.genbank_file,
            'database_name': args.database_name,
        })
        logging.info("Completed create_snpeff_database")
    except Exception as e:
        tracker.record_step("Build SnpEff database", "snpEff",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in create_snpeff_database: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting sam2consensus_test2_ivar")
        sam2consensus_test2_ivar(['--input_dir', args.output_dir, '--reference_fasta', args.reference_fasta, '--output_dir', args.output_dir])
        tracker.record_step("Legacy ivar consensus (pre-existing step)", "ivar", {
            'note': 'This is the original sam2consensus_test2_ivar module',
        })
        logging.info("Completed sam2consensus_test2_ivar")
    except Exception as e:
        tracker.record_step("Legacy ivar consensus", "ivar",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in sam2consensus_test2_ivar: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # Build variant_calling_consensus args
    vcc_args = [
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', args.database_name,
    ]
    if args.primer_bed:
        vcc_args.extend(['--primer_bed', args.primer_bed])

    try:
        logging.info("Starting variant_calling_consensus")
        variant_calling_consensus(vcc_args)

        # Record primer trimming status (after step has actually run)
        if args.primer_bed:
            tracker.record_step("Primer trimming", "ivar trim", {
                'primer_bed': args.primer_bed,
                'min_quality': 20,
                'min_length': 30,
                'sliding_window': 4,
                'include_reads_no_primer': True,
            })
        else:
            tracker.record_step("Primer trimming", "ivar trim", {
                'primer_bed': 'NOT PROVIDED',
            }, status="skipped",
               notes="No --primer_bed argument provided. Primer sequences were NOT "
                     "removed from aligned reads. Variants at primer binding sites "
                     "may be masked by primer sequence. Provide a BED file with "
                     "--primer_bed to enable this step.")
        tracker.record_step("Coverage calculation", "samtools depth", {
            'min_mapping_quality': 20,
            'min_base_quality': 20,
        })
        tracker.record_step("Consensus generation", "samtools mpileup + ivar consensus", {
            'mpileup_max_depth': 10000,
            'mpileup_min_base_quality': 20,
            'mpileup_min_mapping_quality': 20,
            'ivar_min_quality': 20,
            'ivar_min_frequency_threshold': 0.5,
            'ivar_min_depth': 10,
            'ivar_ambiguous_char': 'N',
        })
        tracker.record_step("Variant calling", "gatk HaplotypeCaller", {
            'ploidy': 1,
            'standard_min_confidence': 30,
            'min_base_quality_score': 20,
        })
        tracker.record_step("VCF filtering", "gatk VariantFiltration + SelectVariants", {
            'filter_QD': '< 2.0',
            'filter_FS': '> 60.0 (strand bias)',
            'filter_MQ': '< 40.0',
            'filter_DP': '< 10',
            'select': 'PASS variants only',
        })
        tracker.record_step("Variant annotation", "snpEff + SnpSift", {
            'database': args.database_name,
            'input': 'filtered PASS variants',
        })
        logging.info("Completed variant_calling_consensus")
    except Exception as e:
        tracker.record_step("Variant calling + consensus", "multiple",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in variant_calling_consensus: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting summarize_result")
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir])
        tracker.record_step("Coverage summary", "summarize_result", {
            'metrics': 'mean, median, min, max, std, %genome at 1x/10x/30x/100x',
        })
        logging.info("Completed summarize_result")
    except Exception as e:
        tracker.record_step("Coverage summary", "summarize_result",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in summarize_result: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting summarize_snpEff")
        summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir])
        tracker.record_step("SnpEff summary", "summarize_snpEff", {})
        logging.info("Completed summarize_snpEff")
    except Exception as e:
        tracker.record_step("SnpEff summary", "summarize_snpEff",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in summarize_snpEff: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # Write provenance reports
    tracker.write_json()
    tracker.write_report()

    logging.info("Pipeline completed successfully!")
    logging.info(f"Provenance report: {os.path.join(args.output_dir, 'provenance_report.txt')}")
    logging.info(f"Provenance JSON:   {os.path.join(args.output_dir, 'provenance.json')}")

if __name__ == "__main__":
    main()
