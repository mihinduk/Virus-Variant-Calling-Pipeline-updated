#!/usr/bin/env python3
"""Extract protein sequences from consensus FASTAs with mutation tracking."""
import argparse
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def read_fasta(filepath):
    """Read a single-sequence FASTA file. Returns (header, sequence)."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
            else:
                seq_parts.append(line)
    return header, ''.join(seq_parts)


def translate(nuc_seq):
    """Translate nucleotide sequence to protein. N-containing codons -> X."""
    protein = []
    for i in range(0, len(nuc_seq) - 2, 3):
        codon = nuc_seq[i:i + 3].upper()
        if 'N' in codon:
            protein.append('X')
        else:
            protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def write_fasta_seq(f, header, sequence, line_width=60):
    """Write a single FASTA entry with wrapped sequence."""
    f.write(f'>{header}\n')
    for i in range(0, len(sequence), line_width):
        f.write(sequence[i:i + line_width] + '\n')


def parse_config_tsv(filepath):
    """Parse config review TSV for mature peptide annotations."""
    proteins = []
    in_annotations = False
    header_fields = None

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if '[TRANSCRIPT_ANNOTATIONS]' in line:
                in_annotations = True
                continue
            if not in_annotations:
                continue
            parts = line.split('\t')
            if header_fields is None:
                header_fields = parts
                continue
            if len(parts) < 5 or not parts[0].strip():
                continue

            row = dict(zip(header_fields, parts))
            feature_type = row.get('feature_type', '').strip()
            protein_name = row.get('protein_name', '').strip()

            # Skip non-mat_peptide and ancC
            if feature_type != 'mat_peptide':
                continue
            if protein_name == 'ancC':
                continue

            location = row.get('location', '').strip()
            if '..' not in location:
                continue

            start_str, end_str = location.split('..')
            start = int(start_str)
            end = int(end_str)

            proteins.append({
                'feature_id': row.get('feature_id', '').strip(),
                'protein_name': protein_name,
                'product': row.get('product', '').strip(),
                'start': start,
                'end': end,
            })

    return proteins


def main():
    parser = argparse.ArgumentParser(
        description='Extract protein sequences from consensus FASTAs')
    parser.add_argument('--consensus_dir', required=True,
                        help='Directory containing consensus .fa files')
    parser.add_argument('--config_tsv', required=True,
                        help='Config review TSV with transcript annotations')
    parser.add_argument('--reference', required=True,
                        help='Reference FASTA file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for protein FASTAs')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Parse protein definitions
    proteins = parse_config_tsv(args.config_tsv)
    if not proteins:
        logging.error("No mature peptide annotations found in config TSV")
        sys.exit(1)
    logging.info(f"Found {len(proteins)} mature peptides: "
                 f"{', '.join(p['protein_name'] for p in proteins)}")

    # Read reference
    ref_header, ref_seq = read_fasta(args.reference)
    ref_accession = ref_header.split()[0] if ref_header else 'REF'

    # Translate reference proteins
    ref_proteins = {}
    for prot in proteins:
        nuc = ref_seq[prot['start'] - 1:prot['end']]
        nt_len = len(nuc)
        if nt_len % 3 != 0:
            logging.warning(f"Reference {prot['protein_name']}: length {nt_len} "
                            f"not divisible by 3, truncating")
            nuc = nuc[:nt_len - (nt_len % 3)]
        ref_proteins[prot['protein_name']] = translate(nuc)

    # Find consensus FASTAs
    consensus_files = sorted([
        f for f in os.listdir(args.consensus_dir)
        if f.endswith('.fa') or f.endswith('.fasta')
    ])
    if not consensus_files:
        logging.error(f"No .fa/.fasta files found in {args.consensus_dir}")
        sys.exit(1)
    logging.info(f"Found {len(consensus_files)} consensus files")

    # Open output FASTA files and write reference sequences first
    fasta_handles = {}
    for prot in proteins:
        pname = prot['protein_name']
        out_path = os.path.join(args.output_dir, f"{pname}.fasta")
        fh = open(out_path, 'w')
        write_fasta_seq(fh, f"REF|{pname}|{ref_accession}", ref_proteins[pname])
        fasta_handles[pname] = fh

    # Summary data
    summary_rows = []

    # Process each sample
    for cf in consensus_files:
        sample = cf.replace('.fa', '').replace('.fasta', '')
        _, cons_seq = read_fasta(os.path.join(args.consensus_dir, cf))

        for prot in proteins:
            pname = prot['protein_name']
            ref_prot = ref_proteins[pname]

            # Extract and translate
            nuc = cons_seq[prot['start'] - 1:prot['end']]
            nt_len = len(nuc)
            if nt_len % 3 != 0:
                logging.warning(f"{sample} {pname}: length {nt_len} "
                                f"not divisible by 3, truncating")
                nuc = nuc[:nt_len - (nt_len % 3)]
            sample_prot = translate(nuc)

            # Calculate coverage (non-X positions)
            aa_len = len(sample_prot)
            x_count = sample_prot.count('X')
            callable_count = aa_len - x_count
            coverage_pct = round(callable_count / aa_len * 100, 1) if aa_len > 0 else 0

            # Find mutations (skip X positions)
            mutations = []
            for i in range(min(len(ref_prot), len(sample_prot))):
                if sample_prot[i] == 'X':
                    continue
                if sample_prot[i] != ref_prot[i]:
                    mutations.append(f"{ref_prot[i]}{i + 1}{sample_prot[i]}")

            mut_str = ','.join(mutations) if mutations else 'none'

            # Write FASTA entry
            header = f"{sample}|{pname}|mutations:{mut_str}|coverage:{coverage_pct}%"
            write_fasta_seq(fasta_handles[pname], header, sample_prot)

            # Summary row
            summary_rows.append({
                'SAMPLE': sample,
                'PROTEIN': pname,
                'AA_LENGTH': aa_len,
                'COVERAGE_PCT': coverage_pct,
                'N_MUTATIONS': len(mutations),
                'MUTATIONS': mut_str,
            })

    # Close FASTA handles
    for fh in fasta_handles.values():
        fh.close()

    # Write summary TSV
    summary_path = os.path.join(args.output_dir, 'protein_summary.tsv')
    summary_header = ['SAMPLE', 'PROTEIN', 'AA_LENGTH', 'COVERAGE_PCT',
                      'N_MUTATIONS', 'MUTATIONS']
    with open(summary_path, 'w') as f:
        f.write('\t'.join(summary_header) + '\n')
        for row in summary_rows:
            f.write('\t'.join(str(row[h]) for h in summary_header) + '\n')

    logging.info(f"Summary written to {summary_path}")
    logging.info(f"Protein FASTAs written to {args.output_dir}/")


if __name__ == '__main__':
    main()
