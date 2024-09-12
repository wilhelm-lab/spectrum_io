import argparse
import collections
import csv
import itertools
import logging
import sys
from typing import Dict, List

import numpy as np

logger = logging.getLogger(__name__)

cleavage_sites = {
    "trypsinp": (["K", "R"], []),
    "trypsin": (["K", "R"], ["P"]),
    "chymotrypsin": (["F", "W", "Y", "L"], ["P"]),
    "proteinasek": (["A", "E", "F", "I", "L", "T", "V", "W", "Y"], []),
    "elastase": (["L", "V", "A", "G"], ["P"]),
    "lys-c": (["K"], ["P"]),
    "arg-c": (["R"], ["P"]),
    "glu-c": (["E"], ["P"]),
    "v8-de": (["N", "D", "E", "Q"], ["P"]),
    "no_enzyme": ([], []),
}


def main(args):
    """Main to run digestion."""
    args = parse_args(args)
    # python digest.py --enzyme trypsinp --cleavages 0 --fasta /home/matthewt/data/MaxLFQ_benchmark/ups.fasta
    # --prosit_input ~/data/Prosit_test/ups_prosit_input.csv
    if args.prosit_input:
        writer = get_tsv_writer(args.prosit_input, delimiter=",")
        writer.writerow("modified_sequence,collision_energy,precursor_charge,fragmentation".split(","))
        prosit_input_file_with_proteins = args.prosit_input.replace(".csv", "_with_proteins.csv")
        writer_with_proteins = get_tsv_writer(prosit_input_file_with_proteins, delimiter=",")
        writer_with_proteins.writerow(
            "modified_sequence,collision_energy,precursor_charge,protein,fragmentation".split(",")
        )

        for peptide, proteins in get_peptide_to_protein_map(
            args.fasta,
            db=args.db,
            digestion=args.digestion,
            min_len=args.min_length,
            max_len=args.max_length,
            enzyme=args.enzyme,
            miscleavages=args.cleavages,
            methionine_cleavage=True,
            special_aas=list(args.special_aas),
            use_hash_key=False,
        ).items():
            if valid_prosit_peptide(peptide):
                for charge in [2, 3, 4]:
                    writer.writerow([peptide, 30, charge, args.fragmentation])
                    writer_with_proteins.writerow([peptide, 30, charge, proteins[0], args.fragmentation])

    if args.peptide_protein_map:
        with open(args.peptide_protein_map + ".params.txt", "w") as f:
            f.write(" ".join(sys.argv))
        writer = get_tsv_writer(args.peptide_protein_map, delimiter="\t")

        for peptide, proteins in get_peptide_to_protein_map(
            args.fasta,
            db="concat",
            digestion=args.digestion,
            min_len=args.min_length,
            max_len=args.max_length,
            enzyme=args.enzyme,
            miscleavages=args.cleavages,
            methionine_cleavage=True,
            special_aas=list(args.special_aas),
            use_hash_key=False,
        ).items():
            writer.writerow([peptide, ";".join(proteins)])

    if args.ibaq_map:
        writer = get_tsv_writer(args.ibaq_map, delimiter="\t")

        num_peptides_per_protein = get_num_ibaq_peptides_per_protein(args)

        for protein, num_peptides in num_peptides_per_protein.items():
            writer.writerow([protein, num_peptides])


def valid_prosit_peptide(peptide: str) -> bool:
    """Check if peptide is valid."""
    return len(peptide) <= 30 and "U" not in peptide and "X" not in peptide and "*" not in peptide


def parse_args(call):
    """Parse args."""
    apars = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument(
        "--fasta",
        metavar="F",
        required=False,
        help="""Fasta file used as input
                                                    """,
    )

    apars.add_argument(
        "--prosit_input",
        metavar="M",
        required=False,
        help="""Path to file where to write the prosit input file.
                                                    """,
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write mapping from peptides to all its proteins to
                                                         the specified file.
                                                    """,
    )

    apars.add_argument(
        "--ibaq_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write number of peptides per protein to the specified
                                                         file that meet the iBAQ criteria
                                                         (6 <= pepLen <= 30, no miscleavages).
                                                    """,
    )
    apars.add_argument(
        "--fragmentation", default=None, metavar="M", required=True, help="""Fragmentation (HCD or CID)"""
    )

    apars.add_argument("--db", metavar="M", required=False, default="concat", help="""Target/decoy/concat""")
    add_arguments(apars)
    args = apars.parse_args(call)
    return args


def add_arguments(apars):
    """Add arguments to arg parser."""
    apars.add_argument(
        "-e",
        "--enzyme",
        default="trypsin",
        metavar="E",
        help="""Type of enzyme "no_enzyme","elastase","pepsin",
                                                         "proteinasek","thermolysin","chymotrypsin",
                                                         "lys-n","lys-c","arg-c","asp-n","glu-c","trypsin",
                                                         "trypsinp".
                                                    """,
    )

    apars.add_argument(
        "-c",
        "--cleavages",
        default=2,
        metavar="C",
        type=int,
        help="""Number of allowed miss cleavages used in the search
                                                         engine (Only valid when using option -F).
                                                    """,
    )

    apars.add_argument(
        "-l",
        "--min-length",
        default=7,
        metavar="L",
        type=int,
        help="""Minimum peptide length allowed used in the search
                                                         engine (Only valid when using option -F).
                                                    """,
    )

    apars.add_argument(
        "-t",
        "--max-length",
        default=60,
        metavar="L",
        type=int,
        help="""Maximum peptide length allowed used in the search
                                                         engine (Only valid when using option -F).
                                                    """,
    )

    apars.add_argument(
        "--special-aas",
        default="KR",
        metavar="S",
        help="""Special AAs that MaxQuant uses for decoy generation.
                                                    """,
    )

    apars.add_argument(
        "--digestion",
        default="full",
        metavar="D",
        help="""Digestion mode ('full', 'semi' or 'none').
                                                    """,
    )


def write_protein_to_gene_map(fasta_file, output_file):
    """Write protein to gene map."""
    writer = csv.writer(open(output_file, "w"), delimiter="\t")
    for protein_name, _ in read_fasta_tide(fasta_file, db="target"):
        protein_id = protein_name.split("|")[1]
        gene_id = protein_name.split("|")[2].split(" ")[0]
        writer.writerow([protein_id, gene_id])


parse_until_first_space = lambda x: x.split(" ")[0]
parseprotein_name_func = lambda x: " ".join(x.split(" OS=")[0].split(" ")[1:])
parse_gene_name_func = lambda x: x.split(" GN=")[1].split(" ")[0] if "GN=" in x else ""


def parse_uniprot_id(fasta_id):
    """Parse uniprot id."""
    protein_id = parse_until_first_space(fasta_id)
    if "|" in protein_id:
        return protein_id.split("|")[1]
    else:
        return protein_id


def read_fasta_proteins(
    file_path,
    db="concat",
    parse_id=parse_until_first_space,
    parseprotein_name=parseprotein_name_func,
    parse_gene_name=parse_gene_name_func,
):
    """Read fasta proteins."""
    with open(file_path) as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if db in ["target", "concat"]:
                    yield (parse_id(line[1:]), parseprotein_name(line[1:]), parse_gene_name(line[1:]), line[1:])

                if db in ["decoy", "concat"]:
                    yield ("REV__" + parse_id(line[1:]), "", "REV__" + parse_gene_name(line[1:]), "")


def get_protein_annotations(fasta_file, parse_id):
    """Get protein annotations."""
    protein_annotations = dict()
    if not fasta_file:
        return protein_annotations

    for protein, protein_name, gen_name, fasta_header in read_fasta_proteins(fasta_file, parse_id=parse_id):
        if protein not in protein_annotations:
            protein_annotations[protein] = (protein_name, gen_name, fasta_header)

    return protein_annotations


def has_gene_names(protein_annotations) -> bool:
    """Check if gene names are provided."""
    counts = sum(1 for _, gen_name, _ in protein_annotations.values() if len(gen_name) > 0)
    return counts / len(protein_annotations) > 0.5


def read_fasta_tide(file_path, db: str = "target", parse_id=parse_until_first_space):
    """Read fasta tide."""
    read_fasta_maxquant(file_path, db, parse_id, special_aas=[], decoy_prefix="decoy_")


def read_fasta_maxquant(
    file_path, db: str = "target", parse_id=parse_until_first_space, special_aas=None, decoy_prefix: str = "REV__"
):
    """Read fasta maxquant."""
    if special_aas is None:
        special_aas = ["K", "R"]
    if db not in ["target", "decoy", "concat"]:
        sys.exit("unknown db mode: %s" % db)

    hasspecial_aas = len(special_aas) > 0
    with open(file_path) as fp:
        line = next(fp).rstrip()
        name = parse_id(line[1:])
        sequence_lines: List[str] = []
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                sequence = "".join(sequence_lines)
                if db in ["target", "concat"]:
                    yield name, sequence
                if db in ["decoy", "concat"]:
                    rev_sequence = sequence[::-1]
                    if hasspecial_aas:
                        rev_sequence = swap_special_aas(rev_sequence, special_aas)
                    yield decoy_prefix + name, rev_sequence
                name = parse_id(line[1:])
                sequence_lines = []
                continue
            sequence_lines.append(line)


read_fasta = read_fasta_maxquant


# e.g. special_aas = ['R', 'K'] transforms ABCKDEFRK into ABKCDERKF
def swap_special_aas(seq, special_aas):
    """Swaps the special AAs with its preceding amino acid, as is done in MaxQuant."""
    seq = list(seq)
    for i in range(1, len(seq)):
        if seq[i] in special_aas:
            swap_positions(seq, i, i - 1)
    seq = "".join(seq)
    return seq


def swap_positions(seq, pos1, pos2):
    """Swaps positions."""
    seq[pos1], seq[pos2] = seq[pos2], seq[pos1]


def get_protein_ids(file_path):
    """Get protein ids."""
    protein_ids = list()
    for protein_id, _ in read_fasta(file_path):
        protein_ids.append(protein_id)
    return set(protein_ids)


def get_protein_sequences(file_path, parse_id):
    """Get protein sequences."""
    protein_sequences = dict()
    for protein_id, protein_seq in read_fasta(file_path, db="concat", parse_id=parse_id):
        if protein_id not in protein_sequences:  # keep only first sequence per identifier
            protein_sequences[protein_id] = protein_seq
    return protein_sequences


def filter_fasta_file(fasta_file, filteredfasta_file, proteins):
    """Filter fasta file."""
    with open(filteredfasta_file, "w") as f:
        for prot, seq in read_fasta(fasta_file):
            if prot in proteins:
                f.write(">" + prot + "\n" + seq + "\n")


def get_peptides(
    fasta_file,
    db="concat",
    min_len=6,
    max_len=50,
    pre=None,
    not_post=None,
    digestion="full",
    miscleavages=0,
    methionine_cleavage=True,
):
    """Get peptides."""
    if pre is None:
        pre = ["K", "R"]
    if not_post is None:
        not_post = ["P"]
    for _, seq in read_fasta(fasta_file, db):
        yield from get_digested_peptides(
            seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionine_cleavage
        )


def get_digested_peptides(
    seq,
    min_len=6,
    max_len=50,
    pre=None,
    not_post=None,
    digestion="full",
    miscleavages=0,
    methionine_cleavage=True,
):
    """Get digested peptides."""
    if pre is None:
        pre = ["K", "R"]
    if not_post is None:
        not_post = ["P"]
    if digestion == "none":
        yield from non_specific_digest(seq, min_len, max_len)
    elif digestion == "semi":
        yield from semi_specific_digest(seq, min_len, max_len, pre, not_post, miscleavages, methionine_cleavage)
    else:
        yield from full_digest(seq, min_len, max_len, pre, not_post, miscleavages, methionine_cleavage)


def non_specific_digest(seq, min_len, max_len):
    """Generate non specific digestion."""
    len_s = len(seq)
    for i in range(len_s + 1):
        for j in range(i + min_len, min(len_s + 1, i + max_len + 1)):
            if j <= len_s:
                yield seq[i:j]


def semi_specific_digest(seq, min_len, max_len, pre, not_post, miscleavages, methionine_cleavage):
    """Generate semi-specific digestion."""
    len_s, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    for i in range(len_s + 1):
        is_cleavage_site = seq[min([len_s - 1, i])] in pre and seq[min([len_s - 1, i + 1])] not in not_post
        ismethionine_cleavage_site = i == 0 and methionine_cleavage
        if i == len_s or is_cleavage_site or ismethionine_cleavage_site:
            # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
            start = starts[0]
            for j in range(start, min([i + 1, len_s])):
                len_p = min([i, len_s - 1]) - j + 1
                if length_accepted(len_p):
                    yield (seq[j : i + 1])
            starts.append(i + 1)
            methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
            if len(starts) > miscleavages + 1 + methionine_cleaved or i == len_s:
                starts = starts[1 + methionine_cleaved :]
        else:  # peptides with non enzymatic C-terminal
            for start in starts:
                len_p = i - start + 1
                if length_accepted(len_p) and i + 1 not in starts:
                    yield (seq[start : i + 1])


def full_digest(seq, min_len, max_len, pre, not_post, miscleavages, methionine_cleavage):
    """Generate full digestion."""
    len_s, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    cleavage_sites = [0] if methionine_cleavage else []
    cleavage_sites.extend([i for i in range(len_s) if seq[i] in pre and seq[min([len_s - 1, i + 1])] not in not_post])
    cleavage_sites.append(len_s)
    for i in cleavage_sites:
        for start in starts:
            len_p = i - start + 1
            if length_accepted(len_p):
                yield (seq[start : i + 1])
        starts.append(i + 1)
        methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
        if len(starts) > miscleavages + 1 + methionine_cleaved:
            starts = starts[1 + methionine_cleaved :]


def get_peptide_to_protein_map(
    fasta_file,
    db="concat",
    min_len=6,
    max_len=52,
    enzyme: str = "trypsin",
    digestion="full",
    miscleavages=2,
    methionine_cleavage=True,
    use_hash_key=False,
    special_aas=None,
    parse_id=parse_until_first_space,
):
    """Get peptide to protein map."""
    pre, not_post = cleavage_sites[enzyme]

    if pre is None:
        pre = ["K", "R"]
    if not_post is None:
        not_post = ["P"]
    if special_aas is None:
        special_aas = ["K", "R"]
    peptide_to_protein_map = collections.defaultdict(list)
    protein_to_seq_map = dict()
    for protein_idx, (protein, seq) in enumerate(read_fasta(fasta_file, db, parse_id, special_aas=special_aas)):
        if (protein_idx + 1) % 10000 == 0:
            logger.info(f"Digesting protein {protein_idx + 1}")
        seen_peptides = set()
        protein_to_seq_map[protein] = seq
        # for peptide in digestfast.get_digested_peptides(seq, min_len, max_len, pre, not_post, digestion,
        # miscleavages, methionine_cleavage):
        for peptide in get_digested_peptides(
            seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionine_cleavage
        ):
            if use_hash_key:
                hash_key = peptide[:6]
            else:
                hash_key = peptide
            if hash_key not in seen_peptides:
                seen_peptides.add(hash_key)
                peptide_to_protein_map[hash_key].append(protein)

    if use_hash_key:
        return (peptide_to_protein_map, protein_to_seq_map)
    else:
        return peptide_to_protein_map


def get_peptide_to_protein_map_from_file(peptide_to_protein_map_file, use_hash_key=False):
    """Get peptide to protein map from file."""
    if use_hash_key:
        logger.warning("Hash key not supported yet, continuing without hash key...")
        use_hash_key = False
    peptide_to_protein_map = collections.defaultdict(list)
    reader = get_tsv_reader(peptide_to_protein_map_file)
    for i, row in enumerate(reader):
        if (i + 1) % 1000000 == 0:
            logger.info(f"Processing peptide  {i + 1}")

        peptide, proteins = row[0], row[1].split(";")
        if use_hash_key:
            sys.exit("Hash key not supported yet...")
            hash_key = peptide[:6]
        else:
            hash_key = peptide
        for protein in proteins:
            peptide_to_protein_map[hash_key].append(protein)
    return peptide_to_protein_map


def get_proteins(peptide_to_protein_map, peptide):
    """Get proteins."""
    peptide = peptide  # .replace("I", "L")
    if len(peptide_to_protein_map) == 2:
        hash_key = peptide[:6]
        proteins = list()
        if hash_key in peptide_to_protein_map[0]:
            for protein in peptide_to_protein_map[0][hash_key]:
                # TODO: This does not work correctly for full or partial digestion, since we might find the peptide
                # with the wrong number of enzymatic terminals
                if peptide in peptide_to_protein_map[1][protein]:
                    proteins.append(protein)
            proteins = sorted(proteins)
        return proteins
    else:
        return peptide_to_protein_map.get(peptide, [])


def get_all_proteins(peptide_to_protein_map):
    """Get all proteins."""
    seen_proteins = set()
    if len(peptide_to_protein_map) == 2:
        for _, proteins in peptide_to_protein_map[0].items():
            for protein in proteins:
                if protein not in seen_proteins:
                    seen_proteins.append(protein)
    else:
        for _, proteins in peptide_to_protein_map.items():
            for protein in proteins:
                if protein not in seen_proteins:
                    seen_proteins.append(protein)
    return list(seen_proteins)


def get_ibaq_peptide_to_protein_map(args):
    """Get ibaq peptide to protein map."""
    pre, not_post = cleavage_sites(args.enzyme)
    return get_peptide_to_protein_map(
        args.fasta,
        db="concat",
        digestion="full",
        min_len=max([6, args.min_length]),
        max_len=min([30, args.max_length]),
        enzyme=args.enzyme,
        miscleavages=0,
        methionine_cleavage=False,
        special_aas=list(args.special_aas),
    )


def get_num_ibaq_peptides_per_protein(args):
    """Get number of ibaq peptides per protein."""
    peptide_to_protein_map_ibaq = get_ibaq_peptide_to_protein_map(args)
    return get_num_peptides_per_protein(peptide_to_protein_map_ibaq)


def get_num_peptides_per_protein(peptide_to_protein_map):
    """Get number of peptides per protein."""
    num_peptides_per_protein = collections.defaultdict(int)
    for _, proteins in peptide_to_protein_map.items():
        for protein in proteins:
            num_peptides_per_protein[protein] += 1
    return num_peptides_per_protein


def is_enzymatic(aa1, aa2, pre=None, not_post=None, methionine_cleavage: bool = True) -> bool:
    """Check if is enzymatic."""
    if pre is None:
        pre = ["K", "R"]
    if not_post is None:
        not_post = ["P"]
    return aa1 == "-" or aa2 == "-" or (aa1 in pre and aa2 not in not_post) or (methionine_cleavage and aa1 == "M")


def has_miscleavage(seq, pre=None, not_post=None):
    """Check if peptide sequence has miscleavage."""
    if pre is None:
        pre = ["K", "R"]
    if not_post is None:
        not_post = ["P"]
    for i in range(len(seq) - 1):
        if is_enzymatic(seq[i], seq[i + 1], pre, not_post):
            return True
    return False


def get_tsv_reader(filename, delimiter="\t"):
    """Get tsv reader."""
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.reader(open(filename, newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.reader(open(filename, "rb"), delimiter=delimiter)


def get_tsv_writer(filename, delimiter="\t"):
    """Get tsv writer."""
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.writer(open(filename, "w", newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.writer(open(filename, "wb"), delimiter=delimiter)


if __name__ == "__main__":
    main(sys.argv[1:])
