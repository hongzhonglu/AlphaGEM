"""Match complete COBRA gene tokens when transferring reference GPRs."""

import re


def replace_gene_token(rule, gene, replacement):
    # COBRA separates gene IDs with whitespace and parentheses. A word boundary
    # would also match YER060W inside the distinct locus YER060W-A.
    pattern = r"(?<![^\s()])" + re.escape(gene) + r"(?![^\s()])"
    return re.sub(pattern, lambda match: replacement, rule)


def has_gene_token(rule, gene):
    return replace_gene_token(rule, gene, "") != rule
