"""Independent Rhea-based annotation path from commit 5048e62.

This module needs only the original Rhea RDF, KEGG-to-Rhea table and
general2specific pickle. It never imports the experimental annotation database
or its resource loader. Scientific rules are unchanged: source weights,
CLEAN > 0.7, direct RDF EC associations (LIMIT 100), and PLM support for generic
reaction children. Repairs are limited to identifier normalization, empty
outputs, duplicate source rows and copying a parent row for each supported leaf.
"""

from pathlib import Path
import math
import pickle
import re

import pandas as pd
from rdflib import Graph, URIRef


ROOT = Path(__file__).resolve().parents[1]
SCORE_COLUMNS = ["reaction", "rhea", "anno_A", "anno_B", "anno_C", "anno_D", "final_score"]
SOURCE_COLUMNS = {"clean": "anno_A", "eggnog_mapper": "anno_B",
                  "deepec": "anno_C", "plmsearch": "anno_D"}
SOURCE_WEIGHTS = {"anno_A": 0.2, "anno_B": 0.3, "anno_C": 0.2, "anno_D": 0.1}


def _gene_id(value):
    value = str(value).strip()
    parts = value.split("|")
    return parts[1] if len(parts) > 2 else value


def _rhea_id(value):
    if value is None or pd.isna(value):
        return ""
    value = str(value).strip().removeprefix("RHEA:")
    if re.fullmatch(r"\d+\.0+", value):
        value = value.split(".", 1)[0]
    return value if value.isdigit() else ""


def _directories(name, work_dir, output_dir, resource_root):
    work_dir = Path(work_dir) if work_dir is not None else ROOT / "working" / name
    output_dir = Path(output_dir) if output_dir is not None else work_dir
    resource_root = Path(resource_root) if resource_root is not None else ROOT
    output_dir.mkdir(parents=True, exist_ok=True)
    return work_dir, output_dir, resource_root


def calculate_score(row):
    """Original once-per-tool weights; no substrate model or extra votes."""
    return round(sum(weight for column, weight in SOURCE_WEIGHTS.items()
                     if pd.notnull(row.get(column))), 1)


def findtargetreaction(g, ec):
    """Use only the original RDF association, including its 100-result limit."""
    query = """
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?x ?reaction WHERE {
            ?x rdfs:subClassOf rh:Reaction .
            ?x rh:ec ?enzyme .
        }
        LIMIT 100
    """
    enzyme = URIRef("http://purl.uniprot.org/enzyme/" + str(ec).strip())
    return [str(row.x).rsplit("/", 1)[-1]
            for row in g.query(query, initBindings={"enzyme": enzyme})]


def nonhome(name, clean_use=True, deepec_use=True, plm_use=True, ecgraph=False,
            *, work_dir=None, output_dir=None, resource_root=None):
    """Write the legacy source-vote table; return it for inspection.

    Inputs are read from work_dir; output_dir only controls output. When using
    isolated output, pass that directory as get_spe_reactions(work_dir=...).
    """
    work_dir, output_dir, resource_root = _directories(name, work_dir, output_dir, resource_root)
    graph = Graph()
    graph.parse(resource_root / "rhea/rhea.rdf", format="xml")
    kegg_table = pd.read_csv(resource_root / "data_available/rhea2kegg_reaction.tsv", sep="\t")
    kegg_to_rhea = kegg_table.set_index("ID")["MASTER_ID"].to_dict()
    records, ec_cache = [], {}

    def add(gene, rhea, source):
        identifier = _rhea_id(rhea)
        if identifier and gene and gene not in {"nan", "None"}:
            records.append({"reaction": gene, "rhea": identifier, "source": source})

    def ec_reactions(ec):
        if ec not in ec_cache:
            ec_cache[ec] = findtargetreaction(graph, ec)
        return ec_cache[ec]

    egg = pd.read_excel(work_dir / f"eggec2{name}.xlsx", dtype=str).fillna("")
    columns = [column for column in egg.columns if not str(column).startswith("Unnamed:")]
    if "Rec" not in egg or not columns:
        raise ValueError("Legacy eggNOG input requires a gene column and Rec column")
    gene_columns = [column for column in columns if str(column) not in {"EC", "Rec"}]
    if not gene_columns:
        raise ValueError("Legacy eggNOG input is missing its gene column")
    for _, row in egg.iterrows():
        gene = _gene_id(row[gene_columns[0]])
        # Preserve the original KEGG-reaction route. An EC-only eggNOG row
        # cannot introduce a new EC vote in the legacy workflow.
        for kegg in str(row["Rec"]).split(","):
            if kegg.strip() in kegg_to_rhea:
                add(gene, kegg_to_rhea[kegg.strip()], "eggnog_mapper")

    if clean_use:
        with (work_dir / f"{name}_homoleft_maxsep.csv").open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                gene = _gene_id(re.split(r"[,\t ]", line.strip(), maxsplit=1)[0])
                for ec, raw_score in re.findall(r"EC:([^,\s/]+)/([0-9.eE+-]+)", line):
                    score = float(raw_score)
                    if math.isfinite(score) and score > 0.7:
                        for rhea in ec_reactions(ec):
                            add(gene, rhea, "clean")

    if deepec_use:
        deep = pd.read_csv(work_dir / f"{name}_deepec_result/DeepECv2_result.txt", sep="\t", dtype=str).fillna("")
        if "sequence_ID" not in deep or len(deep.columns) < 2:
            raise ValueError("Legacy DeepEC input requires sequence_ID and a prediction column")
        for _, row in deep.iterrows():
            for token in str(row.iloc[1]).split(";"):
                token = token.strip()
                if not token.startswith("EC:") or token.endswith("-"):
                    continue
                for rhea in ec_reactions(token[3:]):
                    add(_gene_id(row["sequence_ID"]), rhea, "deepec")

    if plm_use:
        plm = pd.read_excel(work_dir / f"{name}_rheaid_plmsearch.xlsx", dtype=str).fillna("")
        if not {"Query", "Rhea_id"} <= set(plm.columns):
            raise ValueError("Legacy PLM input requires Query and Rhea_id columns")
        for _, row in plm.iterrows():
            add(_gene_id(row["Query"]), row["Rhea_id"], "plmsearch")
    if ecgraph:
        # 5048e62 read this file but never added ECgraph to the evidence table.
        # Keep that behavior explicit instead of inventing a fifth vote.
        pd.read_excel(work_dir / f"{name}_ecgraphsearch.xlsx")

    votes = {}
    for record in records:
        key = record["reaction"], record["rhea"]
        row = votes.setdefault(key, {"reaction": key[0], "rhea": key[1],
                                     **{column: None for column in SOURCE_WEIGHTS}})
        row[SOURCE_COLUMNS[record["source"]]] = record["source"]
    rows = [{**row, "final_score": calculate_score(row)} for _, row in sorted(votes.items())]
    result = pd.DataFrame(rows, columns=SCORE_COLUMNS)
    result.to_excel(output_dir / f"{name}multiple_score.xlsx", index=False)
    return result


def filter_plmsearch(merged, general2spe, filter_num, done=True):
    """Original PLM-supported generic-to-leaf policy with stable ID types."""
    if not 0 <= filter_num <= sum(SOURCE_WEIGHTS.values()):
        raise ValueError("Legacy ensemble cutoff must be between 0 and 0.8")
    merged = merged.copy()
    for column in SCORE_COLUMNS:
        if column not in merged:
            merged[column] = pd.Series(index=merged.index, dtype="object")
    merged["rhea"] = merged["rhea"].map(_rhea_id)
    merged = merged.loc[merged["rhea"] != ""]
    mapping = {_rhea_id(parent): [_rhea_id(child) for child in children]
               for parent, children in general2spe.items()}
    supported = {(str(row["reaction"]), row["rhea"])
                 for _, row in merged.iterrows() if row["anno_D"] == "plmsearch"}
    results = []
    for _, original in merged.iterrows():
        children = mapping.get(original["rhea"])
        if children is None:
            results.append(original.copy())
            continue
        for child in children:
            if not child or (done and (str(original["reaction"]), child) not in supported):
                continue
            row = original.copy()
            row["rhea"] = child
            results.append(row)
    if not results:
        return pd.DataFrame(columns=SCORE_COLUMNS)
    result = pd.DataFrame(results).drop_duplicates()
    # groupby.first combines the source votes exactly as the original did.
    result = result.groupby(["reaction", "rhea"], as_index=False).first()
    result["final_score"] = result.apply(calculate_score, axis=1)
    return result.loc[result["final_score"] >= filter_num, SCORE_COLUMNS].reset_index(drop=True)


def get_spe_reactions(name, filter_num=0.4, *, work_dir=None, output_dir=None, resource_root=None):
    """Read legacy scores and write Rhea/GPR output using the original pickle."""
    work_dir, output_dir, resource_root = _directories(name, work_dir, output_dir, resource_root)
    merged = pd.read_excel(work_dir / f"{name}multiple_score.xlsx").drop_duplicates()
    # This is the repository's original trusted local asset, not user input.
    with (resource_root / "data_available/general2specific.pkl").open("rb") as handle:
        general2spe = pickle.load(handle)
    resolved = filter_plmsearch(merged, general2spe, filter_num, done=True)
    if resolved.empty:
        result = pd.DataFrame(columns=["rhea", "reaction"])
    else:
        result = resolved.groupby("rhea")["reaction"].apply(
            lambda genes: " or ".join(dict.fromkeys(str(gene) for gene in genes))
        ).reset_index()
    # Keep the original Excel index column for historical downstream readers.
    result.to_excel(output_dir / f"{name}gpr_score.xlsx")
    return result
