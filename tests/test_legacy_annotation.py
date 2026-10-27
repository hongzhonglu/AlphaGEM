"""Original Rhea annotation must run without experimental database assets."""

import ast
from pathlib import Path
import pickle

import pandas as pd
import pytest
from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDFS

from src import legacy_annotation as legacy


@pytest.fixture
def legacy_inputs(tmp_path):
    resources, work, output = tmp_path / "resources", tmp_path / "working", tmp_path / "isolated"
    for directory in (resources / "rhea", resources / "data_available", work, output):
        directory.mkdir(parents=True, exist_ok=True)
    rh = Namespace("http://rdf.rhea-db.org/")
    graph = Graph()
    for identifier, ec in (("10000", "1.1.1.1"), ("10001", "1.1.1.1"), ("20000", "2.2.2.2"), ("30000", "3.3.3.3")):
        node = rh[identifier]
        graph.add((node, RDFS.subClassOf, rh.Reaction))
        graph.add((node, rh.ec, URIRef("http://purl.uniprot.org/enzyme/" + ec)))
    # Merely carrying EC without being a direct Reaction subclass is insufficient.
    graph.add((rh["99999"], rh.ec, URIRef("http://purl.uniprot.org/enzyme/1.1.1.1")))
    graph.serialize(destination=str(resources / "rhea/rhea.rdf"), format="xml")
    pd.DataFrame({"ID": ["R00001"], "MASTER_ID": [10000]}).to_csv(
        resources / "data_available/rhea2kegg_reaction.tsv", sep="\t", index=False)
    with (resources / "data_available/general2specific.pkl").open("wb") as handle:
        pickle.dump({"10000": ["20000", "30000"]}, handle)
    name = "fixture"
    pd.DataFrame({0: ["gene1", "ec_only"], "EC": ["1.1.1.1", "3.3.3.3"],
                  "Rec": ["R00001", "-"]}).to_excel(work / f"eggec2{name}.xlsx")
    (work / f"{name}_homoleft_maxsep.csv").write_text(
        "sp|gene1|gene1,EC:1.1.1.1/0.8\n"
        "sp|gene1|gene1,EC:1.1.1.1/0.8\n"
        "sp|cutoff|cutoff,EC:3.3.3.3/0.7\n"
    )
    deep_dir = work / f"{name}_deepec_result"
    deep_dir.mkdir()
    pd.DataFrame({"sequence_ID": ["sp|gene1|gene1", "sp|partial|partial"],
                  "prediction": ["EC:1.1.1.1", "EC:3.3.3.-"]}).to_csv(
        deep_dir / "DeepECv2_result.txt", sep="\t", index=False)
    pd.DataFrame({"Query": ["sp|gene1|gene1", "sp|gene1|gene1", "sp|other|other"],
                  "Rhea_id": [20000.0, 30000.0, 20000.0]}).to_excel(work / f"{name}_rheaid_plmsearch.xlsx")
    return name, resources, work, output, graph


def test_rdf_query_matches_original_without_namespace_accident(legacy_inputs):
    *_, graph = legacy_inputs
    # The historical query relied on namespace bindings from the RDF file.
    # Supply those explicitly here and compare the exact query semantics.
    original_query = """
        SELECT DISTINCT ?x ?reaction WHERE {
            ?x rdfs:subClassOf rh:Reaction .
            ?x rh:ec <http://purl.uniprot.org/enzyme/1.1.1.1> .
        } LIMIT 100
    """
    expected = {str(row.x).rsplit("/", 1)[-1] for row in graph.query(
        original_query, initNs={"rh": Namespace("http://rdf.rhea-db.org/"), "rdfs": RDFS})}
    assert set(legacy.findtargetreaction(graph, "1.1.1.1")) == expected == {"10000", "10001"}
    assert legacy.findtargetreaction(graph, "9.9.9.9") == []


def test_old_assets_only_and_supported_leaves_keep_all_votes(legacy_inputs):
    name, resources, work, output, _ = legacy_inputs
    assert not (resources / "GENERATE_POOL_NEW").exists()
    scores = legacy.nonhome(name, work_dir=work, output_dir=output, resource_root=resources)
    generic = scores.loc[(scores["reaction"] == "gene1") & (scores["rhea"] == "10000")]
    assert len(generic) == 1 and generic.iloc[0]["final_score"] == 0.7
    assert "cutoff" not in set(scores["reaction"])  # strictly > 0.7
    assert "ec_only" not in set(scores["reaction"])  # original eggNOG KEGG route
    assert "partial" not in set(scores["reaction"])
    resolved = legacy.get_spe_reactions(name, work_dir=output, output_dir=output, resource_root=resources)
    assert dict(zip(resolved["rhea"], resolved["reaction"])) == {
        "10001": "gene1", "20000": "gene1", "30000": "gene1",
    }
    assert not (work / f"{name}multiple_score.xlsx").exists()
    assert not (work / f"{name}gpr_score.xlsx").exists()
    assert (output / f"{name}multiple_score.xlsx").exists()
    assert (output / f"{name}gpr_score.xlsx").exists()


def test_float_rhea_ids_and_independent_parent_copies():
    frame = pd.DataFrame([
        {"reaction": "g1", "rhea": 10000.0, "anno_A": "clean", "anno_B": None, "anno_C": "deepec", "anno_D": None},
        {"reaction": "g1", "rhea": "20000.0", "anno_D": "plmsearch"},
        {"reaction": "g1", "rhea": "30000", "anno_D": "plmsearch"},
        {"reaction": "g2", "rhea": "40000", "anno_D": "plmsearch"},
    ])
    result = legacy.filter_plmsearch(frame, {"10000": ["20000", "30000", "40000"]}, 0.4)
    assert result["rhea"].tolist() == ["20000", "30000"]
    assert result["final_score"].tolist() == [0.5, 0.5]
    assert result["reaction"].tolist() == ["g1", "g1"]
    assert frame.iloc[0]["rhea"] == 10000.0  # caller frame unchanged


def test_no_supported_leaf_returns_valid_empty_output(legacy_inputs):
    name, resources, _, output, _ = legacy_inputs
    frame = pd.DataFrame([{"reaction": "g1", "rhea": 10000.0, "anno_A": "clean", "anno_B": None,
                           "anno_C": "deepec", "anno_D": None, "final_score": 0.4}])
    frame.to_excel(output / f"{name}multiple_score.xlsx", index=False)
    result = legacy.get_spe_reactions(name, work_dir=output, resource_root=resources)
    assert result.empty and list(result.columns) == ["rhea", "reaction"]
    written = pd.read_excel(output / f"{name}gpr_score.xlsx")
    assert written.empty and {"rhea", "reaction"} <= set(written.columns)


def test_all_annotation_sources_empty_or_disabled(legacy_inputs):
    name, resources, work, output, _ = legacy_inputs
    pd.DataFrame(columns=[0, "EC", "Rec"]).to_excel(work / f"eggec2{name}.xlsx")
    scores = legacy.nonhome(name, False, False, False, work_dir=work, output_dir=output, resource_root=resources)
    assert scores.empty and list(scores.columns) == legacy.SCORE_COLUMNS
    assert legacy.get_spe_reactions(name, work_dir=output, resource_root=resources).empty


def test_module_has_no_experimental_database_imports():
    tree = ast.parse(Path(legacy.__file__).read_text())
    imports = {node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)}
    imports.update(alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names)
    assert not any(any(blocked in value for blocked in ("resource_loader", "annotation_database", "ec_specificity"))
                   for value in imports)


def test_legacy_weight_contract():
    assert legacy.calculate_score({"anno_A": "clean", "anno_B": "eggnog_mapper", "anno_C": "deepec", "anno_D": "plmsearch"}) == 0.8
    assert legacy.calculate_score({"anno_B": "eggnog_mapper"}) == 0.3
