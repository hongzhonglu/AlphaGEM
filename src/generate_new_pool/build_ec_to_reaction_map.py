import re
import pickle
from typing import Dict, List, Set, Tuple

from read_mnxref_reacs import read_mnxref_reacs


EC_FULL_REGEX = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def is_full_ec(ec: str) -> bool:
    return bool(EC_FULL_REGEX.match(ec))


def build_ec_to_reaction_map(
    ttl_path: str,
    include_generic: bool = True,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    从 MNXref.ttl 解析反应并返回两个映射：
    - ec_to_reacs_full: 完整 EC (如 1.1.1.1) -> [MNXR...]
    - ec_to_reacs_generic: 含通配/不完整 (含 n 或 - 等) -> [MNXR...]
    """
    reac_ann = read_mnxref_reacs(ttl_path)

    ec_to_full: Dict[str, Set[str]] = {}
    ec_to_generic: Dict[str, Set[str]] = {}

    for mnxr, entry in reac_ann.items():
        ecs = entry.get("ec") or []
        for ec in ecs:
            ec_norm = ec.strip()
            if not ec_norm:
                continue
            if is_full_ec(ec_norm):
                ec_to_full.setdefault(ec_norm, set()).add(mnxr)
            else:
                if include_generic:
                    ec_to_generic.setdefault(ec_norm, set()).add(mnxr)

    ec_to_reacs_full = {k: sorted(list(v)) for k, v in ec_to_full.items()}
    ec_to_reacs_generic = {k: sorted(list(v)) for k, v in ec_to_generic.items()}
    return ec_to_reacs_full, ec_to_reacs_generic


if __name__ == "__main__":
    ttl_path = "/data_available/MNXref.ttl"
    full_map, generic_map = build_ec_to_reaction_map(ttl_path, include_generic=True)

    with open("/data_available/ec_to_mnxr_full.pkl", "wb") as f:
        pickle.dump(full_map, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open("/data_available/ec_to_mnxr_generic.pkl", "wb") as f:
        pickle.dump(generic_map, f, protocol=pickle.HIGHEST_PROTOCOL)

    # 简短控制台汇总
    print(f"完整 EC 数量: {len(full_map)}，含不完整/通配 EC 数量: {len(generic_map)}")

