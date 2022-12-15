import argparse
import csv
import json
from pathlib import Path


def main():
    args = init_cli()

    path_out_tree = args.output_dir / "endpoint_tree.json"
    path_out_info = args.output_dir / "endpoints_info.json"

    info = gather_info(
        args.definitions,
        args.metaresults_ukbb,
        args.metaresults_est
    )
    output_info(info, path_out_info)

    core_endpoints = get_core_endpoints(info)
    tree = build_tree(
        args.correlations,
        core_endpoints,
        args.subset_threshold
    )
    output_tree(
        tree,
        path_out_tree
    )


def init_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--definitions",
        help="path to FinnGen case and control endpoint definitions (CSV)",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-c", "--correlations",
        help="path to the FinnGen endpoint-endpoint correlations (CSV)",
        required=True,
        type=Path
    )
    default_threshold = 0.5
    parser.add_argument(
        "-t", "--subset-threshold",
        help=f"threshold value (between 0.0–1.0) to consider a subset between two endpoints (default: {default_threshold})",
        default=default_threshold,
        type=float
    )
    parser.add_argument(
        "-u", "--metaresults-ukbb",
        help="path to the FinnGen+UKBB metaresults (TSV)",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-e", "--metaresults-est",
        help="path to the FinnGen+EST metaresults (TSV)",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-o", "--output-dir",
        help="path to output directory",
        required=True,
        type=Path
    )

    args = parser.parse_args()

    return args


def gather_info(path_definitions, path_meta_ukbb, path_meta_est):
    endpoints = {}

    with open(path_definitions) as fd:
        reader = csv.DictReader(fd)

        for row in reader:
            endpoint = row["NAME"]

            endpoints[endpoint] = {
                "endpoint": endpoint,
                "longname": row["LONGNAME"],
                "core_endpoint": row["CORE_ENDPOINTS"],

                # placeholders, will be set from other files later on
                "uk_meta_analysed": "no",
                "est_meta_analysed": "no"
            }

    set_meta_analysis(endpoints, "uk_meta_analysed", path_meta_ukbb)
    set_meta_analysis(endpoints, "est_meta_analysed", path_meta_est)

    return endpoints


def set_meta_analysis(endpoints, name_meta, path_meta):
    with open(path_meta) as fd:
        reader = csv.DictReader(fd, dialect=csv.excel_tab)

        for row in reader:
            endpoint = row["phenocode"]
            endpoints[endpoint][name_meta] = "yes"


def output_info(data, path):
    with open(path, "x") as fd:
        values = list(data.values())
        json.dump(values, fd)


def get_core_endpoints(info):
    cores = set()

    for endpoint, data in info.items():
        if data["core_endpoint"] == "yes":
            cores.add(endpoint)

    return cores


def build_tree(path_correlations, core_endpoints, threshold):
    tree = []

    with open(path_correlations) as fd:
        reader = csv.DictReader(fd)

        for row in reader:
            if row["endpoint_a"] == row["endpoint_b"]:
                continue

            if not row["endpoint_a"] in core_endpoints:
                continue

            if not row["endpoint_b"] in core_endpoints:
                continue

            parent = row["endpoint_a"]
            child = row["endpoint_b"]
            ratio_shared_of_b = float(row["ratio_shared_of_b"])
            case_overlap = float(row["case_overlap_percent"])

            # Check if B is a complete subset of A
            if ratio_shared_of_b == 1.0:
                tree.append({
                    "parent": parent,
                    "child": child,
                    "subsets": True,
                    "case_overlap": case_overlap
                })

            # Check if B is a "partial subset" of A
            elif case_overlap >= threshold:
                tree.append({
                    "parent": parent,
                    "child": child,
                    "subsets": False,
                    "case_overlap": case_overlap
                })

    return tree


def output_tree(tree, path):
    with open(path, "x") as fd:
        json.dump(tree, fd)


if __name__ == "__main__":
    main()
