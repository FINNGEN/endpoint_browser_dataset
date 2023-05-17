import argparse
import csv
import json
import sys
from datetime import date
from math import isnan
from pathlib import Path


DATASET_NAME = "FG_pheno-r10_geno-r9"  # will be used to name output files


def main():
    log_info("hello.")
    args = init_cli()

    path_out_tree = name_output(args.output_dir / "endpoint_tree")
    path_out_info = name_output(args.output_dir / "endpoints_info")
    path_out_code_cases = name_output(args.output_dir / "endpoint_code_cases")

    log_info("gathering info for all endpoints")
    info = gather_info(
        args.definitions,
        args.basic_stats,
        args.metaresults_ukbb,
        args.metaresults_est,
        args.genetic_correlations
    )
    output_info(info, path_out_info)

    log_info("building endpoint tree from case overlap")
    tree = build_tree(
        args.correlations,
        args.subset_threshold
    )
    output_tree(
        tree,
        path_out_tree
    )

    log_info("finding N cases / registry-code / endpoint")
    code_dir = args.upset_plots / "name_tables_json"
    code_cases = find_code_cases(code_dir)
    output_code_cases(code_cases, path_out_code_cases)

    log_info("bye.")


def log_info(message):
    print(f"INFO {message}", file=sys.stderr)


def init_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--definitions",
        help="path to FinnGen case and control endpoint definitions (CSV)",
        required=True,
        type=Path
    )
    parser.add_argument(
        "-b", "--basic-stats",
        help="path to FinnGen Risteys basic stats (JSON)",
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
        "-g", "--genetic-correlations",
        help="path to the FinnGen genetic correlations (CSV)",
        required=True,
        type=Path
    )
    parser.add_argument(
        "--upset-plots",
        help="path to the FinnGen upset plots directory",
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


def name_output(base_name):
    return f"{base_name}__{DATASET_NAME}__{date.today()}.json"


def gather_info(
        path_definitions,
        path_basic_stats,
        path_meta_ukbb,
        path_meta_est,
        path_genetic_correlations
):
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

    set_ncases(endpoints, path_basic_stats)

    set_meta_analysis(endpoints, "uk_meta_analysed", path_meta_ukbb)
    set_meta_analysis(endpoints, "est_meta_analysed", path_meta_est)

    gws_hits = get_gws_hits(path_genetic_correlations)
    for endpoint in endpoints:
        hits = gws_hits.get(endpoint)
        endpoints[endpoint]["gwas_hits"] = hits

    return endpoints


def set_ncases(endpoints, path_stats):
    with open(path_stats) as fd:
        stats = json.load(fd)

    for endpoint, data in stats["stats"].items():
        n_cases = {
            "n_cases_all": data["nindivs_all"],
            "n_cases_female": data["nindivs_female"],
            "n_cases_male": data["nindivs_male"]
        }

        # Check for individual-level data
        for nn in n_cases.values():
            assert nn is None or isnan(nn) or nn == 0 or nn >= 5, f"FAIL individual-level data detected: {nn=} in {n_cases=} for {endpoint=}"

        endpoints[endpoint].update(n_cases)


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


def get_gws_hits(path):
    gws_hits = {}

    with open(path) as fd:
        reader = csv.DictReader(fd)

        for row in reader:
            endpoint = row["pheno1"]
            hits = row["n_gwsig_1"]

            gws_hits[endpoint] = hits

    return gws_hits


def build_tree(path_correlations, threshold):
    tree = []

    with open(path_correlations) as fd:
        reader = csv.DictReader(fd)

        for row in reader:
            if row["endpoint_a"] == row["endpoint_b"]:
                continue

            parent = row["endpoint_a"]
            child = row["endpoint_b"]
            ratio_shared_of_b = float(row["ratio_shared_of_b"])
            case_overlap = float(row["jaccard_index"])

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


def find_code_cases(path):
    out = {}

    data_files = path.glob("name_table_*.json")

    for fpath in data_files:
        with open(fpath) as fd:
            data = json.load(fd)

        endpoint = fpath.name.lstrip("name_table_").rstrip(".json")

        code_cases = {}
        for record in data:
            code = record["Tag"]
            ncases = record["Frequency_people"]

            assert ncases == "<5" or int(ncases) == 0  or int(ncases) >= 5, \
                f"FAIL individual-level data detected: {ncases=} for {endpoint=} in {code=}"

            code_cases[code] = ncases

        out[endpoint] = code_cases

    return out


def output_code_cases(data, path):
    with open(path, "x") as fd:
        json.dump(data, fd)


if __name__ == "__main__":
    main()
