import argparse
import csv
import json
from pathlib import Path


def main():
    args = init_cli()

    path_out_tree = args.output_dir / "endpoint_tree.json"
    path_out_info = args.output_dir / "endpoints_info.json"

    endpoint_tree = build_tree(args.definitions)
    output_tree(endpoint_tree, path_out_tree)

    info = gather_info(
        args.definitions,
        args.metaresults_ukbb,
        args.metaresults_est
    )
    output_info(info, path_out_info)


def init_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--definitions",
        help="path to FinnGen case and control endpoint definitions (CSV)",
        required=True,
        type=Path
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


def build_tree(path_definitions):
    tree = {}

    def empty():
        return {
            "parents": set(),
            "children": set()
        }

    with open(path_definitions) as fd:
        reader = csv.DictReader(fd)

        for row in reader:
            endpoint = row["NAME"]
            include = split(row["INCLUDE"], "|")

            data = tree.get(endpoint, empty())

            updated_children = include.union(data["children"])

            # Go through each included endpoint to mark the current one as parent
            for child in include:
                child_data = tree.get(child, empty())
                updated_parents = set([endpoint]).union(child_data["parents"])
                tree[child] = {
                    "parents": updated_parents,
                    "children": child_data["children"]
                }


            tree[endpoint] = {
                "parents": data["parents"],
                "children": updated_children
            }

    return tree


def split(text, delimiter):
    """Custom text splitting returning empty set instead of [""]."""
    elements = text.split(delimiter)

    if elements == [""]:
        return set()
    else:
        return set(elements)


def output_tree(tree, path):
    out = {}

    for endpoint, data in tree.items():
        parents = list(data["parents"])
        children = list(data["children"])
        out[endpoint] = {
            "parents": parents,
            "children": children
        }

    with open(path, "x") as fd:
        json.dump(out, fd)


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


if __name__ == "__main__":
    main()
