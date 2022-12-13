import argparse
import csv
import json
from pathlib import Path


def main():
    args = init_cli()

    path_out_tree = args.output_dir / "endpoint_tree.json"

    endpoint_tree = build_tree(args.definitions)
    output_tree(endpoint_tree, path_out_tree)


def init_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--definitions",
        help="path to FinnGen case and control endpoint definitions (CSV)",
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


if __name__ == "__main__":
    main()
