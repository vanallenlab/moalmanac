import argparse
import collections
import difflib
import json
import pathlib
import re
import sys

import pandas as pd

VOLATILE_KEYS = {"execution_runtime", "environment"}
DATE_LINE = re.compile(r"Generated on:.*")
BASE64_IMAGE = re.compile(r"data:image/png;base64,[A-Za-z0-9+/=]+")


def list_files(folder):
    """
    Lists the names of the files directly within a folder.

    Args:
        folder (str | pathlib.Path): Folder to list.

    Returns:
        set[str]: File names, without directories.
    """
    return {p.name for p in pathlib.Path(folder).iterdir() if p.is_file()}


def row_hashes(df):
    """
    Hashes each row of a data frame, ignoring the index.

    Args:
        df (pandas.DataFrame): Data frame to hash.

    Returns:
        list[int]: One hash per row, in row order.
    """
    return pd.util.hash_pandas_object(df, index=False).tolist()


def unmatched_rows(df1, df2):
    """
    Finds the rows of each data frame that have no identical row in the other,
    counting duplicate rows separately.

    Args:
        df1 (pandas.DataFrame): First data frame.
        df2 (pandas.DataFrame): Second data frame, with the same columns.

    Returns:
        tuple[pandas.DataFrame, pandas.DataFrame]: The unmatched rows of each.
    """
    hashes1, hashes2 = row_hashes(df1), row_hashes(df2)
    remaining = collections.Counter(hashes2)
    keep1 = []
    for h in hashes1:
        if remaining[h] > 0:
            remaining[h] -= 1
            keep1.append(False)
        else:
            keep1.append(True)
    remaining = collections.Counter(hashes1)
    keep2 = []
    for h in hashes2:
        if remaining[h] > 0:
            remaining[h] -= 1
            keep2.append(False)
        else:
            keep2.append(True)
    return df1[keep1], df2[keep2]


def describe_row(row, columns):
    """
    Describes a row by its values in the first few columns.

    Args:
        row (pandas.Series): The row.
        columns (list[str]): Column names, in order.

    Returns:
        str: The row's values in up to three columns, such as "gene='BRAF', ...".
    """
    return ", ".join(f"{c}={row[c]!r}" for c in columns[:3])


def compare_tsv(path1, path2, ignore_columns):
    """
    Compares two tab-delimited files, ignoring row order and any ignored columns.
    Reports column, row count, row order, and cell differences separately.
    Rows identical in both files are set aside first; if the same number of rows
    remain on each side, they are sorted and paired to report differing cells,
    otherwise a few unmatched rows from each side are shown.

    Args:
        path1 (str | pathlib.Path): First file.
        path2 (str | pathlib.Path): Second file.
        ignore_columns (list[str]): Column names to exclude from the comparison.

    Returns:
        list[str]: Descriptions of differences; empty if the files match.
    """
    df1 = pd.read_csv(path1, sep="\t", dtype=str, keep_default_na=False)
    df2 = pd.read_csv(path2, sep="\t", dtype=str, keep_default_na=False)
    messages = []

    ignored = [c for c in ignore_columns if c in df1.columns or c in df2.columns]
    only1 = [c for c in df1.columns if c not in df2.columns and c not in ignored]
    only2 = [c for c in df2.columns if c not in df1.columns and c not in ignored]
    if only1:
        messages.append(f"columns only in folder 1: {only1}")
    if only2:
        messages.append(f"columns only in folder 2: {only2}")
    columns = [c for c in df1.columns if c in df2.columns and c not in ignored]
    if [c for c in df1.columns if c in columns] != [
        c for c in df2.columns if c in columns
    ]:
        messages.append("shared columns are in a different order")
    df1 = df1[columns].reset_index(drop=True)
    df2 = df2[columns].reset_index(drop=True)

    if df1.shape[0] != df2.shape[0]:
        messages.append(f"row count differs: {df1.shape[0]} vs. {df2.shape[0]}")
    if not columns:
        return messages

    left1, left2 = unmatched_rows(df1, df2)
    if left1.empty and left2.empty:
        if not df1.equals(df2):
            messages.append("same content, different row order")
        return messages

    messages.append(
        f"rows only in folder 1: {len(left1)}; only in folder 2: {len(left2)}",
    )
    left1 = left1.sort_values(columns, kind="stable").reset_index(drop=True)
    left2 = left2.sort_values(columns, kind="stable").reset_index(drop=True)
    if len(left1) == len(left2):
        for i in range(min(len(left1), 10)):
            row1, row2 = left1.iloc[i], left2.iloc[i]
            for column in columns:
                if row1[column] != row2[column]:
                    messages.append(
                        f"changed row ({describe_row(row1, columns)}) {column}: "
                        f"{row1[column]!r} vs. {row2[column]!r}",
                    )
        if len(left1) > 10:
            messages.append(f"{len(left1) - 10} more changed rows")
    else:
        for label, left in (("folder 1", left1), ("folder 2", left2)):
            for i in range(min(len(left), 5)):
                messages.append(
                    f"only in {label}: {describe_row(left.iloc[i], columns)}",
                )
            if len(left) > 5:
                messages.append(f"only in {label}: {len(left) - 5} more rows")
    return messages


def normalize_json(obj, roots):
    """
    Recursively removes volatile keys and replaces datasource root paths with
    a placeholder.

    Args:
        obj (typing.Any): Parsed JSON value.
        roots (list[str]): Root path strings to replace with "<root>".

    Returns:
        typing.Any: The normalized value.
    """
    if isinstance(obj, dict):
        return {
            k: normalize_json(v, roots)
            for k, v in obj.items()
            if k not in VOLATILE_KEYS
        }
    if isinstance(obj, list):
        return [normalize_json(v, roots) for v in obj]
    if isinstance(obj, str):
        for root in roots:
            if root:
                obj = obj.replace(root, "<root>")
        return obj
    return obj


def json_paths(a, b, prefix=""):
    """
    Recursively lists the paths at which two parsed JSON values differ.
    Two NaN floats are treated as equal.

    Args:
        a (typing.Any): First parsed JSON value.
        b (typing.Any): Second parsed JSON value.
        prefix (str): Path to the current values (default: "").

    Returns:
        list[str]: A description of each difference; empty if the values match.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for key in sorted(set(a) | set(b)):
            if key not in a:
                out.append(f"{prefix}/{key}: only in folder 2")
            elif key not in b:
                out.append(f"{prefix}/{key}: only in folder 1")
            else:
                out += json_paths(a[key], b[key], f"{prefix}/{key}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{prefix}: list length {len(a)} vs. {len(b)}"]
        out = []
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            out += json_paths(x, y, f"{prefix}[{i}]")
        return out
    if isinstance(a, float) and isinstance(b, float) and a != a and b != b:
        return []
    return [] if a == b else [f"{prefix}: {str(a)[:80]!r} vs. {str(b)[:80]!r}"]


def compare_json(path1, path2):
    """
    Compares two execution JSON files after normalizing volatile content.

    Args:
        path1 (str | pathlib.Path): First file.
        path2 (str | pathlib.Path): Second file.

    Returns:
        list[str]: Descriptions of differences; empty if the files match.
    """
    d1 = json.loads(pathlib.Path(path1).read_text())
    d2 = json.loads(pathlib.Path(path2).read_text())
    roots = [(d.get("input_datasources") or {}).get("root", "") for d in (d1, d2)]
    return json_paths(normalize_json(d1, roots), normalize_json(d2, roots))


def normalize_html_line(line):
    """
    Removes the generated date and embedded base64 images from a line of report HTML.

    Args:
        line (str): A line of HTML.

    Returns:
        str: The normalized line, truncated to 300 characters.
    """
    line = DATE_LINE.sub("Generated on:", line)
    return BASE64_IMAGE.sub("data:image/png;base64,<stripped>", line)[:300]


def compare_html(path1, path2):
    """
    Compares two report HTML files after normalizing dates and embedded images.

    Args:
        path1 (str | pathlib.Path): First file.
        path2 (str | pathlib.Path): Second file.

    Returns:
        list[str]: A summary line followed by up to 40 lines of unified diff;
            empty if the files match.
    """
    lines = []
    for p in (path1, path2):
        text = pathlib.Path(p).read_text(errors="replace")
        lines.append([normalize_html_line(ln) for ln in text.splitlines()])
    if lines[0] == lines[1]:
        return []
    diff = list(
        difflib.unified_diff(
            lines[0],
            lines[1],
            "folder1",
            "folder2",
            lineterm="",
            n=0,
        ),
    )
    changed = sum(1 for d in diff if d[:1] in "+-" and d[:3] not in ("+++", "---"))
    return [f"differs ({changed} lines)", *[d[:160] for d in diff[:40]]]


def compare_png(path1, path2, check_size):
    """
    Compares the byte sizes of two figures, only when check_size is True.

    Args:
        path1 (str | pathlib.Path): First file.
        path2 (str | pathlib.Path): Second file.
        check_size (bool): Whether to compare byte sizes at all.

    Returns:
        list[str]: A description of the size difference; empty if sizes match
            or the check is disabled.
    """
    # Figures are not byte-reproducible between runs; presence only by default
    if not check_size:
        return []
    s1, s2 = pathlib.Path(path1).stat().st_size, pathlib.Path(path2).stat().st_size
    return [] if s1 == s2 else [f"size differs: {s1} vs. {s2} bytes"]


def compare_folders(folder1, folder2, ignore_columns, png_size=False):
    """
    Compares every file in two output folders, using a comparison suited to each
    file type.
    Files with other extensions, such as logs, are not compared.

    Args:
        folder1 (str | pathlib.Path): First output folder.
        folder2 (str | pathlib.Path): Second output folder.
        ignore_columns (list[str]): Column names to exclude from table comparisons.
        png_size (bool): Also compare figure byte sizes (default: False).

    Returns:
        dict[str, list[str]]: Difference descriptions keyed by file name; an
            empty list means the file matches.
    """
    folder1, folder2 = pathlib.Path(folder1), pathlib.Path(folder2)
    files1, files2 = list_files(folder1), list_files(folder2)
    results = {}
    for name in sorted(files1 | files2):
        if name not in files2:
            results[name] = ["only in folder 1"]
        elif name not in files1:
            results[name] = ["only in folder 2"]
        elif name.endswith(".txt"):
            results[name] = compare_tsv(folder1 / name, folder2 / name, ignore_columns)
        elif name.endswith(".moalmanac-execution.json"):
            results[name] = compare_json(folder1 / name, folder2 / name)
        elif name.endswith(".report.html"):
            results[name] = compare_html(folder1 / name, folder2 / name)
        elif name.endswith(".png"):
            results[name] = compare_png(folder1 / name, folder2 / name, png_size)
        # example.log and other files are intentionally not compared
    return results


def print_results(results):
    """
    Prints whether each file matched, with its differences if not.

    Args:
        results (dict[str, list[str]]): Difference descriptions keyed by file name.
    """
    for name, messages in results.items():
        if not messages:
            print(f"  identical: {name}")
            continue
        print(f"  DIFFERS: {name}")
        for message in messages:
            print(f"    {message}")


def main(folder1, folder2, ignore_columns=(), png_size=False):
    """
    Compares two output folders and prints the results.

    Args:
        folder1 (str): First output folder.
        folder2 (str): Second output folder.
        ignore_columns (list[str]): Column names to exclude from table
            comparisons (default: ()).
        png_size (bool): Also compare figure byte sizes (default: False).

    Returns:
        int: 1 if any file differs, otherwise 0.
    """
    print(f"Comparing {folder1} and {folder2}")
    results = compare_folders(folder1, folder2, ignore_columns, png_size)
    print_results(results)
    return 1 if any(results.values()) else 0


def read_manifest(path):
    """
    Reads a manifest and resolves its run folders relative to the manifest's folder.

    Args:
        path (str | pathlib.Path): Manifest JSON, written by the regression run script.

    Returns:
        dict[str, pathlib.Path]: Output folder keyed by run name.
    """
    path = pathlib.Path(path)
    runs = json.loads(path.read_text())["runs"]
    return {name: path.parent / folder for name, folder in runs.items()}


def main_all(manifest1, manifest2, ignore_columns=(), png_size=False):
    """
    Compares all four example runs listed in two manifests and prints a summary.

    Args:
        manifest1 (str): First manifest JSON, written by the regression run script.
        manifest2 (str): Second manifest JSON.
        ignore_columns (list[str]): Column names to exclude from table
            comparisons (default: ()).
        png_size (bool): Also compare figure byte sizes (default: False).

    Returns:
        int: 1 if any run differs or is listed in only one manifest, otherwise 0.
    """
    runs1, runs2 = read_manifest(manifest1), read_manifest(manifest2)
    summary = []
    for name in list(runs1) + [n for n in runs2 if n not in runs1]:
        if name not in runs2:
            summary.append((name, "only in first manifest"))
            continue
        if name not in runs1:
            summary.append((name, "only in second manifest"))
            continue
        print(f"=== {name}")
        results = compare_folders(runs1[name], runs2[name], ignore_columns, png_size)
        print_results(results)
        n = sum(1 for m in results.values() if m)
        summary.append((name, "identical" if n == 0 else f"{n} file(s) differ"))
    print("\nSummary")
    for name, status in summary:
        print(f"  {name:<22} {status}")
    return 0 if all(s == "identical" for _, s in summary) else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Molecular Oncology Almanac output regression test",
        description="Compare all four runs from two manifests, or two output folders",
    )
    parser.add_argument(
        "--folder1",
        "-f1",
        required=True,
        help="Output folder, or manifest json to compare all four runs",
    )
    parser.add_argument(
        "--folder2",
        "-f2",
        required=True,
        help="Output folder, or manifest json to compare all four runs",
    )
    parser.add_argument(
        "--ignore-columns",
        nargs="*",
        default=[],
        help="Columns to ignore in TSVs",
    )
    parser.add_argument(
        "--png-size",
        action="store_true",
        help="Also compare PNG byte sizes (noisy)",
    )
    args = parser.parse_args()
    paths = [pathlib.Path(args.folder1), pathlib.Path(args.folder2)]
    for path in paths:
        if not path.exists():
            parser.error(f"{path} does not exist")
    manifests = [path.is_file() and path.suffix == ".json" for path in paths]
    if any(manifests) and not all(manifests):
        parser.error("pass either two output folders or two manifest json files")
    if all(manifests):
        code = main_all(args.folder1, args.folder2, args.ignore_columns, args.png_size)
    else:
        code = main(args.folder1, args.folder2, args.ignore_columns, args.png_size)
    sys.exit(code)
