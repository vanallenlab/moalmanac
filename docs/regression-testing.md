Before opening a Pull Request that should not change outputs, run the four example executions before and after your changes and compare them. All commands run from inside `moalmanac/`.

Run the following before making any code changes. Four runs of moalmanac are executed with `run_4x_for_output_regression_test.py`, using all inputs ("full") and no inputs ("empty"), each with all function toggles on and then off.

```bash
python run_4x_for_output_regression_test.py --label pre
```

After making changes, repeat this with `--label post`.

```bash
python run_4x_for_output_regression_test.py --label post
```

Next, compare all four pairs with `run_output_regression_test.py`. An exit code of 1 is returned if anything differs.

```bash
python run_output_regression_test.py \
    -f1 <date>-regression-pre.json \
    -f2 <date>-regression-post.json
```

## Generating outputs (`run_4x_for_output_regression_test.py`)

All arguments are optional.

Each run writes to `<output-root>/<date>-example-outputs-<run>-<label>/`, where `<run>` is one of `full-all-enabled`, `empty-all-enabled`, `full-all-disabled`, or `empty-all-disabled`. The script also writes a manifest, `<output-root>/<date>-regression-<label>.json`, which lists the four output folders relative to the manifest, so the manifest can be passed to `run_output_regression_test.py` from any working directory.

Optional arguments:

| Argument | Type | Default | Description |
|---|---|---|---|
| `--label` | `pre` or `post` | `pre` | Label appended to each output folder and the manifest name |
| `--output-root` | path | `.` | Folder in which the output folders and manifest are created |
| `--date` | `YYYY-MM-DD` | today | Date prefixed to each output folder and the manifest name |
| `--overwrite` | flag | off | Delete existing output folders for this date and label before running. Without it, the script stops if any output folder already contains files, because files left from an earlier run would hide outputs that are no longer written |

## Comparing outputs (`run_output_regression_test.py`)

Pass two manifests to compare all four runs, or two output folders to compare a single run.

Required arguments:

| Argument | Type | Description |
| --- | --- | --- |
| `--folder1`, `-f1` | path | Manifest from the pre run, or an output folder |
| `--folder2`, `-f2` | path | Manifest from the post run, or an output folder |

Optional arguments:

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--ignore-columns` | one or more column names | none | TSV columns to leave out of the comparison, such as columns your change adds on purpose |
| `--png-size` | flag | off | Also compare figure byte sizes. Figures are not byte-reproducible between runs, so expect false positives |

`run_output_regression_test.py` performs different comparisons depending on the type of file:

- `*.txt` tables are compared by content, ignoring row order. Rows identical in both files are set aside; the remaining rows are paired to show changed cells, or listed if their counts differ. Row-order-only changes, row count changes and column differences are reported separately.
- `*.moalmanac-execution.json` is deep-diffed, ignoring runtime, environment and the datasource root path.
- `*.report.html` is compared after removing the "Generated on" date and base64 image payloads; a capped unified diff is printed.
- `*.png` figures are compared by presence only by default, because they are not byte-reproducible between runs.
- `*.log` files are not compared.

Expected differences (for example a version bump in `config.ini`) should be listed in the Pull Request notes.
