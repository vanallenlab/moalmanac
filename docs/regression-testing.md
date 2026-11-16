Before opening a Pull Request that should not change outputs, run the four example executions before and after your changes and compare them. All commands run from inside `moalmanac/`.

```bash
# 1. On the base code: four runs (full/empty inputs x function toggles on/off)
python run_4x_for_output_regression_test.py --label pre

# 2. Make changes, then repeat
python run_4x_for_output_regression_test.py --label post

# 3. Compare all four pairs; exit code is 1 if anything differs
python run_output_regression_test.py --all \
    -f1 <date>-regression-pre.json \
    -f2 <date>-regression-post.json
```

Options:

- `run_4x_for_output_regression_test.py`: 
	- `--label {pre,post}`
	- `--output-root DIR` (default `.`)
	- `--date YYYY-MM-DD`. 
	- Writes `<root>/<date>-regression-<label>.json`, a manifest listing the four output folders.
- `run_output_regression_test.py`: 
	- compare a single pair of folders with `-f1 FOLDER -f2 FOLDER`
	- `--ignore-columns COL ...` ignores expected new TSV columns
	- `--png-size` also compares figure byte sizes.

What is compared:

- `*.txt` tables are sorted and compared by content; row-order-only changes, row count changes, column differences and per-cell differences are reported separately.
- `*.moalmanac-execution.json` is deep-diffed, ignoring runtime, environment and the datasource root path.
- `*.report.html` is compared after removing the "Generated on" date and base64 image payloads; a capped unified diff is printed.
- `*.png` figures are compared by presence only by default, because they are not byte-reproducible between runs.
- `*.log` files are not compared.

Expected differences (for example a version bump in `config.ini`) should be listed in the Pull Request notes.
