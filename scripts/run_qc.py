"""Run the QC rules on a CSV file and save the findings.

Run from the repository root:

    python -m scripts.run_qc data/example/example_plates.csv
"""

import argparse
from pathlib import Path

import pandas as pd

from assay_qc.rules import QCConfig, run_qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="CSV file with plate data")
    parser.add_argument("--out", type=Path, default=Path("findings.csv"))
    parser.add_argument("--cv-threshold", type=float, default=QCConfig.replicate_cv_pct)
    args = parser.parse_args()

    data = pd.read_csv(args.data)
    findings = run_qc(data, QCConfig(replicate_cv_pct=args.cv_threshold))
    findings.to_csv(args.out, index=False)

    print(findings[["finding_id", "plate_id", "rule", "well", "compound",
                    "concentration_uM", "observed", "threshold"]].to_string(index=False))
    print(f"\n{len(findings)} findings saved to {args.out}")


if __name__ == "__main__":
    main()
