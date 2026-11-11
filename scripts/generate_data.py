"""Generate a synthetic dataset and save it as two CSV files.

Run from the repository root:

    python -m scripts.generate_data --n-plates 6 --seed 20 --out-dir data/example --name example
"""

import argparse
from pathlib import Path

from assay_qc.generator import generate_dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-plates", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20)
    parser.add_argument("--out-dir", type=Path, default=Path("data/example"))
    parser.add_argument("--name", default="example")
    args = parser.parse_args()

    data, truth = generate_dataset(n_plates=args.n_plates, seed=args.seed)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    data_path = args.out_dir / f"{args.name}_plates.csv"
    truth_path = args.out_dir / f"{args.name}_truth.csv"
    data.to_csv(data_path, index=False)
    truth.to_csv(truth_path, index=False)

    print(f"Wrote {len(data)} wells from {data['plate_id'].nunique()} plates to {data_path}")
    print(f"Wrote {len(truth)} injected issues to {truth_path}")


if __name__ == "__main__":
    main()
