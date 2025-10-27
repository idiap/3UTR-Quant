# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------

import argparse
import os

import yaml
from scanalysis.loader import read_bustools_counts


def load_mtx(path_to_matrices, treatments):
    dfs_matrices = {}
    for treatment, ids in treatments.items():
        dfs_treatment = {}
        for id in ids:
            dfs_treatment[id] = read_bustools_counts(
                path=os.path.join(path_to_matrices, id)
            )
        dfs_matrices[treatment] = dfs_treatment

    return dfs_matrices


def main():
    FLAGS = argparse.ArgumentParser(description="Length analysis for bulk RNAseq data.")
    FLAGS.add_argument("--config_file", help="Location of config file")

    args = FLAGS.parse_args()
    with open(args.config_file, "r") as stream:
        config = yaml.safe_load(stream)

    dfs_matrices = load_mtx(
        path_to_matrices=config["path_to_matrices"], treatments=config["treatments"]
    )

    # TODO: continue
    if len(dfs_matrices) == 0:
        print("Oh matrix is empty.")


if __name__ == "__main__":
    main()
