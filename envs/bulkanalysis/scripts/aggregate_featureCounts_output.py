
# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------
import pandas as pd
import argparse
from bulkanalysis.loader import bulk_aggregate_featureCounts_outputs




if __name__ == '__main__':
    FLAGS = argparse.ArgumentParser(description='Length analysis for bulk RNAseq data.')
    FLAGS.add_argument('-f', 
                       '--files', 
                       nargs='+',
                      help='<Required> Files to aggregate',
                      required=True)
    FLAGS.add_argument('-n', 
                       '--names', 
                       nargs='+',
                      help='<Required> Names of the files',
                      required=True)
    FLAGS.add_argument('-s', '--save_path', help='<Required> Path to save', required=True)
    args = FLAGS.parse_args()

    print(args.files)
    print(args.save_path)

    df = bulk_aggregate_featureCounts_outputs(args.files, args.names, args.save_path)


