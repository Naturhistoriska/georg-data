#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


DTYPES = snakemake.params.dtypes
DISTINCT_COLUMN_SET = snakemake.params.distinct_column_set


frame = pd.read_table(
    snakemake.input[0], usecols=list(DTYPES.keys()), dtype=DTYPES)


frame_sweden = (
    frame[frame.country.eq('Sweden')]
    .dropna(subset=['decimalLatitude', 'decimalLongitude'])
    .drop_duplicates(DISTINCT_COLUMN_SET)
)

frame_sweden.to_csv(snakemake.output[0], sep='\t', index=False)
