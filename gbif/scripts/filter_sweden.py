#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Select specific columns, filter Swedish records, remove non-WGS84
records and remove duplicate entries.
"""

import pandas as pd
import numpy as np


DTYPES = snakemake.params.dtypes
DISTINCT_COLUMN_SET = snakemake.params.distinct_column_set


def remove_whitespace(s):
    """Remove whitespace in string"""
    try:
        s = "".join(s.split())
    except AttributeError:
        pass
    return s


gbif = pd.read_table(
    snakemake.input[0], usecols=list(DTYPES.keys()), dtype=DTYPES)

include_mask = gbif.country.eq('Sweden')

if 'geodeticDatum' in gbif.columns:
    geodetic_datum = gbif.geodeticDatum.str.lower().apply(remove_whitespace)
    include_mask = (
        include_mask &
        (geodetic_datum.isnull() | geodetic_datum.str.contains('wgs84')))
else:
    gbif['geodeticDatum'] = np.nan

if 'coordinateUncertaintyInMeters' not in gbif.columns:
    gbif['coordinateUncertaintyInMeters'] = np.nan

gbif_filtered = (
    gbif[include_mask]
    .dropna(subset=['decimalLatitude', 'decimalLongitude'])
    .drop_duplicates(DISTINCT_COLUMN_SET)
)

gbif_filtered.to_csv(snakemake.output[0], sep='\t', index=False)
