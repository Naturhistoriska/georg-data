#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

dtypes = {
    'District': str,
    'Province': str,
    'Country': str,
    'Län': str,
    'Kommun': str,
    'Latitude': float,
    'Longitude': float,
    'Continent': str,
    'ID': int,
    'xmax': 'Int64',
    'xmin': 'Int64',
    'ymax': 'Int64',
    'ymin': 'Int64',
    'precision': 'Int64',
    'RT90centrE': 'Int64',
    'RT90centrN': 'Int64',
    'Code': str,
}

addendum_vh = (
    pd.read_table(snakemake.input[0], dtype=dtypes)
    .dropna(subset=['District', 'Latitude', 'Longitude'])
)

addendum_json_vh = addendum_vh.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False,), axis=1)

addendum_georg = pd.DataFrame({
    'coordinateUncertaintyInMeters': addendum_vh.precision,
    'locationDisplayLabel': addendum_vh.District})
addendum_json_georg = addendum_georg.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False,), axis=1)

pelias_frame = pd.DataFrame(index=addendum_vh.index)
pelias_frame = pelias_frame.assign(
    source='swe-virtual-herbarium',
    layer='socken',
    lat=addendum_vh.Latitude,
    lon=addendum_vh.Longitude,
    name=addendum_vh.District,
    addendum_json_vh=addendum_json_vh,
    addendum_json_georg=addendum_json_georg
)
pelias_frame.to_csv(snakemake.output[0], index=False)
