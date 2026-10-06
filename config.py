"""Pipeline configuration.

Every path is derived from this file's own location, so the checkout can be
moved and re-run anywhere. Only the large inputs (L1000 matrix, GEO series
matrices) sit outside the checkout; point HCC_HF_RAW at whatever directory
holds them.

    output/        numbers written by the pipeline
    figures/       figures
    tables/        tables
    data/          small annotation inputs committed with the code
    large_inputs/  L1000 release files and discovery GEO matrices ($HCC_HF_RAW)
    downloads/     validation cohorts fetched by 01_download_geo.sh ($HCC_HF_VAL)

No absolute paths; one SEED for every stochastic step.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE

OUT = os.path.join(ROOT, 'output')
FIG = os.path.join(ROOT, 'figures')
TBL = os.path.join(ROOT, 'tables')

RAW = os.environ.get('HCC_HF_RAW', os.path.join(ROOT, 'large_inputs'))
DATA = RAW
VAL = os.environ.get('HCC_HF_VAL', os.path.join(ROOT, 'downloads'))

HUB = os.path.join(ROOT, 'data', 'repurposing_drugs_20200324.txt')
L1000 = os.path.join(RAW, 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx')

SEED = 20260925

for _d in (OUT, FIG, TBL):
    os.makedirs(_d, exist_ok=True)
