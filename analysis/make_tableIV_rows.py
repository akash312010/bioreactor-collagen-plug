"""Generate Table IV rows from headline_numbers.json at three significant figures (half-up).
Usage: python make_tableIV_rows.py <ROOT>/headline_numbers.json > tableIV_rows.tex"""
import json, math, sys
from decimal import Decimal, ROUND_HALF_UP
H = json.load(open(sys.argv[1]))['tracer']['runs']
def sf(x, n=3):
    e = math.floor(math.log10(abs(x))); q = Decimal(1).scaleb(e - n + 1)
    return format(Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP), 'f')
for k, v in H.items():
    print(f"{k} & {v['geom']} & {sf(v['Pe'])} & {sf(v['t_ideal_ms'])} & {sf(v['tau_ms'])} & {sf(v['t99_out'])} & {sf(v['t99_V'])}\\\\")
