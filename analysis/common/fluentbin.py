import re, numpy as np
HDR = re.compile(rb"\((30\d\d|20\d\d|3300|2300) \(([0-9a-fA-F \-]+)\)\s*\(")
def _sections(b):
    for m in HDR.finditer(b):
        yield int(m.group(1)), m.group(2).split(), m.end()
def read_cas(fn):
    b = open(fn, 'rb').read(); nodes = None; faces = {}
    for sid, h, st in _sections(b):
        if sid == 3010:
            zone, f, l, typ, nd = [int(x, 16) for x in h]
            n = l-f+1; nodes = np.frombuffer(b, '<f8', n*nd, st).reshape(n, nd)
        elif sid == 3013 or sid == 2013:
            zone, f, l, bc, typ = [int(x, 16) for x in h]
            n = l-f+1; assert typ == 2, typ
            faces[zone] = (bc, f, np.frombuffer(b, '<i4', n*4, st).reshape(n, 4).astype(np.int64))
    names = dict((int(z), (t.decode(), nm.decode())) for z, t, nm in re.findall(rb"\(39 \((\d+) ([a-z-]+) ([A-Za-z0-9_\-:.]+)", b))
    if not names:
        names = dict((int(z), (t.decode(), nm.decode())) for z, t, nm in re.findall(rb"\(45 \((\d+) ([a-z-]+) ([A-Za-z0-9_\-:.]+)", b))
    return nodes, faces, names
def read_dat(fn):
    b = open(fn, 'rb').read(); out = {}
    for sid, h, st in _sections(b):
        if sid != 3300: continue
        fid, zone, size, nt, nph, f, l = [int(x) for x in h]
        n = l-f+1; out[(fid, zone)] = (f, np.frombuffer(b, '<f8', n*size, st).reshape(n, size) if size > 1 else np.frombuffer(b, '<f8', n, st))
    return out
