#!/usr/bin/env python3
"""GSE75379 sayım tablosunu DESeq2 için hazırlar.

Girdi : yazarların sayım dosyası (GSE75379_count.txt.gz) + örnek bilgisi (metadata CSV)
Çıktı : <out>/counts.csv  (satır: gen, sütun: örnek)
        <out>/coldata.csv (örnek bilgisi: inek, durum, çeyrek)
Sayılar tam sayı değilse (örn. RPKM) DESeq2 için uygun olmadığından durur.
"""
import argparse, os, re, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--counts', required=True)
ap.add_argument('--meta', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--cows', default='', help='virgülle ayrılmış inek numaraları (boşsa hepsi)')
a = ap.parse_args()

meta = pd.read_csv(a.meta)
if a.cows:
    keep = [int(x) for x in a.cows.split(',')]
    meta = meta[meta.cow.isin(keep)]
print(f'{len(meta)} örnek, {meta.cow.nunique()} inek seçildi')

df = pd.read_csv(a.counts, sep=None, engine='python', index_col=0)
print('sayım tablosu:', df.shape[0], 'satır x', df.shape[1], 'sütun')

def norm(s):
    return re.sub(r'[^A-Z0-9]', '-', str(s).upper())

cols = {norm(c): c for c in df.columns}
chosen = {}
for _, r in meta.iterrows():
    found = None
    for key in (r['title'], r['GSM'], r['Run']):
        if norm(key) in cols:
            found = cols[norm(key)]; break
    if found is None:
        hits = [c for n, c in cols.items() if norm(r['title']) in n.split('--')[-1:] or n.endswith('-' + norm(r['title']))]
        if len(hits) == 1:
            found = hits[0]
    if found is None:
        print(f"HATA: '{r['title']}' ({r['GSM']}) sayım dosyasında bulunamadı.")
        print('Dosyadaki ilk sütun adları:', list(df.columns[:40]))
        sys.exit(1)
    chosen[r['title']] = found

X = df[[chosen[t] for t in meta.title]].copy()
X.columns = list(meta.title)
X = X.apply(pd.to_numeric, errors='coerce')
if X.isna().any().any():
    print('HATA: sayım tablosunda sayı olmayan değerler var'); sys.exit(1)
if (X < 0).any().any():
    print('HATA: negatif değerler var'); sys.exit(1)

dev = (X - X.round()).abs().max().max()
print('örnek başına toplam okuma (milyon):')
print((X.sum() / 1e6).round(2).to_string())
if dev > 1e-6:
    print(f'HATA: değerler tam sayı değil (en büyük sapma {dev:.3f}). Bu tablo RPKM/normalize olabilir; '
          'DESeq2 ham sayım ister. Salmon yoluna dönmemiz gerekir.')
    sys.exit(2)
X = X.round().astype(int)

os.makedirs(a.out, exist_ok=True)
X.to_csv(os.path.join(a.out, 'counts.csv'))
cd = meta[['title', 'cow', 'condition', 'quarter', 'time_h']].rename(columns={'title': 'sample'})
cd.to_csv(os.path.join(a.out, 'coldata.csv'), index=False)
print('yazıldı:', a.out, '-> counts.csv, coldata.csv')
