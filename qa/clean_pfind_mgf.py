"""Normalize MGF to the supplied pFind title convention and intensity top N."""
import argparse
import json
import logging
from pathlib import Path
import re
import time
import numpy as np
import pandas as pd


def records(path):
    headers, peaks = None, []
    with Path(path).open(encoding='utf-8-sig', buffering=1024*1024) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == 'BEGIN IONS':
                if headers is not None:
                    raise ValueError('Nested BEGIN IONS')
                headers, peaks = {}, []
            elif line == 'END IONS':
                if headers is None:
                    raise ValueError('Unexpected END IONS')
                data = np.fromstring(' '.join(peaks), sep=' ')
                if len(data) % 2:
                    raise ValueError('Invalid peak pair')
                yield headers, data.reshape(-1, 2)
                headers = None
            elif headers is not None:
                if '=' in line:
                    key, value = line.split('=', 1)
                    if key in headers:
                        raise ValueError(f'Duplicate header {key}')
                    headers[key] = value
                else:
                    peaks.append(line)
        if headers is not None:
            raise ValueError('Truncated spectrum')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--mapping', required=True)
    p.add_argument('--raw-name', required=True)
    p.add_argument('--top', type=int, default=250)
    a = p.parse_args()
    source, target = Path(a.input), Path(a.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == target.resolve():
        raise ValueError('Input and output must differ')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s',
                        handlers=[logging.FileHandler(target.with_suffix('.log'), encoding='utf-8'), logging.StreamHandler()])
    mapping = pd.read_csv(a.mapping, sep='\t').set_index('spectrum_id', drop=False)
    titles, old_titles, counts_before, counts_after, spectrum_ids = [], [], [], [], []
    t0 = time.time()
    with target.open('w', encoding='ascii', newline='\r\n', buffering=1024*1024) as f:
        for h, peaks in records(source):
            if 'SCANS' in h:
                sid = int(h['SCANS'])
            else:
                match = re.search(r'\.(\d+)\.\1\.(\d+)\.0\.dta$', h['TITLE'])
                if not match:
                    raise ValueError(f'Cannot extract spectrum ID from {h["TITLE"]}')
                sid = int(match[1])
            charge = int(h['CHARGE'].rstrip('+'))
            row = mapping.loc[sid]
            if charge != int(row.charge) or abs(float(h['PEPMASS'].split()[0])-row.precursor_mz) > 1e-7:
                raise ValueError(f'Precursor mapping mismatch at {sid}')
            good = np.isfinite(peaks).all(axis=1) & (peaks > 0).all(axis=1)
            peaks = peaks[good]
            # Full sort gives deterministic tie handling: intensity descending, then m/z.
            ix = np.lexsort((peaks[:, 0], -peaks[:, 1]))[:a.top]
            selected = peaks[ix]
            selected = selected[np.argsort(selected[:, 0], kind='stable')]
            title = f'{a.raw_name}.{sid}.{sid}.{charge}.0.dta'
            f.write(f'BEGIN IONS\nTITLE={title}\nCHARGE={charge}+\n'
                    f'RTINSECONDS={float(h["RTINSECONDS"]):.6f}\nPEPMASS={h["PEPMASS"]}\n')
            f.write(''.join(f'{mz:.6f} {it:.3f}\n' for mz, it in selected))
            f.write('END IONS\n\n')
            titles.append(title)
            old_titles.append(h['TITLE'])
            spectrum_ids.append(sid)
            counts_before.append(len(peaks))
            counts_after.append(len(selected))
            if len(titles) % 10000 == 0:
                logging.info('Cleaned %s / %s spectra', len(titles), len(mapping))
    assert len(titles) == len(mapping) and len(set(titles)) == len(titles)
    assert len(set(spectrum_ids)) == len(mapping)
    exported = mapping.loc[spectrum_ids].copy()
    exported['original_title'] = old_titles
    exported['pfind_title'] = titles
    exported['peaks_before'] = counts_before
    exported['peaks_after'] = counts_after
    exported.to_csv(target.with_suffix('.mapping.tsv'), sep='\t', index=False)
    # Independently parse the completed output and check its formatted fields and peak order.
    checked = 0
    pattern = re.compile(re.escape(a.raw_name) + r'\.(\d+)\.\1\.(\d+)\.0\.dta$')
    for h, peaks in records(target):
        match = pattern.fullmatch(h['TITLE'])
        assert match and int(match[2]) == int(h['CHARGE'].rstrip('+'))
        assert h['TITLE'] == titles[checked]
        assert len(peaks) == counts_after[checked] <= a.top
        assert (np.diff(peaks[:, 0]) >= 0).all()
        assert np.isfinite(peaks).all() and (peaks > 0).all()
        assert abs(float(h['PEPMASS'])-exported.iloc[checked].precursor_mz) < 1e-7
        checked += 1
    assert checked == len(mapping)
    summary = dict(input=str(source), output=str(target), spectra=checked,
                   title_format='raw_name.scan.scan.charge.0.dta', synthetic_unique_scan_numbers=True,
                   top_n=a.top, peak_selection='highest intensity; m/z ascending for intensity ties',
                   peak_order='m/z ascending', newline='CRLF',
                   peaks_before=sum(counts_before), peaks_after=sum(counts_after),
                   min_peaks=min(counts_after), max_peaks=max(counts_after),
                   charge_counts=exported.charge.value_counts().sort_index().to_dict(),
                   output_bytes=target.stat().st_size, elapsed_seconds=time.time()-t0,
                   validation='PASS', original_precursor_mz_charge_rt_retained=True,
                   fragment_deconvolution=False)
    target.with_suffix('.summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    logging.info('DONE %s', json.dumps(summary))


if __name__ == '__main__':
    main()
