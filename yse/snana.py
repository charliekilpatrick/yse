#!/usr/bin/env python3
"""Parse YSE (+ZTF) light curves in SNANA-style multi-epoch text format."""

# Original reader: Konstantin Malanchev & Patrick Aleo; path handling and bugfixes here.

from __future__ import annotations

import dataclasses
import os
from pathlib import Path

import pandas as pd
from astropy.table import Table

REDSHIFT_UNKNOWN = -99.0

# Passbands to drop when keep_ztf=False. YSE DR1 often uses X/Y for survey bands; extend this if
# your build tags ZTF-specific rows with dedicated passband names.
ZTF_BANDS: frozenset[str] = frozenset()


@dataclasses.dataclass
class Observation:
    MJD: float
    PASSBAND: str
    FLUX: float
    FLUXERR: float
    MAG: float
    MAGERR: float
    PHOTFLAG: str


def read_YSE_ZTF_snana_dir(
    dir_name: str | os.PathLike[str],
    keep_ztf: bool = True,
) -> tuple[list[str], list[dict], list[pd.DataFrame]]:
    """Load every ``*.dat`` SNANA file in *dir_name*.

    Parameters
    ----------
    dir_name
        Directory containing ``*.snana.dat`` (or similar) files—not the process working directory.
    keep_ztf
        If False, rows whose ``PASSBAND`` is in :data:`ZTF_BANDS` are removed. If ``ZTF_BANDS`` is
        empty, no rows are removed (configure bands for your data release).
    """
    root = Path(dir_name).expanduser().resolve()
    dat_files = sorted(p.name for p in root.glob("*.dat"))
    snid_list: list[str] = []
    meta_list: list[dict] = []
    frames: list[pd.DataFrame] = []

    for fname in dat_files:
        file_path = root / fname
        meta: dict = {}
        lc: list[Observation] = []
        with file_path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                try:
                    if line.startswith("SNID: "):
                        _, snid = line.split()
                        meta["object_id"] = snid
                        meta["original_object_id"] = snid
                    if line.startswith("RA: "):
                        _, ra, _ = line.split()
                        meta["ra"] = float(ra)
                    if line.startswith("DECL: "):
                        _, decl, _ = line.split()
                        meta["dec"] = float(decl)
                    if line.startswith("MWEBV: "):
                        _, mwebv, _, _mwebv_error, *_ = line.split()
                        meta["mwebv"] = float(mwebv)
                    if line.startswith("REDSHIFT_FINAL: "):
                        try:
                            _, redshift, _, _redshift_error, _z_type, _z_frame = line.split()
                        except ValueError:
                            redshift = -99
                            _redshift_error = -99
                            _z_type = "NaN"
                            _z_frame = "HELIO"
                        meta["redshift"] = float(redshift)
                        meta["redshift_err"] = float(_redshift_error)
                        meta["redshift_type"] = str(_z_type.split("(")[1].split(",")[0])
                        meta["redshift_frame"] = str(_z_frame.split(")")[0])
                    if line.startswith("PHOTO_Z: "):
                        try:
                            _, photoz, _, _photoz_error, _, _ = line.split()
                        except ValueError:
                            photoz = -99
                            _photoz_error = -99
                        meta["photo_z"] = float(photoz)
                        meta["photoz_err"] = float(_photoz_error)
                    if line.startswith("SN_OFFSET_TO_VETTED_HOST_GALAXY_CENTER: "):
                        try:
                            _, sn_offset, _ = line.split()
                        except ValueError:
                            sn_offset = -99.000
                        meta["sn_offset"] = float(sn_offset)
                    if line.startswith("VETTED_HOST_GALAXY_NAME: "):
                        try:
                            _, host_gal_name_cat, host_gal_name_id, host_gal_name_source = line.split()
                            host_gal_name = str(host_gal_name_cat) + " " + str(host_gal_name_id)
                        except ValueError:
                            host_gal_name = "None (or error)"
                            host_gal_name_source = "(SIMBAD,SDSS)"
                        meta["host_gal_name"] = host_gal_name
                        meta["host_gal_name_source"] = str(host_gal_name_source)
                    if line.startswith("VETTED_HOST_GALAXY_REDSHIFT: "):
                        try:
                            _, hostz, _, _hostz_error, _hostz_type, _hostz_frame = line.split()
                        except ValueError:
                            hostz = -99
                            _hostz_error = -99
                            _hostz_type = "NaN"
                            _hostz_frame = "HELIO"
                        meta["host_gal_z"] = float(hostz)
                        meta["host_gal_z_err"] = float(_hostz_error)
                        meta["host_gal_z_type"] = str(_hostz_type.split("(")[1].split(",")[0])
                        meta["host_gal_z_frame"] = str(_hostz_frame.split(")")[0])
                    if line.startswith("SEARCH_PEAKMJD: "):
                        _, pkmjd = line.split()
                        meta["peakmjd"] = float(pkmjd)
                    if line.startswith("HOST_LOGMASS: "):
                        _, host_logmass, _, host_logmass_error = line.split()
                        meta["host_logmass"] = float(host_logmass)
                    if line.startswith("PEAK_ABS_MAG: "):
                        _, pkabsmag = line.split()
                        try:
                            meta["peak_abs_mag"] = float(pkabsmag)
                        except ValueError:
                            meta["peak_abs_mag"] = str(pkabsmag)
                    if line.startswith("SPEC_CLASS: "):
                        try:
                            _, sn, spec_subtype = line.split()
                            meta["transient_spec_class"] = str(sn + spec_subtype)
                        except ValueError:
                            _, spec_subtype = line.split()
                            meta["transient_spec_class"] = str(spec_subtype)
                    if line.startswith("SPEC_CLASS_BROAD: "):
                        try:
                            _, sn, subtype = line.split()
                            meta["spectype_3class"] = str(sn + subtype)
                        except ValueError:
                            _, subtype = line.split()
                            meta["spectype_3class"] = str(subtype)
                    if line.startswith("PARSNIP_PRED: "):
                        try:
                            _, sn, p_pred = line.split()
                            meta["parsnip_pred_class"] = str(sn + p_pred)
                        except ValueError:
                            _, p_pred = line.split()
                            meta["parsnip_pred_class"] = str(p_pred)
                    if line.startswith("PARSNIP_CONF: "):
                        _, p_conf = line.split()
                        meta["parsnip_pred_conf"] = str(p_conf)
                    if line.startswith("PARSNIP_S1: "):
                        _, s1, _, s1_error = line.split()
                        try:
                            meta["parsnip_s1"] = float(s1)
                            meta["parsnip_s1_err"] = float(s1_error)
                        except ValueError:
                            meta["parsnip_s1"] = str(s1)
                            meta["parsnip_s1_err"] = str(s1_error)
                    if line.startswith("PARSNIP_S2: "):
                        _, s2, _, s2_error = line.split()
                        try:
                            meta["parsnip_s2"] = float(s2)
                            meta["parsnip_s2_err"] = float(s2_error)
                        except ValueError:
                            meta["parsnip_s2"] = str(s2)
                            meta["parsnip_s2_err"] = str(s2_error)
                    if line.startswith("PARSNIP_S3: "):
                        _, s3, _, s3_error = line.split()
                        try:
                            meta["parsnip_s3"] = float(s3)
                            meta["parsnip_s3_err"] = float(s3_error)
                        except ValueError:
                            meta["parsnip_s3"] = str(s3)
                            meta["parsnip_s3_err"] = str(s3_error)
                    if line.startswith("SUPERPHOT_PRED: "):
                        try:
                            _, sn, s_pred = line.split()
                            meta["superphot_pred_class"] = str(sn + s_pred)
                        except ValueError:
                            _, s_pred = line.split()
                            meta["superphot_pred_class"] = str(s_pred)
                    if line.startswith("SUPERPHOT_CONF: "):
                        _, s_conf = line.split()
                        meta["superphot_pred_conf"] = str(s_conf)
                    if line.startswith("SUPERRAENN_PRED: "):
                        try:
                            _, sn, sr_pred = line.split()
                            meta["superraenn_pred_class"] = str(sn + sr_pred)
                        except ValueError:
                            _, sr_pred = line.split()
                            meta["superraenn_pred_class"] = str(sr_pred)
                    if line.startswith("SUPERRAENN_CONF: "):
                        _, sr_conf = line.split()
                        meta["superraenn_pred_conf"] = str(sr_conf)
                    if line.startswith("SET_ZTF_FP: "):
                        _, ztf_fp = line.split()
                        try:
                            meta["ztf_zeropoint"] = float(ztf_fp)
                        except ValueError:
                            meta["ztf_zeropoint"] = str(ztf_fp)
                    if line.startswith("PEAK_SNR: "):
                        _, pkSNR = line.split()
                        meta["peakSNR"] = float(pkSNR)
                    if line.startswith("MAX_MJD_GAP(days): "):
                        _, max_mjd_gap = line.split()
                        meta["max_mjd_gap"] = float(max_mjd_gap)
                    if line.startswith("NOBS_BEFORE_PEAK: "):
                        _, nobs_before_peak = line.split()
                        meta["nobs_before_peak"] = int(nobs_before_peak)
                    if line.startswith("NOBS_TO_PEAK: "):
                        _, nobs_to_peak = line.split()
                        meta["nobs_to_peak"] = int(nobs_to_peak)
                    if line.startswith("NOBS_AFTER_PEAK: "):
                        _, nobs_after_peak = line.split()
                        meta["nobs_after_peak"] = int(nobs_after_peak)
                    if line.startswith("SEARCH_PEAKMAG: "):
                        _, pkmag = line.split()
                        meta["peakmag"] = float(pkmag)
                    if line.startswith("SEARCH_PEAKFLT: "):
                        _, pkflt = line.split()
                        meta["peakflt"] = str(pkflt)
                    if line.startswith("PEAKMAG_YSE-r/ZTF-r(Y): "):
                        _, pkmag_rY = line.split()
                        meta["peakmag_rY"] = float(pkmag_rY)
                    if line.startswith("PEAKFLT_YSE-r/ZTF-r(Y): "):
                        _, pkflt_rY = line.split()
                        meta["peakflt_rY"] = str(pkflt_rY)
                    if line.startswith("FILTERS: "):
                        _, pbs = line.split()
                        meta["passbands"] = str(pbs)
                    if line.startswith("NOBS_wZTF: ") or line.startswith("NOBS_AFTER_MASK: "):
                        _, desired_nobs = line.split()
                        meta["num_points"] = int(desired_nobs)
                        continue
                except ValueError as e:
                    print(e)
                    print(meta.get("object_id", fname))
                    raise

                if not line.startswith("OBS: "):
                    continue

                _obs, mjd, flt, _field, fluxcal, fluxcalerr, mag, magerr, _flag = line.split()
                lc.append(
                    Observation(
                        MJD=float(mjd),
                        PASSBAND=str(flt),
                        FLUX=float(fluxcal),
                        FLUXERR=float(fluxcalerr),
                        MAG=float(mag),
                        MAGERR=float(magerr),
                        PHOTFLAG=str(_flag),
                    )
                )

        meta.setdefault("mwebv", 0.0)
        if len(lc) != meta["num_points"]:
            raise AssertionError(f"{fname}: {len(lc)} obs != num_points {meta['num_points']}")
        rows = [dataclasses.asdict(o) for o in lc if keep_ztf or o.PASSBAND not in ZTF_BANDS]
        table = Table(rows)
        snid_list.append(meta["object_id"])
        meta_list.append(meta)
        frames.append(table.to_pandas())

    return snid_list, meta_list, frames
