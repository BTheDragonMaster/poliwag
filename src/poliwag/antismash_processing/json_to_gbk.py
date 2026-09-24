#!/usr/bin/env python3
"""
Convert antiSMASH JSON results file(s) into per-region GenBank (.gbk)
files, mirroring the naming convention and content of the files
antiSMASH itself writes to its "regions" output.
"""
import argparse
import datetime
import json
import re
import warnings
from pathlib import Path

from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature, Reference
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

LOCATION_PART_RE = re.compile(r"([<>]?)(\d+)")


def parse_simple_location(text, strand_required=True):
    """Parse a Biopython-style location fragment such as '[0:63029](+)'
    or '[0:63029]' (used in references, no strand) into a FeatureLocation.
    """
    text = text.strip()
    m = re.match(r"^\[([<>]?\d+):([<>]?\d+)\](?:\(([+-])\))?$", text)
    if not m:
        raise ValueError(f"Unrecognised location fragment: {text!r}")
    start_txt, end_txt, strand_txt = m.groups()

    def to_pos(txt, is_start):
        fuzzy, num = LOCATION_PART_RE.match(txt).groups()
        num = int(num)
        if fuzzy == "<":
            from Bio.SeqFeature import BeforePosition
            return BeforePosition(num)
        if fuzzy == ">":
            from Bio.SeqFeature import AfterPosition
            return AfterPosition(num)
        return num

    start = to_pos(start_txt, True)
    end = to_pos(end_txt, False)
    strand = {"+": 1, "-": -1}.get(strand_txt, None)
    return FeatureLocation(start, end, strand)


def parse_location(text):
    """Parse a full antiSMASH/Biopython location string, including
    compound 'join{...}' / 'order{...}' locations."""
    text = text.strip()
    m = re.match(r"^(join|order)\{(.+)\}$", text)
    if m:
        operator, inner = m.groups()
        parts = [parse_simple_location(p) for p in inner.split(", ")]
        return CompoundLocation(parts, operator=operator.lower())
    return parse_simple_location(text)


def location_bounds(location):
    """Return (start, end) ints for a (possibly compound) location."""
    return int(location.start), int(location.end)


def build_features(raw_features):
    features = []
    for rf in raw_features:
        loc = parse_location(rf["location"])
        qualifiers = {}
        for key, values in rf["qualifiers"].items():
            # Preserve bare "/flag" qualifiers (JSON null -> Python None)
            if values is None:
                qualifiers[key] = [None]
            else:
                qualifiers[key] = list(values)
        feat = SeqFeature(loc, type=rf["type"], qualifiers=qualifiers)
        features.append(feat)
    return features


def build_references(raw_references):
    refs = []
    for rr in raw_references:
        ref = Reference()
        ref.authors = rr.get("authors", "")
        ref.consrtm = rr.get("consrtm", "")
        ref.title = rr.get("title", "")
        ref.journal = rr.get("journal", "")
        ref.medline_id = rr.get("medline_id", "")
        ref.pubmed_id = rr.get("pubmed_id", "")
        ref.comment = rr.get("comment", "")
        locs = []
        for loc_txt in rr.get("location", []):
            locs.append(parse_simple_location(loc_txt))
        ref.location = locs
        refs.append(ref)
    return refs


def build_full_record(raw_record):
    """Build a full-length Bio.SeqRecord from one 'records[]' entry of the
    antiSMASH JSON, keeping original (whole-record) coordinates."""
    seq = Seq(raw_record["seq"]["data"])
    rec = SeqRecord(
        seq,
        id=raw_record["id"],
        name=raw_record["name"],
        description=raw_record["description"],
    )
    rec.annotations = dict(raw_record.get("annotations", {}))
    if "references" in rec.annotations:
        rec.annotations["references"] = build_references(rec.annotations["references"])
    rec.dbxrefs = list(raw_record.get("dbxrefs", []))
    rec.features = build_features(raw_record.get("features", []))
    return rec


def slice_region(full_record, area_start, area_end, antismash_version, run_date):
    """Slice out [area_start:area_end) from full_record, keeping only
    features that are fully contained in that window, and shifting all
    coordinates so the region starts at 1."""
    sub_seq = full_record.seq[area_start:area_end]

    sub_features = []
    for feat in full_record.features:
        f_start, f_end = location_bounds(feat.location)
        if f_start >= area_start and f_end <= area_end:
            new_feat = SeqFeature(
                feat.location + (-area_start),
                type=feat.type,
                qualifiers=feat.qualifiers,
            )
            sub_features.append(new_feat)

    sub_rec = SeqRecord(
        sub_seq,
        id=full_record.id,
        name=full_record.name,
        description=full_record.description,
    )
    # Keep the original record-level annotations (source organism,
    # taxonomy, references, etc.) as antiSMASH does for region records.
    sub_rec.annotations = dict(full_record.annotations)
    sub_rec.dbxrefs = list(full_record.dbxrefs)
    sub_rec.features = sub_features

    whole_record = (area_start == 0 and area_end == len(full_record.seq))
    comment_lines = [
        "##antiSMASH-Data-START##",
        f"Version     :: {antismash_version}",
        f"Run date    :: {run_date}",
    ]
    if not whole_record:
        comment_lines.append(
            "NOTE :: This is a single region extracted from a larger record!"
        )
        comment_lines.append(f"Orig. start :: {area_start}")
        comment_lines.append(f"Orig. end   :: {area_end}")
    else:
        comment_lines.append(
            "NOTE :: This is a single region extracted from a larger record!"
        )
        comment_lines.append(f"Orig. start :: {area_start}")
        comment_lines.append(f"Orig. end   :: {area_end}")
    comment_lines.append("##antiSMASH-Data-END##")
    sub_rec.annotations["comment"] = "\n".join(comment_lines)

    return sub_rec


def is_qualifying_nrps(area):
    """True if this area is a "true" NRPS BGC: its product list contains
    the exact rule name "NRPS" (not merely "NRPS-like"). Hybrid regions
    (e.g. products = ["NRPS", "T1PKS"]) still qualify, since "NRPS" is
    one of their products."""
    return "NRPS" in area.get("products", [])


def convert_one_file(json_path, out_dir, only_nrps=False):
    """Convert a single antiSMASH JSON file into per-region .gbk files
    written directly into out_dir. Returns the list of paths written."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(json_path) as fh:
        data = json.load(fh)

    antismash_version = data.get("version", "unknown")
    run_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    written = []
    for raw_record in data.get("records", []):
        full_record = build_full_record(raw_record)
        areas = raw_record.get("areas", [])
        # idx tracks the original antiSMASH region_number (1-based position
        # among ALL areas in this record), so filenames/numbering stay
        # consistent with antiSMASH's own output even when some regions
        # are skipped by a filter.
        for idx, area in enumerate(areas, start=1):
            products = area.get("products", [])
            if only_nrps and not is_qualifying_nrps(area):
                print(f"  Skipping {full_record.id} region{idx:03d} "
                      f"({','.join(products) or 'no product'}) - "
                      f"not a qualifying NRPS BGC")
                continue

            sub_rec = slice_region(
                full_record, area["start"], area["end"],
                antismash_version, run_date,
            )
            safe_id = full_record.id.replace(".", "_")
            out_name = f"{safe_id}_region{idx:03d}.gbk"
            out_path = out_dir / out_name
            with open(out_path, "w") as out_fh:
                SeqIO.write(sub_rec, out_fh, "genbank")
            written.append(out_path)
            print(f"  Wrote {out_path}  ({area['start']}-{area['end']}, "
                  f"{len(sub_rec.features)} features, product(s): "
                  f"{','.join(products)})")
    return written


def convert(input_dir, output_dir, only_nrps=False):
    """Iterate over every .json file in input_dir and convert each one.
    Each input file's regions are written into their own subfolder of
    output_dir (named after the JSON file, minus extension), so that
    identically-numbered regions from different antiSMASH runs never
    collide."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    json_files = sorted(input_dir.glob("*.json"))
    if not json_files:
        print(f"No .json files found in {input_dir}")
        return []

    all_written = []
    for json_path in json_files:
        print(f"Processing {json_path.name}")
        sub_out_dir = output_dir / json_path.stem
        written = convert_one_file(json_path, sub_out_dir, only_nrps=only_nrps)
        all_written.extend(written)

    print(f"\nDone: {len(all_written)} region file(s) written from "
          f"{len(json_files)} input file(s) to {output_dir}")
    return all_written


def build_parser():
    parser = argparse.ArgumentParser(
        description="Convert a folder of antiSMASH JSON results files "
                     "into per-region GenBank (.gbk) files."
    )
    parser.add_argument(
        "-i", "--input", required=True, dest="input_dir", type=Path,
        metavar="DIR",
        help="folder containing one or more antiSMASH .json result files",
    )
    parser.add_argument(
        "-o", "--output", required=True, dest="output_dir", type=Path,
        metavar="DIR",
        help="folder to write .gbk files to (a subfolder is created per "
             "input .json file); created if it doesn't already exist",
    )
    parser.add_argument(
        "--only-nrps", action="store_true",
        help="only write regions whose product list contains the exact "
             "rule name 'NRPS' (excludes NRPS-like-only regions; hybrid "
             "NRPS regions, e.g. NRPS/T1PKS, are still included)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = build_parser()
    convert(args.input_dir, args.output_dir, only_nrps=args.only_nrps)