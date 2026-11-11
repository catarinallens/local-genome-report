#!/usr/bin/env python3
"""
genome_lab.py — your private, local genome lab.

Everything runs on your own computer. Your DNA never leaves it; only public
research databases are downloaded.

ONE-TIME SET-UP (in Terminal, inside the conda environment):
    conda activate genome
    python3 genome_lab.py setup                 # reference genome + 1000 Genomes panel (~1-3 h)
    python3 genome_lab.py chip  MyHeritage_raw_dna_data.zip
    caffeinate -i python3 genome_lab.py impute MyHeritage_full_sequencing.cram   # several hours

EVERY TIME YOU COME BACK:
    python3 genome_lab.py run                   # refresh databases + papers, build the report

OTHER COMMANDS:
    python3 genome_lab.py lookup rs7903146      # everything known about one variant
    python3 genome_lab.py watch "Migraine" migraine
    python3 genome_lab.py status

All files live in ~/genome-lab (change with the GENOME_LAB_HOME variable).

Educational tool, not a medical test. Imputed genotypes are statistical
estimates; rare variants in particular need clinical confirmation.
"""

import argparse
import bisect
import csv
import datetime
import glob
import gzip
import hashlib
import html
import io
import json
import math
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time
import zipfile

VERSION = "1.0"
HOME = os.environ.get("GENOME_LAB_HOME", os.path.expanduser("~/genome-lab"))


def P(*parts):
    return os.path.join(HOME, *parts)


URLS = {
    "fasta": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz",
    "panel_dir": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
    "map": "https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/maps/genetic_maps.b37/chr{c}.b37.gmap.gz",
    "clinvar": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
    "gwas_dir": "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/",
    "europepmc": "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
    "chain": "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz",
}

CHROM_IDX = {str(i): i for i in range(1, 23)}
CHROM_IDX.update({"X": 23, "Y": 24, "MT": 25, "M": 25, "23": 23, "24": 24, "25": 25, "26": 25})


def gkey(chrom, pos):
    """Single integer key for a GRCh37 position."""
    c = CHROM_IDX.get(str(chrom).replace("chr", ""))
    if c is None or pos is None:
        return None
    return c * 1_000_000_000 + int(pos)

GENETIC_TERMS = ('GWAS OR "genome-wide" OR polymorphism OR polymorphisms OR "genetic variant" '
                 'OR "genetic variants" OR SNP OR SNPs OR "polygenic" OR "Mendelian randomization" '
                 'OR heritability OR "risk locus" OR "risk loci"')

DEFAULT_CONFIG = {
    "watch": {},   # conditions to follow closely; add with:  python3 genome_lab.py watch "Migraine" migraine
    "pmax": 5e-8,
    "chromosomes": [str(i) for i in range(1, 23)],
    "threads": 0,
    "papers_first_lookback_days": 365,
    "refresh_days": 7,
    "explorer_min_signals": 10,
    "sex": None,
    "ancestry": None,
    "urls": {},
}

# GRCh37 X outside the pseudo-autosomal regions; panel restricted to female samples
X_REGION = "X:2699521-154931043"
# Traits whose genetic effects differ substantially between women and men, so other-sex-only studies are left out.
SEX_DIFFERENT_TRAITS = re.compile(r"testosterone|estradiol|oestradiol|sex hormone|shbg|androgen|sexual dimorphism|"
                                  r"waist|hip circumference|fat distribution|body fat percentage|balding|hair loss|"
                                  r"facial hair|beard|voice|age at voice")
FEMALE_WORDS = re.compile(r"\b(female|females|women|woman|girls|mothers|maternal)\b")
MALE_WORDS = re.compile(r"\b(male|males|men|man|boys|fathers|paternal)\b")
MALE_ONLY_TRAITS = re.compile(r"prostat|testic|testis|erectile|male[- ]pattern|sperm|male infertility|"
                              r"hypospadias|cryptorchidism")
FEMALE_ONLY_TRAITS = re.compile(r"endometri|ovar|uter|cervical cancer|cervix|menarche|menopaus|menstrua|"
                                r"pregnan|preeclampsia|pre-eclampsia|polycystic ovary|female infertility|"
                                r"gestation|miscarriage|breastfeed")


EUROPEAN_WORDS = re.compile(r"\b(european|british|finnish|icelandic|dutch|swedish|danish|norwegian|german|"
                            r"italian|spanish|portuguese|french|irish|scottish|estonian|polish|greek|sardinian|"
                            r"croatian|hungarian|swiss|austrian|belgian|caucasian|white)\b")
OTHER_ANCESTRY_WORDS = re.compile(r"\b(african|asian|chinese|japanese|korean|han|hispanic|latino|latin|"
                                  r"indian|pakistani|bangladeshi|thai|vietnamese|filipino|taiwanese|malay|"
                                  r"native|hawaiian|arab|iranian|turkish|mexican|brazilian|peruvian|"
                                  r"greater middle eastern|oceanian|aboriginal)\b")


def study_ancestry(row):
    """'eur' (only European-ancestry samples), 'other' (none), 'mixed', or '' if not stated."""
    t = (row.get("INITIAL SAMPLE SIZE") or "").lower()
    e, o = bool(EUROPEAN_WORDS.search(t)), bool(OTHER_ANCESTRY_WORDS.search(t))
    return "eur" if e and not o else "other" if o and not e else "mixed" if e and o else ""


def study_sex(row):
    """'f' / 'm' if a catalog row comes from a one-sex study, else ''."""
    t = ((row.get("DISEASE/TRAIT") or "") + " | " + (row.get("INITIAL SAMPLE SIZE") or "")).lower()
    f, m = bool(FEMALE_WORDS.search(t)), bool(MALE_WORDS.search(t))
    return "f" if f and not m else "m" if m and not f else ""


def chromosomes(cfg):
    chs = list(cfg["chromosomes"])
    if cfg.get("sex") == "female" and "X" not in chs:
        chs.append("X")
    return chs

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
LD_WINDOW = 500_000
MIN_CONF = 0.9
STARS = [("practice guideline", 4), ("expert panel", 3),
         ("multiple submitters, no conflicts", 2), ("single submitter", 1), ("conflicting", 1)]


# =============================================================================
# Small helpers
# =============================================================================

def say(msg=""):
    print(msg, flush=True)


def die(msg):
    sys.exit("\n" + msg.strip() + "\n")


def today():
    return datetime.date.today().isoformat()


def load_config():
    path = P("config.json")
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    if os.path.exists(path):
        with open(path) as f:
            user = json.load(f)
        cfg.update(user)
    else:
        os.makedirs(HOME, exist_ok=True)
        save_config(cfg)
    if not cfg.get("threads"):
        cfg["threads"] = max(1, (os.cpu_count() or 2))
    return cfg


def save_config(cfg):
    os.makedirs(HOME, exist_ok=True)
    out = {k: v for k, v in cfg.items()}
    with open(P("config.json"), "w") as f:
        json.dump(out, f, indent=2)


def url(cfg, key):
    return cfg.get("urls", {}).get(key, URLS[key])


def need(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        die("Missing programs: " + ", ".join(missing) + "\n\n"
            "Make sure you ran:  conda activate genome\n"
            "If they're still missing, install them with:\n"
            "  conda install -y -n genome -c conda-forge -c bioconda "
            "\"samtools>=1.20\" \"bcftools>=1.20\" glimpse-bio=2.0.1")


def run(cmd, quiet=True):
    """Runs a command; on failure prints the tail of its error output."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE if quiet else None,
                       stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        tail = "\n".join((r.stderr or "").strip().splitlines()[-15:])
        die("This step failed:\n  " + " ".join(cmd) + "\n\n" + tail)
    return r.stdout if quiet else ""


def curl_text(u, params=None):
    cmd = ["curl", "-sSL", "--fail", "--retry", "3", "--max-time", "120"]
    if params:
        cmd += ["--get"]
        for k, v in params.items():
            cmd += ["--data-urlencode", "%s=%s" % (k, v)]
    cmd.append(u)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr.strip() or "download failed")
    return r.stdout


def download(u, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    say("  downloading " + u)
    r = subprocess.run(["curl", "-L", "--fail", "--retry", "5", "-C", "-",
                        "--progress-bar", "-o", tmp, u])
    if r.returncode not in (0, 33):  # 33 = server doesn't support resume, file complete
        raise RuntimeError("download failed: " + u)
    os.replace(tmp, dest)
    return dest


def intact(path):
    """True if a downloaded .gz/.zip file is complete and uncorrupted."""
    if path.endswith(".zip"):
        try:
            with zipfile.ZipFile(path) as z:
                return z.testzip() is None
        except zipfile.BadZipFile:
            return False
    if path.endswith(".gz"):
        # some published files carry harmless padding after the data ("trailing garbage ignored")
        return intact_gzip(path)
    return True


CHUNK = 16 * 1024 * 1024


def remote_size(u):
    """Size in bytes if the server supports partial downloads, else None."""
    if not u.startswith("http"):
        return None
    r = subprocess.run(["curl", "-sIL", "--max-time", "60", "--retry", "3", u], capture_output=True, text=True)
    size, ranges = None, False
    for line in r.stdout.splitlines():
        k, _, v = line.partition(":")
        if k.strip().lower() == "content-length" and v.strip().isdigit():
            size = int(v.strip())
        if k.strip().lower() == "accept-ranges" and "bytes" in v.lower():
            ranges = True
    return size if ranges and size and size > 2 * CHUNK else None


def fetch_piece(u, path, start, end, tries=8):
    """Downloads bytes start..end (inclusive) into a file of their own, retrying until complete."""
    want = end - start + 1
    for t in range(tries):
        r = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "600", "--retry", "3",
                            "-r", "%d-%d" % (start, end), "-o", path, u], capture_output=True)
        if r.returncode == 0 and os.path.exists(path) and os.path.getsize(path) == want:
            return
        time.sleep(min(60, 5 * (t + 1)))
    raise RuntimeError("couldn't download part of " + u)


def bad_bgzf_offsets(path):
    """Byte offsets of damaged blocks in a BGZF (bgzip) file; None if it isn't BGZF."""
    import struct
    import zlib
    bad = []
    with open(path, "rb") as f:
        head = f.read(18)
        if len(head) < 18 or head[:4] != b"\x1f\x8b\x08\x04" or head[12:14] != b"BC":
            return None
        f.seek(0)
        pos, total = 0, os.path.getsize(path)
        while pos < total:
            h = f.read(18)
            if len(h) < 18 or h[:4] != b"\x1f\x8b\x08\x04" or h[12:14] != b"BC":
                bad.append(pos)
                break  # can't find the next block reliably; re-fetch from here
            bsize = struct.unpack("<H", h[16:18])[0] + 1
            body = f.read(bsize - 18)
            if len(body) != bsize - 18:
                bad.append(pos)
                break
            data, crc, isize = body[:-8], body[-8:-4], body[-4:]
            try:
                out = zlib.decompress(data, -15)
                if struct.pack("<I", zlib.crc32(out) & 0xffffffff) != crc or struct.pack("<I", len(out)) != isize:
                    bad.append(pos)
            except zlib.error:
                bad.append(pos)
            pos += bsize
    return bad


def download_pieces(u, dest, size):
    """Downloads a big file in pieces; damaged pieces are fetched again instead of the whole file."""
    d = dest + ".pieces"
    os.makedirs(d, exist_ok=True)
    n = (size + CHUNK - 1) // CHUNK
    say("  downloading %s (%d MB, in %d pieces)" % (u, size // 1_000_000, n))
    for rnd in range(4):
        done = 0
        for i in range(n):
            pth = os.path.join(d, "%05d" % i)
            want = min(CHUNK, size - i * CHUNK)
            if not (os.path.exists(pth) and os.path.getsize(pth) == want):
                fetch_piece(u, pth, i * CHUNK, i * CHUNK + want - 1)
            done += 1
            if done % 8 == 0 or done == n:
                sys.stdout.write("\r  %d%%" % (100 * done // n))
                sys.stdout.flush()
        sys.stdout.write("\n")
        with open(dest + ".part", "wb") as out:
            for i in range(n):
                with open(os.path.join(d, "%05d" % i), "rb") as f:
                    shutil.copyfileobj(f, out, 1 << 22)
        bad = bad_bgzf_offsets(dest + ".part")
        if bad is None:
            ok = intact_gzip(dest + ".part")
            bad_pieces = set() if ok else set(range(n))
        else:
            bad_pieces = {o // CHUNK for o in bad} | {o // CHUNK + 1 for o in bad if o // CHUNK + 1 < n}
        if not bad_pieces:
            os.replace(dest + ".part", dest)
            shutil.rmtree(d)
            return dest
        say("  %d damaged piece(s); fetching them again…" % len(bad_pieces))
        for i in bad_pieces:
            os.remove(os.path.join(d, "%05d" % i))
        os.remove(dest + ".part")
    raise RuntimeError("download keeps arriving damaged: " + u)


def intact_gzip(path):
    r = subprocess.run(["gzip", "-t", path], capture_output=True, text=True)
    return r.returncode == 0 or (r.returncode == 2 and "trailing garbage" in (r.stderr or ""))


def download_checked(u, dest):
    """Download and verify integrity. Big files go in pieces so damage is repaired piece by piece."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        size = remote_size(u)
        if size:
            return download_pieces(u, dest, size)
        for attempt in (1, 2, 3):
            try:
                download(u, dest)
            except RuntimeError:
                time.sleep(10)
                continue
            if intact(dest):
                return dest
            say("  the download was damaged; trying again…")
            os.remove(dest)
    except RuntimeError:
        pass
    die("Download keeps failing for " + u + "\nThe connection may be unstable — try again on a steady wifi "
        "connection. Finished work is kept, so just run the same command again.")


def fmt_dur(sec):
    sec = int(sec)
    if sec < 90:
        return "%ds" % sec
    if sec < 5400:
        return "%dmin" % (sec // 60)
    return "%dh%02d" % (sec // 3600, (sec % 3600) // 60)


def to_float(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


# =============================================================================
# Database of your genotypes
# =============================================================================

def db():
    os.makedirs(HOME, exist_ok=True)
    con = sqlite3.connect(P("genome.db"))
    for table, col in (("seq", "k"), ("chip", "k")):  # older layouts get rebuilt
        cols = [r[1] for r in con.execute("PRAGMA table_info(%s)" % table)]
        if cols and col not in cols:
            con.execute("DROP TABLE %s" % table)
    # positions are GRCh37, encoded with gkey()
    con.execute("CREATE TABLE IF NOT EXISTS seq (k INTEGER PRIMARY KEY, ref TEXT, alt TEXT,"
                " gt TEXT, ds INTEGER, conf INTEGER)")
    con.execute("CREATE TABLE IF NOT EXISTS chip (rs INTEGER PRIMARY KEY, k INTEGER, gt TEXT)")
    con.execute("CREATE INDEX IF NOT EXISTS chip_k ON chip(k)")
    con.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT)")
    return con


def meta_get(con, k, default=None):
    row = con.execute("SELECT v FROM meta WHERE k=?", (k,)).fetchone()
    return row[0] if row else default


def meta_set(con, k, v):
    con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (k, str(v)))
    con.commit()


class Genome:
    """Looks up your genotypes: sequencing (imputed) first, chip as fallback."""

    def __init__(self):
        self.con = db()
        self.n_seq = self.con.execute("SELECT COUNT(*) FROM seq").fetchone()[0]
        self.n_chip = self.con.execute("SELECT COUNT(*) FROM chip").fetchone()[0]

    def get_many(self, keys):
        """keys: GRCh37 position keys -> {key: {"seq": (ref, alt, gt, ds, conf), "chip": "AG"}}"""
        out = {}
        ks = sorted({k for k in keys if k})
        for i in range(0, len(ks), 900):
            part = ks[i:i + 900]
            q = ",".join("?" * len(part))
            if self.n_seq:
                for k, ref, alt, gt, ds, conf in self.con.execute(
                        "SELECT * FROM seq WHERE k IN (%s)" % q, part):
                    out.setdefault(k, {})["seq"] = (ref, alt, gt, ds / 100.0, conf / 100.0)
            if self.n_chip:
                for k, gt in self.con.execute("SELECT k, gt FROM chip WHERE k IN (%s)" % q, part):
                    out.setdefault(k, {})["chip"] = gt
        return out

    def chip_keys(self, rsids):
        """rsID -> GRCh37 key, from the chip file's own positions."""
        out = {}
        nums = sorted({int(r[2:]) for r in rsids if r.startswith("rs") and r[2:].isdigit()})
        for i in range(0, len(nums), 900):
            part = nums[i:i + 900]
            for rs, k in self.con.execute(
                    "SELECT rs, k FROM chip WHERE rs IN (%s)" % ",".join("?" * len(part)), part):
                out["rs%d" % rs] = k
        return out


def seq_letters(ref, alt, gt):
    alle = [ref, alt]
    try:
        return "/".join(alle[int(x)] for x in gt.replace("|", "/").split("/"))
    except (ValueError, IndexError):
        return gt


def rec_text(rec):
    if not rec:
        return None
    if rec.get("seq"):
        ref, alt, gt, ds, conf = rec["seq"]
        return seq_letters(ref, alt, gt)
    if rec.get("chip"):
        return rec["chip"][0] + "/" + rec["chip"][1]
    return None


def allele_dose(rec, allele, flip=False):
    """Copies (0-2, can be fractional) of `allele` you carry.
    flip=True means the allele is given on the opposite strand (inverted region between builds).
    Returns (dose, genotype_text, source, confidence, note) or None."""
    if not rec or allele not in COMP:
        return None
    if flip:
        allele = COMP[allele]
    s = rec.get("seq")
    if s:
        ref, alt, gt, ds, conf = s
        a = allele
        if a not in (ref, alt) and COMP.get(a) in (ref, alt):
            a = COMP[a]
        dose = ds if a == alt else (2 - ds) if a == ref else None
        if dose is not None:
            note = "" if conf >= MIN_CONF else "low imputation confidence"
            return dose, seq_letters(ref, alt, gt), "seq", conf, note
    c = rec.get("chip")
    if c:
        note = ""
        if set(c) <= {allele, COMP[allele]}:
            if COMP[allele] in c and allele not in c and allele in "AT" and COMP[allele] in "AT" or \
                    (allele + COMP[allele]) in ("CG", "GC") and set(c) == {COMP[allele]}:
                note = "A/T or C/G variant: strand can't be verified"
        elif allele not in c and COMP[allele] in c:
            allele = COMP[allele]  # the study reported the other DNA strand
        return float(c.count(allele)), c[0] + "/" + c[1], "chip", None, note
    return None


class Liftover:
    """Minimal GRCh38 -> GRCh37 liftover from an Ensembl chain file."""

    def __init__(self, path):
        self.starts, self.blocks = {}, {}
        with gzip.open(path, "rt") as f:
            t_name = None
            for line in f:
                p = line.split()
                if not p:
                    continue
                if p[0] == "chain":
                    t_name, t_pos = p[2], int(p[5])
                    q_name, q_size, q_strand, q_pos = p[7], int(p[8]), p[9], int(p[10])
                    continue
                size = int(p[0])
                self.blocks.setdefault(t_name, []).append((t_pos, t_pos + size, q_name, q_pos, q_strand, q_size))
                if len(p) == 3:
                    t_pos += size + int(p[1])
                    q_pos += size + int(p[2])
        for c in self.blocks:
            self.blocks[c].sort()
            self.starts[c] = [b[0] for b in self.blocks[c]]

    def lift(self, chrom, pos):
        """1-based GRCh38 position -> (chrom, pos, flipped) on GRCh37, or None."""
        if not chrom or pos is None:
            return None
        chrom = str(chrom).replace("chr", "")
        st = self.starts.get(chrom)
        if not st:
            return None
        p0 = int(pos) - 1
        i = bisect.bisect_right(st, p0) - 1
        for j in range(i, max(i - 4, -1), -1):
            ts, te, qn, qs, strand, qsize = self.blocks[chrom][j]
            if ts <= p0 < te:
                q0 = qs + (p0 - ts)
                if strand == "-":
                    return qn, qsize - q0, True
                return qn, q0 + 1, False
        return None


def get_liftover(cfg):
    path = P("db", "GRCh38_to_GRCh37.chain.gz")
    if not os.path.exists(path):
        download_checked(url(cfg, "chain"), path)
    return Liftover(path)


# =============================================================================
# setup: reference genome, genetic maps, 1000 Genomes reference panel
# =============================================================================

def fasta_path():
    return P("ref", "human_g1k_v37.fasta")


def cmd_setup(a):
    cfg = load_config()
    need(["curl", "samtools", "bcftools", "GLIMPSE2_chunk", "GLIMPSE2_split_reference",
          "GLIMPSE2_phase", "GLIMPSE2_ligate"])
    free = shutil.disk_usage(os.path.dirname(HOME) if not os.path.exists(HOME) else HOME).free
    say("Genome lab folder: %s  (%.0f GB free)" % (HOME, free / 1e9))
    if free < 15e9:
        say("  Warning: under 15 GB free. The lab needs about 11 GB for good, 13 GB while imputing.")
    os.makedirs(P("ref"), exist_ok=True)

    fa = fasta_path()
    if not os.path.exists(fa + ".fai"):
        say("\n[1/2] Reference genome (GRCh37, ~900 MB download)")
        gz = fa + ".gz"
        if not os.path.exists(fa):
            if not os.path.exists(gz):
                download_checked(url(cfg, "fasta"), gz)
            say("  unpacking…")
            with open(fa + ".tmp", "wb") as dst:
                r = subprocess.run(["gzip", "-dc", gz], stdout=dst, stderr=subprocess.PIPE, text=True)
            if r.returncode not in (0, 2):
                die("Couldn't unpack the reference genome:\n" + r.stderr)
            os.replace(fa + ".tmp", fa)
            os.remove(gz)
        run(["samtools", "faidx", fa])
    say("\n[1/2] Reference genome: ready")

    say("\n[2/2] Genetic maps")
    for c in chromosomes(cfg):
        m = P("ref", "maps", "chr%s.b37.gmap.gz" % c)
        if not os.path.exists(m):
            download_checked(url(cfg, "map").format(c=c), m)
    say("  ready")
    say("\nSet-up complete. Next:  caffeinate -i python3 genome_lab.py impute YOUR_FILE.cram")


def panel_bcf(c, cfg, listing):
    """1000 Genomes panel for one chromosome, two-allele SNPs only (downloaded on demand)."""
    t = str(cfg["threads"])
    bcf = P("panel", "chr%s.bcf" % c)
    if os.path.exists(bcf + ".csi"):
        return bcf
    m = re.search(r'href="(ALL\.chr%s\.phase3[^"]*genotypes\.vcf\.gz)"' % re.escape(c), listing)
    fname = m.group(1) if m else ("ALL.chr%s.phase3_shapeit2_mvncall_integrated_%s.20130502.genotypes.vcf.gz"
                                  % (c, "v1c" if c == "X" else "v5b"))
    src = P("panel", "chr%s.src.vcf.gz" % c)
    if not os.path.exists(src):
        say("  chromosome %s: downloading the 1000 Genomes reference panel" % c)
        download_checked(url(cfg, "panel_dir") + fname, src)
    say("  chromosome %s: preparing the reference panel" % c)
    extra = []
    sheet = P("panel", "samples.panel")
    if not os.path.exists(sheet):
        download(url(cfg, "panel_dir") + "integrated_call_samples_v3.20130502.ALL.panel", sheet)
    if c == "X":
        # women's X chromosomes only (two copies, like the other chromosomes), outside the X/Y-shared ends
        females = [l.split()[0] for l in open(sheet) if len(l.split()) >= 4 and l.split()[3] == "female"]
        with open(P("panel", "females.txt"), "w") as f:
            f.write("\n".join(females) + "\n")
        extra = ["-S", P("panel", "females.txt"), "--force-samples", "-t", X_REGION]
    run(["bcftools", "view", "-m2", "-M2", "-v", "snps", "--threads", t] + extra
        + ["-Ob", "-o", bcf + ".tmp", src])
    os.replace(bcf + ".tmp", bcf)
    run(["bcftools", "index", "-f", bcf])
    os.remove(src)
    return bcf


def refpop_path(c):
    return P("refpop", "chr%s.bcf" % c)


def keep_refpop(c, bcf, t):
    """Keep a compact copy of the 1000 Genomes reference people (common SNPs only) for
    percentiles and the ancestry check. Roughly 3 GB for the whole genome."""
    out = refpop_path(c)
    if os.path.exists(out + ".csi"):
        return
    os.makedirs(P("refpop"), exist_ok=True)
    run(["bcftools", "view", "-q", "0.01:minor", "--threads", t, "-Ob", "-o", out + ".tmp", bcf])
    os.replace(out + ".tmp", out)
    run(["bcftools", "index", "-f", out])


# =============================================================================
# impute: your CRAM + reference panel -> your genotypes
# =============================================================================

def contig_md5(fa, contig):
    p = subprocess.Popen(["samtools", "faidx", fa, contig], stdout=subprocess.PIPE)
    h = hashlib.md5()
    next(p.stdout)  # header
    for line in p.stdout:
        h.update(line.strip().upper())
    p.wait()
    return h.hexdigest()


def cmd_impute(a):
    cfg = load_config()
    need(["samtools", "bcftools", "GLIMPSE2_chunk", "GLIMPSE2_split_reference", "GLIMPSE2_phase", "GLIMPSE2_ligate"])
    cram = os.path.abspath(os.path.expanduser(a.cram))
    if not os.path.exists(cram):
        die("Can't find " + cram)
    fa = fasta_path()
    if not os.path.exists(fa + ".fai") or not all(
            os.path.exists(P("ref", "maps", "chr%s.b37.gmap.gz" % c)) for c in chromosomes(cfg)):
        die("Run set-up first:  python3 genome_lab.py setup")
    t = str(cfg["threads"])
    os.makedirs(P("panel"), exist_ok=True)
    os.makedirs(P("work"), exist_ok=True)
    meta_set(db(), "cram_path", cram)

    if not (os.path.exists(cram + ".crai") or os.path.exists(cram[:-5] + ".crai")):
        say("Indexing your CRAM (one-off)…")
        run(["samtools", "index", cram])

    header = run(["samtools", "view", "-H", cram])
    first = cfg["chromosomes"][0]
    m = re.search(r"@SQ\tSN:%s\t.*?M5:([0-9a-f]{32})" % re.escape(first), header)
    if m:
        say("Checking that the reference genome matches your file…")
        if contig_md5(fa, first) != m.group(1):
            die("The reference genome doesn't match the one your CRAM was made with.\n"
                "Paste this message to Claude so it can adjust the set-up.")
        say("  match")

    chs = chromosomes(cfg)
    todo = [c for c in chs if not os.path.exists(P("work", "chr%s.bcf.ok" % c))]
    sizes = {c: os.path.getsize(P("ref", "maps", "chr%s.b37.gmap.gz" % c)) for c in chs}
    total_size = sum(sizes[c] for c in todo) or 1
    say("Imputing your genome: %d of %d chromosomes to go." % (len(todo), len(chs)))
    say("This takes many hours — best overnight with the Mac plugged in. You can stop any time "
        "with Ctrl+C; rerunning continues where it stopped.\n")
    listing = None
    t0, done_size = time.time(), 0
    for c in todo:
        if listing is None:
            try:
                listing = curl_text(url(cfg, "panel_dir"))
            except Exception:
                listing = ""
        bcf = panel_bcf(c, cfg, listing)
        keep_refpop(c, bcf, t)
        mp = P("ref", "maps", "chr%s.b37.gmap.gz" % c)
        chunks = P("panel", "chunks.chr%s.txt" % c)
        if not os.path.exists(chunks):
            run(["GLIMPSE2_chunk", "--input", bcf, "--region", X_REGION if c == "X" else c, "--map", mp, "--sequential",
                 "--output", chunks + ".tmp", "--threads", t])
            os.replace(chunks + ".tmp", chunks)
        lines = [l.split() for l in open(chunks) if l.strip()]
        d = P("work", "chr%s" % c)
        os.makedirs(d, exist_ok=True)
        outs = []
        for j, cols in enumerate(lines):
            out = os.path.join(d, "%03d.bcf" % j)
            outs.append(out)
            if os.path.exists(out + ".ok"):
                continue
            prefix = os.path.join(d, "ref%03d" % j)
            for old in glob.glob(prefix + "_*.bin"):
                os.remove(old)
            run(["GLIMPSE2_split_reference", "--reference", bcf, "--map", mp,
                 "--input-region", cols[2], "--output-region", cols[3], "--output", prefix, "--threads", t])
            binf = glob.glob(prefix + "_*.bin")
            if len(binf) != 1:
                die("Couldn't prepare reference chunk %d of chromosome %s." % (j + 1, c))
            run(["GLIMPSE2_phase", "--bam-file", cram, "--reference", binf[0],
                 "--fasta", fa, "--output", out, "--threads", t, "--ind-name", "me"])
            os.remove(binf[0])
            open(out + ".ok", "w").close()
            say("  chromosome %-2s  part %d of %d" % (c, j + 1, len(lines)))
        lst = os.path.join(d, "list.txt")
        with open(lst, "w") as f:
            f.write("\n".join(outs) + "\n")
        lig = P("work", "chr%s.bcf" % c)
        run(["GLIMPSE2_ligate", "--input", lst, "--output", lig, "--threads", t])
        if not os.path.exists(lig + ".csi"):
            run(["bcftools", "index", "-f", lig])
        open(lig + ".ok", "w").close()
        shutil.rmtree(d)
        for f in (bcf, bcf + ".csi", chunks):
            if os.path.exists(f):
                os.remove(f)
        done_size += sizes[c]
        el = time.time() - t0
        left = el / done_size * (total_size - done_size)
        say("  chromosome %s finished — about %s to go\n" % (c, fmt_dur(left)) if done_size < total_size
            else "  chromosome %s finished\n" % c)

    say("Combining chromosomes…")
    final = P("my_genome.imputed.bcf")
    # per-chromosome results are kept (small), so chromosomes can be added later without redoing the rest
    run(["bcftools", "concat", "--threads", t, "-Ob", "-o", final]
        + [P("work", "chr%s.bcf" % c) for c in chs])
    run(["bcftools", "index", "-f", final])
    build_seq_db(final)
    say("\nDone. Your imputed genome is in %s\nNext:  python3 genome_lab.py run" % final)


def build_seq_db(bcf):
    say("Building your genotype database (about 10-30 min)…")
    con = db()
    con.execute("DROP TABLE IF EXISTS seq")
    con.execute("CREATE TABLE seq (k INTEGER PRIMARY KEY, ref TEXT, alt TEXT, gt TEXT, ds INTEGER, conf INTEGER)")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA journal_mode=OFF")
    hdr = run(["bcftools", "view", "-h", bcf])
    af = "%INFO/RAF" if "##INFO=<ID=RAF," in hdr else "."
    p = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t" + af + "[\t%GT\t%DS\t%GP]\n", bcf],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    batch, n, skipped = [], 0, 0
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        try:
            ds = int(round(float(f[6]) * 100))
            conf = int(round(max(float(x) for x in f[7].split(",")) * 100))
            raf = float(f[4]) if f[4] not in (".", "") else 0.5
        except ValueError:
            continue
        # Very rare panel variants you almost certainly don't carry add size, not information.
        if raf < 0.001 or raf > 0.999:
            has_minor = ds >= 50 if raf < 0.5 else ds <= 150
            if not has_minor:
                skipped += 1
                continue
        k = gkey(f[0], f[1])
        if k is None:
            continue
        batch.append((k, f[2], f[3], f[5].replace("|", "/"), ds, conf))
        if len(batch) >= 200000:
            con.executemany("INSERT OR REPLACE INTO seq VALUES (?,?,?,?,?,?)", batch)
            n += len(batch)
            batch = []
            if n % 5000000 < 200000:
                say("  %d million variants…" % (n // 1000000))
    if batch:
        con.executemany("INSERT OR REPLACE INTO seq VALUES (?,?,?,?,?,?)", batch)
        n += len(batch)
    p.wait()
    con.commit()
    meta_set(con, "seq_date", today())
    meta_set(con, "seq_count", n)
    say("  %s variants stored" % format(n, ","))


# =============================================================================
# chip: MyHeritage / 23andMe / Ancestry raw data
# =============================================================================

def chip_header_method(path):
    """The '##method=' line some raw-data files carry (MyHeritage writes 'Low-pass Whole Genome Sequencing')."""
    try:
        if path.endswith(".zip"):
            z = zipfile.ZipFile(path)
            member = next(n for n in z.namelist() if n.lower().endswith((".csv", ".txt")))
            f = io.TextIOWrapper(z.open(member), encoding="utf-8", errors="replace")
        else:
            f = open(path, encoding="utf-8", errors="replace")
        with f:
            for i, line in enumerate(f):
                if line.lower().startswith("##method="):
                    return line.split("=", 1)[1].strip()
                if i > 40:
                    break
    except Exception:
        pass
    return ""


_CHIP_KIND = [None]


def chip_is_lowpass():
    """True when the 'chip' file was itself derived from low-pass sequencing (not an independent chip)."""
    if _CHIP_KIND[0] is None:
        con = db()
        m = meta_get(con, "chip_method")
        if m is None:
            name = meta_get(con, "chip_file", "")
            m = ""
            for d in (os.getcwd(), os.path.expanduser("~/Downloads")):
                if name and os.path.exists(os.path.join(d, name)):
                    m = chip_header_method(os.path.join(d, name))
                    break
            if m:
                meta_set(con, "chip_method", m)
        _CHIP_KIND[0] = "low-pass" in (m or "").lower() or "sequencing" in (m or "").lower()
    return _CHIP_KIND[0]


def chip_word():
    return "MyHeritage file" if chip_is_lowpass() else "chip"


def read_chip_file(path):
    gts = {}
    if path.endswith(".zip"):
        z = zipfile.ZipFile(path)
        member = next(n for n in z.namelist() if n.lower().endswith((".csv", ".txt")))
        f = io.TextIOWrapper(z.open(member), encoding="utf-8", errors="replace")
    elif path.endswith(".gz"):
        f = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        f = open(path, encoding="utf-8", errors="replace")
    with f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = next(csv.reader([line.strip()])) if "," in line else line.split()
            parts = [x.strip().strip('"') for x in parts]
            if len(parts) < 4 or not parts[0].lower().startswith("rs"):
                continue
            gt = (parts[3] + parts[4]) if len(parts) == 5 else parts[3]
            gt = gt.upper()
            rid = parts[0].lower()
            if len(gt) == 2 and all(ch in "ACGT" for ch in gt) and rid[2:].isdigit():
                try:
                    k = gkey(parts[1], int(parts[2]))
                except ValueError:
                    k = None
                gts[int(rid[2:])] = (k, gt)
    return gts


def cmd_chip(a):
    path = os.path.expanduser(a.file)
    if not os.path.exists(path):
        die("Can't find " + path)
    gts = read_chip_file(path)
    if len(gts) < 1000:
        die("Only %d genotypes found — is this the raw data file?" % len(gts))
    con = db()
    con.execute("DELETE FROM chip")
    con.executemany("INSERT OR REPLACE INTO chip VALUES (?,?,?)", ((rs, k, gt) for rs, (k, gt) in gts.items()))
    con.commit()
    meta_set(con, "chip_date", today())
    meta_set(con, "chip_file", os.path.basename(path))
    meta_set(con, "chip_method", chip_header_method(path))
    say("Stored %s chip genotypes. Next:  python3 genome_lab.py run" % format(len(gts), ","))


# =============================================================================
# update: research databases + papers
# =============================================================================

def file_age_days(path):
    if not os.path.exists(path):
        return 1e9
    return (time.time() - os.path.getmtime(path)) / 86400


def cmd_update(a, cfg=None):
    cfg = cfg or load_config()
    os.makedirs(P("db"), exist_ok=True)
    days = 0 if getattr(a, "force", False) else cfg["refresh_days"]

    gw = P("db", "gwas.zip")
    if file_age_days(gw) > days:
        say("Updating GWAS Catalog…")
        names = re.findall(r'href="([^"/]*associations_ontology-annotated[^"/]*\.zip)"',
                           curl_text(url(cfg, "gwas_dir")))
        names = [n for n in names if "split" not in n] or ["gwas-catalog-associations_ontology-annotated-full.zip"]
        download_checked(url(cfg, "gwas_dir") + names[0], gw)
    else:
        say("GWAS Catalog is recent (%d days old)" % file_age_days(gw))

    cv = P("db", "clinvar.txt.gz")
    if file_age_days(cv) > days:
        say("Updating ClinVar…")
        download_checked(url(cfg, "clinvar"), cv)
    else:
        say("ClinVar is recent (%d days old)" % file_age_days(cv))

    update_papers(cfg)


def update_papers(cfg):
    path = P("db", "papers.json")
    store = json.load(open(path)) if os.path.exists(path) else {"papers": {}, "fetched": {}}
    for name, w in cfg["watch"].items():
        q = w.get("papers")
        if not q:
            continue
        since = store["fetched"].get(name)
        if not since:
            since = (datetime.date.today() - datetime.timedelta(days=cfg["papers_first_lookback_days"])).isoformat()
        query = "TITLE_ABS:(%s) AND TITLE_ABS:(%s) AND FIRST_PDATE:[%s TO 3000-12-31]" % (q, GENETIC_TERMS, since)
        cursor, n_new, total = "*", 0, 0
        try:
            for _ in range(10):
                d = json.loads(curl_text(url(cfg, "europepmc"), {
                    "query": query, "format": "json", "resultType": "core",
                    "pageSize": "100", "cursorMark": cursor}))
                res = d.get("resultList", {}).get("result", [])
                for r in res:
                    pid = "%s:%s" % (r.get("source"), r.get("id"))
                    text = (r.get("title") or "") + " " + (r.get("abstractText") or "")
                    rsids = sorted({x.lower() for x in re.findall(r"\brs\d{3,}\b", text, re.I)})
                    if pid not in store["papers"]:
                        n_new += 1
                        store["papers"][pid] = {
                            "title": re.sub(r"<[^>]+>", "", r.get("title") or ""),
                            "authors": r.get("authorString", ""),
                            "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title")
                                       or r.get("bookOrReportDetails", {}).get("publisher", "")
                                       or ("Preprint" if r.get("source") == "PPR" else ""),
                            "date": r.get("firstPublicationDate", ""),
                            "url": "https://europepmc.org/article/%s/%s" % (r.get("source"), r.get("id")),
                            "rsids": rsids, "topics": [name], "first_seen": today(),
                            "preprint": r.get("source") == "PPR",
                        }
                    elif name not in store["papers"][pid]["topics"]:
                        store["papers"][pid]["topics"].append(name)
                total += len(res)
                nxt = d.get("nextCursorMark")
                if not res or not nxt or nxt == cursor:
                    break
                cursor = nxt
            store["fetched"][name] = today()
            say("Papers — %s: %d new" % (name, n_new))
        except Exception as e:
            say("Papers — %s: couldn't reach Europe PMC (%s)" % (name, e))
    with open(path, "w") as f:
        json.dump(store, f)


# =============================================================================
# Analysis
# =============================================================================

def open_gwas():
    path = P("db", "gwas.zip")
    if not os.path.exists(path):
        die("Research databases missing. Run:  python3 genome_lab.py run")
    z = zipfile.ZipFile(path)
    member = max(z.namelist(), key=lambda n: z.getinfo(n).file_size)
    return io.TextIOWrapper(z.open(member), encoding="utf-8", errors="replace")


def odds_ratio(row):
    """The study's odds ratio for the risk allele, or None if it reported something else (a beta)."""
    v = to_float(row.get("OR or BETA"))
    ci = (row.get("95% CI (TEXT)") or "").lower()
    if v is None or v <= 0 or any(w in ci for w in ("increase", "decrease", "unit", " sd", "z-score")):
        return None
    return v


RS_RISK = re.compile(r"^(rs\d+)-([ACGT])$", re.I)


def read_gwas(cfg, lift, want_rs=()):
    """One pass over the catalog.
    Returns (watched rows per trait, explorer rows, total rows, rsID -> GRCh37 (key, flip) for want_rs)."""
    pmax = cfg["pmax"]
    sex = {"female": "f", "male": "m"}.get(cfg.get("sex") or "", "")
    other_only = MALE_ONLY_TRAITS if sex == "f" else FEMALE_ONLY_TRAITS if sex == "m" else None
    european = cfg.get("ancestry") == "european"
    want_rs = set(want_rs)
    rs_pos = {}
    watch = {name: [k.lower() for k in w.get("gwas", [])] for name, w in cfg["watch"].items()}
    watched = {name: [] for name in watch}
    excluded, excluded_where = {}, {}
    explorer = []
    total = 0
    lift_cache = {}
    with open_gwas() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            chrom = (row.get("CHR_ID") or "").strip()
            pos = to_float(row.get("CHR_POS"))
            if want_rs and chrom and pos and ";" not in chrom:
                rid = (row.get("SNPS") or "").strip().lower()
                if rid in want_rs and rid not in rs_pos:
                    l = lift.lift(chrom, pos)
                    if l:
                        rs_pos[rid] = (gkey(l[0], l[1]), l[2])
            p = to_float(row.get("P-VALUE"))
            if p is None or p > pmax:
                continue
            m = RS_RISK.match((row.get("STRONGEST SNP-RISK ALLELE") or "").replace(" ", ""))
            if not m or not chrom or pos is None:
                continue
            rsid, risk = m.group(1).lower(), m.group(2).upper()
            ssex = study_sex(row)
            other_sex = bool(sex and ssex and ssex != sex)
            if other_sex and SEX_DIFFERENT_TRAITS.search(((row.get("DISEASE/TRAIT") or "") + " " +
                                                          (row.get("MAPPED_TRAIT") or "")).lower()):
                continue  # a trait whose genetics genuinely differ between the sexes; a study in the other sex won't transfer
            anc = study_ancestry(row)
            other_pop_ok = european and anc == "other" and refpop_ready()
            if european and anc == "other" and not other_pop_ok:
                t_ = ((row.get("DISEASE/TRAIT") or "") + " | " + (row.get("MAPPED_TRAIT") or "")).lower()
                for name, kws in watch.items():
                    if any(kw in t_ for kw in kws):
                        excluded[name] = excluded.get(name, 0) + 1
                        excluded_where.setdefault(name, {})
                        for w_ in re.findall(OTHER_ANCESTRY_WORDS, (row.get("INITIAL SAMPLE SIZE") or "").lower()):
                            excluded_where[name][w_] = excluded_where[name].get(w_, 0) + 1
                continue  # frequencies from a different population would skew the comparison
            ck = (chrom, pos)
            if ck not in lift_cache:
                l = lift.lift(chrom, pos)
                lift_cache[ck] = (gkey(l[0], l[1]), l[2]) if l else (None, False)
            k, flip = lift_cache[ck]
            raf = to_float(row.get("RISK ALLELE FREQUENCY"))
            orv = odds_ratio(row)
            mapped = (row.get("MAPPED_TRAIT") or "").strip()
            if mapped and not other_pop_ok and not (other_only and other_only.search(mapped.lower())):
                explorer.append((mapped, rsid, risk, p, raf, chrom, pos, k, flip, orv,
                                 (2 if other_sex else 0 if ssex == sex and sex else 1, european and anc == "eur"),
                                 (row.get("95% CI (TEXT)") or "").strip()))
            text = ((row.get("DISEASE/TRAIT") or "") + " | " + mapped).lower()
            for name, kws in watch.items():
                if any(kw in text for kw in kws):
                    watched[name].append({
                        "rsid": rsid, "risk": risk, "p": p, "raf": raf, "chr": chrom, "pos": pos,
                        "k": k, "flip": flip, "or": orv, "own_sex": ssex == sex and bool(sex), "other_sex": other_sex,
                        "own_anc": european and anc == "eur", "other_pop": other_pop_ok,
                        "effect": (row.get("OR or BETA") or "").strip(),
                        "ci": (row.get("95% CI (TEXT)") or "").strip(),
                        "trait": (row.get("DISEASE/TRAIT") or "").strip(),
                        "gene": (row.get("MAPPED_GENE") or "").strip(),
                        "author": (row.get("FIRST AUTHOR") or "").strip(),
                        "year": (row.get("DATE") or "")[:4],
                        "pmid": (row.get("PUBMEDID") or "").strip(),
                        "acc": (row.get("STUDY ACCESSION") or "").strip(),
                    })
    cfg["_excluded"] = {n: (c, max(excluded_where.get(n, {"": 0}).items(), key=lambda x: x[1])[0])
                        for n, c in excluded.items()}
    return watched, explorer, total, rs_pos


def eur_freqs(keys):
    """GRCh37 key -> (ref, alt, European ALT frequency) from the kept 1000 Genomes reference people."""
    out = {}
    by_chr = {}
    for k in keys:
        if k:
            by_chr.setdefault(k // 1_000_000_000, []).append(k % 1_000_000_000)
    names = {v: k for k, v in CHROM_IDX.items() if not k.isdigit() or int(k) <= 22}
    names.update({i: str(i) for i in range(1, 23)})
    for ci, poss in by_chr.items():
        c = names.get(ci)
        path = refpop_path(c) if c else None
        if not path or not os.path.exists(path + ".csi"):
            continue
        reg = P("db", "freq_regions.txt")
        with open(reg, "w") as f:
            f.write("".join("%s\t%d\n" % (c, p) for p in sorted(set(poss))))
        txt = subprocess.run(["bcftools", "query", "-R", reg, "-f", "%POS\t%REF\t%ALT\t%INFO/EUR_AF\n", path],
                             capture_output=True, text=True).stdout
        for line in txt.splitlines():
            p = line.split("\t")
            af = to_float(p[3].split(",")[0]) if len(p) > 3 else None
            if af is not None:
                out[gkey(c, int(p[0]))] = (p[1], p[2], af)
    return out


def use_eur_freqs(records, freqs):
    """Replaces each study's reported frequency with the European one for the same risk allele."""
    for r in records:
        f = freqs.get(r.get("k"))
        if not f:
            if r.get("other_pop"):
                r["raf"] = None  # a frequency from a different population would skew the comparison
            continue
        ref, alt, af = f
        a = COMP[r["risk"]] if r.get("flip") else r["risk"]
        if a not in (ref, alt) and COMP.get(a) in (ref, alt):
            a = COMP[a]
        if a == alt:
            r["raf"], r["raf_eur"] = af, True
        elif a == ref:
            r["raf"], r["raf_eur"] = 1 - af, True


CI_RE = re.compile(r"(\d+(?:\.\d+)?)\s*[-–,]\s*(\d+(?:\.\d+)?)")
OR_PLAUSIBLE = (0.1, 10.0)   # outside this, a common variant's odds ratio is almost always a data-entry error
TAU_PRIOR = 0.10             # typical spread of true log odds ratios for common variants (about ×1.1)


def z_from_p(p):
    """Two-sided p-value -> |z|, also for extremely small p-values."""
    p = max(p, 1e-300)
    if p > 1e-15:
        lo, hi = 0.0, 9.0   # solve erfc(z / sqrt 2) = p (works on older Pythons too)
        for _ in range(60):
            mid = (lo + hi) / 2
            if math.erfc(mid / math.sqrt(2)) > p:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    L = -math.log(p / 2)
    z = math.sqrt(2 * L)
    for _ in range(20):
        z = math.sqrt(max(2 * (L - math.log(z) - 0.5 * math.log(2 * math.pi)), 1e-6))
    return z


def log_or_se(r):
    """Standard error of the study's log odds ratio: from its confidence interval, else from its p-value."""
    b = math.log(r["or"])
    m = CI_RE.search(r.get("ci") or "")
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        if 0 < lo < hi and lo <= r["or"] * 1.02 and hi >= r["or"] * 0.98:
            return (math.log(hi) - math.log(lo)) / (2 * 1.96)
    if b and r.get("p"):
        return abs(b) / z_from_p(r["p"])
    return None


def shrink_effects(counted):
    """Empirical Bayes: pulls each log odds ratio toward zero in proportion to how uncertain it is.
    A precise estimate from a huge study keeps nearly its full size; a noisy one from a small study,
    which is also where the winner's curse hits hardest, is pulled in strongly."""
    pts = []
    for r in counted:
        r["or_used"], r["or_note"] = None, ""
        if not r.get("or"):
            continue
        if not (OR_PLAUSIBLE[0] <= r["or"] <= OR_PLAUSIBLE[1]):
            r["or_note"] = "effect size implausible (likely a catalog error), not used"
            continue
        se = log_or_se(r)
        if not se:
            continue
        pts.append((r, math.log(r["or"]), se))
    if not pts:
        return
    # spread of true effects for this trait (method of moments), pooled with a typical value
    est = sum(b * b - se * se for _, b, se in pts) / len(pts)
    k = 10
    tau2 = max((len(pts) * max(est, 0.0) + k * TAU_PRIOR ** 2) / (len(pts) + k), 0.03 ** 2)
    for r, b, se in pts:
        r["or_used"] = math.exp(b * tau2 / (tau2 + se * se))
        r["se"] = se

def tally(records, genos):
    """records: dicts with rsid/risk/p/raf/chr/pos/k/flip. Independent-signal risk-allele tally."""
    best = {}
    # for the same variant: a study in your ancestry wins, then one in your sex, then the strongest result
    # (a study in your sex, then a mixed one, then one in the other sex)
    rank = lambda r: (0 if r.get("own_anc") else 1, 0 if r.get("own_sex") else 2 if r.get("other_sex") else 1, r["p"])
    for r in records:
        prev = best.get(r["rsid"])
        if prev is None or rank(r) < rank(prev):
            r = dict(r, studies=(prev["studies"] + 1) if prev else 1)
            best[r["rsid"]] = r
        else:
            prev["studies"] += 1
    recs = sorted(best.values(), key=lambda r: r["p"])
    kept = []
    kept_pos = {}  # chrom -> sorted positions of kept signals

    def near_kept(chrom, pos):
        if pos is None:
            return False
        lst = kept_pos.get(chrom)
        if not lst:
            return False
        i = bisect.bisect_left(lst, pos)
        return (i < len(lst) and lst[i] - pos < LD_WINDOW) or (i > 0 and pos - lst[i - 1] < LD_WINDOW)

    for r in recs:
        res = allele_dose(genos.get(r["k"]), r["risk"], r["flip"])
        if res is None:
            r["status"] = "no data for you"
            r["dose"] = None
            continue
        r["dose"], r["gt"], r["src"], r["conf"], r["note"] = res
        if near_kept(r["chr"], r["pos"]):
            r["status"] = "same region as a stronger hit"
            continue
        if r["raf"] is None or not (0 < r["raf"] < 1):
            r["status"] = "no population frequency reported"
        else:
            r["status"] = "counted"
        kept.append(r)
        if r["pos"] is not None:
            bisect.insort(kept_pos.setdefault(r["chr"], []), r["pos"])
    counted = [r for r in kept if r["status"] == "counted"]
    obs = sum(r["dose"] for r in counted)
    exp = sum(2 * r["raf"] for r in counted)
    # Robust score: each signal contributes a capped standardized deviation, so one
    # misreported frequency in the catalog can't dominate the total.
    contrib = []
    for r in counted:
        f = min(max(r["raf"], 0.02), 0.98)
        contrib.append(max(-3.0, min(3.0, (r["dose"] - 2 * f) / math.sqrt(2 * f * (1 - f)))))
    # Relative odds vs a typical person: each signal multiplies the odds by OR^(your copies - typical copies),
    # using each odds ratio after shrinkage by its uncertainty (see shrink_effects).
    shrink_effects(counted)
    log_odds, n_or = 0.0, 0
    for r in counted:
        if r.get("or_used"):
            log_odds += (r["dose"] - 2 * r["raf"]) * math.log(r["or_used"])
            n_or += 1
    z = sum(contrib) / math.sqrt(len(contrib)) if contrib else None
    odds = math.exp(log_odds) if n_or >= 5 else None
    capped = sum(1 for c in contrib if abs(c) >= 3.0)
    suspicious = (z is not None and abs(z) > 6) or (odds is not None and not (1 / 5 < odds < 5)) \
        or (len(contrib) >= 5 and capped >= 0.6 * len(contrib))
    return {
        "records": recs, "n_total": len(recs),
        "n_data": sum(1 for r in recs if r["dose"] is not None),
        "n_counted": len(counted), "observed": obs, "expected": exp,
        "z": z, "odds": odds, "n_or": n_or, "suspicious": suspicious,
        "n_own_sex": sum(1 for r in counted if r.get("own_sex")),
        "n_other_pop": sum(1 for r in counted if r.get("other_pop")),
    }


MEASURE_WORDS = ("measurement", "level", "count", "ratio", "amount", "concentration", "volume",
                 "density", "percentage", "rate", "size", "thickness", "area", "width", "height", "mass",
                 "body shape", "age at", "quantity", "intake", "time")
APPEARANCE_RE = re.compile(r"hair colou?r|eye colou?r|skin colou?r|skin pigment|pigmentation|freckl|"
                           r"suntan|tanning|sunburn|hair shape|curly hair|hair morpholog|balding|"
                           r"facial|nose|chin|eyebrow|lip shape|ear lobe|earlobe|unibrow")
MEASURE_RE = re.compile(r"\b(%s)s?\b" % "|".join(re.escape(w) for w in MEASURE_WORDS))


def explore_traits(explorer, genos, min_signals):
    """Returns (conditions, measurements, suspicious) lists of dicts."""
    by_trait = {}
    for mapped, rsid, risk, p, raf, chrom, pos, k, flip, orv, (sex_rank, own_anc), ci in explorer:
        own_sex, other_sex = sex_rank == 0, sex_rank == 2
        by_trait.setdefault(mapped, []).append({"rsid": rsid, "risk": risk, "p": p, "raf": raf, "chr": chrom,
                                                "pos": pos, "k": k, "flip": flip, "or": orv, "ci": ci,
                                                "other_sex": other_sex,
                                                "own_sex": own_sex, "own_anc": own_anc})
    conditions, measurements, suspicious, appearance = [], [], [], []
    for trait, rows in by_trait.items():
        t = tally(rows, genos)
        if t["n_counted"] < min_signals or t["z"] is None:
            continue
        item = {"trait": trait, "n": t["n_counted"], "obs": t["observed"], "exp": t["expected"],
                "z": t["z"], "odds": t["odds"]}
        share_or = sum(1 for r in rows if r["or"]) / len(rows)
        is_measure = bool(MEASURE_RE.search(trait.lower())) or share_or < 0.5
        if t["suspicious"]:
            suspicious.append(item)
        elif APPEARANCE_RE.search(trait.lower()):
            appearance.append(item)
        elif is_measure or t["odds"] is None:
            measurements.append(item)
        else:
            conditions.append(item)
    conditions.sort(key=lambda x: -x["z"])
    measurements.sort(key=lambda x: -abs(x["z"]))
    suspicious.sort(key=lambda x: -abs(x["z"]))
    appearance.sort(key=lambda x: -abs(x["z"]))
    return conditions, measurements, suspicious, appearance


def stars_for(review):
    review = review.lower()
    for key, s in STARS:
        if key in review:
            return s
    return 0


def clinvar_rows(want_rs=()):
    """Germline pathogenic/likely pathogenic GRCh37 SNVs, plus rsID -> key for want_rs."""
    path = P("db", "clinvar.txt.gz")
    rows, rs_pos = [], {}
    if not os.path.exists(path):
        return rows, rs_pos
    want_rs = set(want_rs)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}

        def col(*names):
            for n in names:
                if n in idx:
                    return idx[n]
            raise KeyError("ClinVar column missing: %s" % (names,))
        c = {k: col(*v) for k, v in {
            "sig": ("ClinicalSignificance", "GermlineClassification"), "rs": ("RS# (dbSNP)",),
            "asm": ("Assembly",), "chrom": ("Chromosome",), "pos": ("PositionVCF",),
            "alt": ("AlternateAlleleVCF",), "ref": ("ReferenceAlleleVCF",),
            "gene": ("GeneSymbol",), "name": ("Name",), "pheno": ("PhenotypeList",),
            "review": ("ReviewStatus",), "origin": ("OriginSimple",), "vid": ("VariationID",)}.items()}
        mx = max(c.values())
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= mx or p[c["asm"]] != "GRCh37":
                continue
            try:
                k = gkey(p[c["chrom"]], int(p[c["pos"]]))
            except ValueError:
                continue
            rs = "rs" + p[c["rs"]] if p[c["rs"]] not in ("-1", "", "-") else ""
            if rs in want_rs:
                rs_pos.setdefault(rs, (k, False))
            sl = p[c["sig"]].lower()
            if "pathogenic" not in sl or "conflicting" in sl or "benign" in sl:
                continue
            if "germline" not in p[c["origin"]].lower():
                continue
            alt, ref = p[c["alt"]].upper(), p[c["ref"]].upper()
            if len(alt) != 1 or len(ref) != 1 or alt not in "ACGT" or ref not in "ACGT":
                continue
            if alt == ref or p[c["name"]].rstrip().endswith("="):
                continue  # entry describes the standard reference version, not a change  # single-letter changes only; insertions/deletions can't be read from this data
            rows.append({"k": k, "rsid": rs, "alt": alt, "ref": ref,
                         "sig": p[c["sig"]], "gene": p[c["gene"]], "name": p[c["name"]],
                         "pheno": p[c["pheno"]].replace("|", "; "),
                         "stars": stars_for(p[c["review"]]), "vid": p[c["vid"]]})
    return rows, rs_pos


def clinvar_hits(rows, genome):
    genos = genome.get_many([r["k"] for r in rows])
    hits = []
    for r in rows:
        g = genos.get(r["k"])
        if not g:
            continue
        s = g.get("seq")
        if s and (s[0], s[1]) != (r["ref"], r["alt"]):
            g = {"chip": g["chip"]} if g.get("chip") else None  # a different change at the same position
        if g and g.get("chip") and not set(g["chip"]) <= {r["ref"], r["alt"]}:
            g = {"seq": g["seq"]} if g.get("seq") else None
        if not g:
            continue
        res = allele_dose(g, r["alt"])
        if res is None:
            continue
        dose, gt, src, conf, note = res
        if dose < 0.5 or (src == "seq" and conf < MIN_CONF):
            continue
        r = dict(r, dose=dose, gt=gt, src=src, conf=conf)
        hits.append(r)
    seen, uniq = set(), []
    for r in sorted(hits, key=lambda r: (-r["stars"], r["gene"])):
        if (r["k"], r["alt"]) not in seen:
            seen.add((r["k"], r["alt"]))
            uniq.append(r)
    return uniq


def chip_vs_seq_qc(genome):
    if not (genome.n_seq and genome.n_chip):
        return None
    agree = total = nonref_agree = nonref_total = 0
    for chip, ref, alt, gt, conf in genome.con.execute(
            "SELECT c.gt, s.ref, s.alt, s.gt, s.conf FROM chip c JOIN seq s ON c.k = s.k"):
        if conf < MIN_CONF * 100 or not set(chip) <= {ref, alt}:
            continue
        seq = seq_letters(ref, alt, gt).replace("/", "")
        total += 1
        same = sorted(chip) == sorted(seq)
        agree += same
        if chip != ref + ref:
            nonref_total += 1
            nonref_agree += same
    if not total:
        return None
    return {"total": total, "agree": agree / total,
            "nonref_total": nonref_total, "nonref_agree": nonref_agree / nonref_total if nonref_total else None}


# =============================================================================
# extras: polygenic scores, ancestry, maternal line, medication response
# =============================================================================

POP_NAMES = {
    "CEU": "Utah, N. & W. European ancestry", "TSI": "Tuscany, Italy", "FIN": "Finland",
    "GBR": "England & Scotland", "IBS": "Spain (Iberian)",
    "YRI": "Yoruba, Nigeria", "LWK": "Luhya, Kenya", "GWD": "Gambia", "MSL": "Mende, Sierra Leone",
    "ESN": "Esan, Nigeria", "ASW": "African Americans, SW USA", "ACB": "African Caribbean, Barbados",
    "MXL": "Mexican ancestry, LA", "PUR": "Puerto Rico", "CLM": "Medellín, Colombia", "PEL": "Lima, Peru",
    "CHB": "Beijing, China", "JPT": "Tokyo, Japan", "CHS": "Southern China", "CDX": "Dai, China",
    "KHV": "Kinh, Vietnam", "GIH": "Gujarati, Houston", "PJL": "Punjabi, Pakistan",
    "BEB": "Bengali, Bangladesh", "STU": "Sri Lankan Tamil, UK", "ITU": "Telugu, UK",
}
SUPER_NAMES = {"EUR": "Europe", "AFR": "Africa", "AMR": "Americas", "EAS": "East Asia", "SAS": "South Asia"}
POP_COLORS = {"EUR": "#4C78A8", "AFR": "#E45756", "AMR": "#F2A93B", "EAS": "#54A24B", "SAS": "#B279A2"}
EUR_COLORS = {"CEU": "#4C78A8", "GBR": "#72B7B2", "FIN": "#B279A2", "TSI": "#F2A93B", "IBS": "#E45756"}
HAPLOGREP_ZIP = "https://github.com/genepi/haplogrep3/releases/download/v3.2.2/haplogrep3-3.2.2-linux.zip"
PGS_API = "https://www.pgscatalog.org/rest/"


def have(tool):
    return shutil.which(tool) is not None


def install_hint(tools):
    return ("Install with:  conda install -y -n genome -c conda-forge -c bioconda " + " ".join(tools))


def sample_sheet():
    """1000 Genomes sample -> (population, super-population, sex)."""
    out = {}
    path = P("panel", "samples.panel")
    if os.path.exists(path):
        for line in open(path):
            p = line.split()
            if len(p) >= 4 and p[0] != "sample":
                out[p[0]] = (p[1], p[2], p[3])
    return out


_PLINK_MEM = [None]


def plink_try(args):
    """Runs plink2, lowering its memory request if the Mac can't spare it. Returns (ok, message)."""
    options = [_PLINK_MEM[0]] if _PLINK_MEM[0] else [4000, 2000, 1200, 700]
    last = ""
    for mem in options:
        r = subprocess.run(["plink2", "--threads", str(load_config()["threads"]), "--memory", str(mem)] + args,
                           capture_output=True, text=True)
        out = (r.stdout or "") + (r.stderr or "")
        if r.returncode == 0:
            _PLINK_MEM[0] = mem
            return True, out
        last = out
        if "Out of memory" not in out:
            break
    tail = "\n".join(l for l in last.strip().splitlines() if l.strip())[-600:]
    return False, tail


def plink(args, quiet=True):
    ok, msg = plink_try(args)
    if not ok:
        die("This step failed:\n  plink2 " + " ".join(args) + "\n\n" + msg)
    return msg


def imputed_ready():
    return os.path.exists(P("my_genome.imputed.bcf.csi"))


def refpop_ready():
    return bool(glob.glob(P("refpop", "chr*.bcf.csi")))


def build_plink_sets(force=False):
    """Converts the reference people and your imputed genome to plink2 format, one chromosome at a
    time so it fits in modest memory. Cached; resumes if interrupted. Returns ({chrom: ref}, {chrom: me})."""
    d = P("plink")
    os.makedirs(d, exist_ok=True)
    for old in glob.glob(os.path.join(d, "ref.*")) + glob.glob(os.path.join(d, "me.*")) + \
            glob.glob(os.path.join(d, "ref_concat.bcf")) + glob.glob(os.path.join(d, "stamp")):
        os.remove(old)  # leftovers from the older all-at-once version
    autos = sorted([c for c in (str(i) for i in range(1, 23)) if os.path.exists(refpop_path(c) + ".csi")], key=int)
    stamp = str(os.path.getmtime(P("my_genome.imputed.bcf")))
    refs, mes = {}, {}
    todo = [c for c in autos if force or not os.path.exists(os.path.join(d, "chr%s.ok" % c))
            or open(os.path.join(d, "chr%s.ok" % c)).read() != stamp]
    if todo:
        say("  preparing the reference people and your genome for scoring, one chromosome at a time "
            "(one-off, about 20-40 minutes)…")
    for c in autos:
        ref, me = os.path.join(d, "ref_chr%s" % c), os.path.join(d, "me_chr%s" % c)
        if c in todo:
            plink(["--bcf", refpop_path(c), "--set-all-var-ids", "@:#", "--rm-dup", "force-first",
                   "--make-pgen", "--out", ref])
            tmp = os.path.join(d, "me_chr%s.tmp.bcf" % c)
            run(["bcftools", "view", "-r", c, "-Ob", "-o", tmp, P("my_genome.imputed.bcf")])
            plink(["--bcf", tmp, "dosage=DS", "--set-all-var-ids", "@:#", "--rm-dup", "force-first",
                   "--make-pgen", "--out", me])
            os.remove(tmp)
            with open(os.path.join(d, "chr%s.ok" % c), "w") as f:
                f.write(stamp)
            say("    chromosome %s ready" % c)
        refs[c], mes[c] = ref, me
    return refs, mes


def merge_pfiles(parts, out):
    """Joins per-chromosome sets (small ones only) into one."""
    parts = [p for p in parts if os.path.exists(p + ".pgen")]
    if len(parts) == 1:
        plink(["--pfile", parts[0], "--make-pgen", "--out", out])
        return out
    lst = out + ".list"
    with open(lst, "w") as f:
        f.write("\n".join(parts) + "\n")
    plink(["--pmerge-list", lst, "pfile", "--make-pgen", "--out", out])
    for f in glob.glob(out + "-merge.*"):
        os.remove(f)
    return out


def keep_file(name, samples):
    path = P("plink", name + ".keep")
    with open(path, "w") as f:
        f.write("#IID\n" + "\n".join(samples) + "\n")
    return path


# ------------------------------------------------------------------ ancestry

def run_pca(refs, mes, keep, tag, n_pcs):
    """PCA on reference people, then projects everyone (you included) the same way."""
    base = P("plink", "pca_" + tag)
    ref_parts, me_parts = [], []
    for c in refs:
        part = "%s_chr%s" % (base, c)
        plink(["--pfile", refs[c], "--keep", keep, "--maf", "0.05", "--thin-count", "12000", "--seed", "1",
               "--indep-pairwise", "200kb", "0.2", "--out", part])
        plink(["--pfile", refs[c], "--keep", keep, "--extract", part + ".prune.in", "--make-pgen",
               "--out", part + "_ref"])
        ok, _ = plink_try(["--pfile", mes[c], "--extract", part + ".prune.in", "--make-pgen", "--out", part + "_me"])
        ref_parts.append(part + "_ref")
        if ok:
            me_parts.append(part + "_me")
    ref = merge_pfiles(ref_parts, base + "_ref")
    me = merge_pfiles(me_parts, base + "_me")
    for f in glob.glob(base + "_chr*"):
        os.remove(f)
    plink(["--pfile", ref, "--freq", "counts", "--pca", str(n_pcs), "allele-wts", "--out", base])
    proj = []
    for who, pf in (("ref", ref), ("me", me)):
        out = base + "_proj_" + who
        plink(["--pfile", pf, "--read-freq", base + ".acount",
               "--score", base + ".eigenvec.allele", "2", "5", "header-read", "no-mean-imputation",
               "variance-standardize", "--score-col-nums", "6-%d" % (5 + n_pcs), "--out", out])
        proj.append(out + ".sscore")
    coords = {}
    for path in proj:
        with open(path) as f:
            head = f.readline().lstrip("#").split()
            i_id = head.index("IID")
            pcs = [i for i, h in enumerate(head) if re.match(r"PC\d+_AVG$", h)]
            for line in f:
                p = line.split()
                coords[p[i_id]] = [float(p[i]) for i in pcs]
    return coords


def ancestry_analysis():
    refs, mes = build_plink_sets()
    sheet = sample_sheet()
    everyone = [s for s in sheet]
    eur = [s for s, v in sheet.items() if v[1] == "EUR"]
    result = {}
    for tag, samples, n in (("world", everyone, 6), ("europe", eur, 4)):
        say("  ancestry: %s comparison…" % tag)
        coords = run_pca(refs, mes, keep_file(tag, samples), tag, n)
        mine = [v for k, v in coords.items() if k not in sheet]
        if not mine:
            continue
        you = mine[0]
        groups = {}
        for s, c in coords.items():
            if s in sheet:
                groups.setdefault(sheet[s][0], []).append(c)
        # distance to each population's centre, scaled by that population's spread
        dists = []
        for pop, cs in groups.items():
            centre = [sum(x[i] for x in cs) / len(cs) for i in range(n)]
            spread = math.sqrt(sum(sum((x[i] - centre[i]) ** 2 for i in range(n)) for x in cs) / len(cs)) or 1
            d = math.sqrt(sum((you[i] - centre[i]) ** 2 for i in range(n)))
            dists.append((d / spread, pop))
        dists.sort()
        result[tag] = {"you": you, "points": {s: c[:2] for s, c in coords.items() if s in sheet},
                       "nearest": dists[:5]}
    with open(P("db", "ancestry.json"), "w") as f:
        json.dump(result, f)
    return result


def svg_scatter(points, sheet, you, colour_of, label_of, title):
    xs = [c[0] for c in points.values()] + [you[0]]
    ys = [c[1] for c in points.values()] + [you[1]]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    W, H, M = 560, 400, 24
    sx = lambda x: M + (x - x0) / ((x1 - x0) or 1) * (W - 2 * M)
    sy = lambda y: H - M - (y - y0) / ((y1 - y0) or 1) * (H - 2 * M)
    out = ["<svg viewBox='0 0 %d %d' role='img' aria-label='%s' style='width:100%%;max-width:%dpx;"
           "background:var(--card);border:1px solid var(--line);border-radius:12px'>" % (W, H, esc(title), W)]
    for s, c in points.items():
        out.append("<circle cx='%.1f' cy='%.1f' r='2.6' fill='%s' fill-opacity='0.55'/>"
                   % (sx(c[0]), sy(c[1]), colour_of(s)))
    out.append("<circle cx='%.1f' cy='%.1f' r='8' fill='none' stroke='var(--fg)' stroke-width='2.5'/>"
               "<circle cx='%.1f' cy='%.1f' r='3' fill='var(--fg)'/>" % (sx(you[0]), sy(you[1]), sx(you[0]), sy(you[1])))
    lx = sx(you[0]) + 12 if sx(you[0]) < W - 60 else sx(you[0]) - 40
    out.append("<text x='%.1f' y='%.1f' font-size='13' font-weight='600' fill='var(--fg)'>you</text>"
               % (lx, sy(you[1]) + 4))
    out.append("</svg>")
    legend = "".join("<span style='margin-right:14px;white-space:nowrap'><span style='display:inline-block;"
                     "width:10px;height:10px;border-radius:50%%;background:%s;margin-right:5px'></span>%s</span>"
                     % (col, esc(lab)) for lab, col in label_of)
    return "".join(out) + "<div class='muted' style='font-size:13px;margin:6px 0 14px'>%s</div>" % legend


# ------------------------------------------------------------------ polygenic scores

def pgs_api(path, params=None):
    return json.loads(curl_text(PGS_API + path, params))


def pgs_candidates(keywords):
    traits = {}
    for kw in keywords:
        try:
            d = pgs_api("trait/search", {"term": kw, "include_children": "1"})
        except Exception:
            continue
        for t in d.get("results", []):
            if any(k in t["label"].lower() for k in keywords):
                traits[t["id"]] = t["label"]
    scores = {}
    for tid in traits:
        try:
            d = pgs_api("score/search", {"trait_id": tid, "include_children": "1", "limit": "250"})
        except Exception:
            continue
        for r in d.get("results", []):
            scores[r["id"]] = r
    return list(scores.values())


def pgs_pick(scores):
    """European-developed, GRCh37-harmonized, preferably genome-wide, most recent."""
    good = []
    for r in scores:
        url37 = ((r.get("ftp_harmonized_scoring_files") or {}).get("GRCh37") or {}).get("positions")
        if not url37:
            continue
        ad = r.get("ancestry_distribution") or {}
        eur = max((ad.get(k, {}).get("dist", {}).get("EUR", 0) for k in ("dev", "gwas")), default=0)
        if eur < 80 or (r.get("variants_number") or 0) < 20:
            continue
        evaluated_eur = ad.get("eval", {}).get("dist", {}).get("EUR", 0) > 0
        method = (r.get("method_name") or "").lower()
        ld_aware = any(k in method for k in ("ldpred", "prs-cs", "prscs", "sbayes", "sdpr", "lassosum", "megaprs",
                                             "clump", "c+t", "penalized", "bigstatsr", "snpnet", "prsice"))
        good.append(((r.get("variants_number") or 0) >= 1000, ld_aware, evaluated_eur,
                     r["publication"].get("date_publication") or "", r.get("variants_number") or 0, r))
    good.sort(key=lambda x: x[:5], reverse=True)
    return good[0][5] if good else None


def pgs_performance(pgs_id):
    """An 'odds per standard deviation' from a European evaluation, if the catalog has one."""
    try:
        d = pgs_api("performance/search", {"pgs_id": pgs_id})
    except Exception:
        return None
    for r in d.get("results", []):
        anc = [s.get("ancestry_broad", "") for s in r.get("sampleset", {}).get("samples", [])]
        if not any("European" in a for a in anc):
            continue
        for m in r.get("performance_metrics", {}).get("effect_sizes", []):
            name = (m.get("name_short") or "") + " " + (m.get("name_long") or "")
            if m.get("estimate") and ("OR" in name or "Odds" in name) and ("SD" in name or "standard" in name.lower()):
                return float(m["estimate"])
    return None


def pgs_scorefile(score):
    """Downloads a harmonized score and writes a plink2 score file. Returns (path, n_variants)."""
    pid = score["id"]
    raw = P("db", "pgs", pid + "_hmPOS_GRCh37.txt.gz")
    out = P("db", "pgs", pid + ".plink.tsv")
    if os.path.exists(out):
        return out, sum(1 for _ in open(out)) - 1
    if not os.path.exists(raw):
        url37 = score["ftp_harmonized_scoring_files"]["GRCh37"]["positions"].replace("ftp://", "https://")
        download_checked(url37, raw)
    out = P("db", "pgs", pid + ".plink.tsv")
    if os.path.exists(out):
        return out, sum(1 for _ in open(out)) - 1
    weight_log = False
    seen, n = set(), 0
    with gzip.open(raw, "rt") as f, open(out, "w") as o:
        o.write("ID\tA1\tW\n")
        header = None
        for line in f:
            if line.startswith("#"):
                if line.lower().startswith("#weight_type=") and line.split("=", 1)[1].strip().upper() in ("OR", "HR"):
                    weight_log = True
                continue
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = {h: i for i, h in enumerate(p)}
                continue
            try:
                c, pos = p[header["hm_chr"]], p[header["hm_pos"]]
                ea, w = p[header["effect_allele"]].upper(), float(p[header["effect_weight"]])
            except (KeyError, ValueError, IndexError):
                continue
            if not c or not pos or len(ea) != 1 or ea not in "ACGT" or c not in CHROM_IDX or c in ("X", "Y", "MT"):
                continue
            vid = "%s:%s" % (c, pos)
            if vid in seen:
                continue
            seen.add(vid)
            o.write("%s\t%s\t%r\n" % (vid, ea, math.log(w) if weight_log and w > 0 else w))
            n += 1
    os.remove(raw)  # the converted file is all that's needed
    return out, n


def score_people(pfiles, scorefile, keep, out, extract=None):
    """Scores people chromosome by chromosome and adds up the totals.
    Returns ({person: score}, {chrom: file listing the variants used})."""
    vals, used = {}, {}
    for c, pf in pfiles.items():
        if extract is not None and c not in extract:
            continue
        o = "%s_chr%s" % (out, c)
        args = ["--pfile", pf]
        if keep:
            args += ["--keep", keep]
        if extract is not None:
            args += ["--extract", extract[c]]
        ok, _ = plink_try(args + ["--score", scorefile, "1", "2", "3", "header", "no-mean-imputation",
                                  "list-variants", "cols=+scoresums", "--out", o])
        if not ok:
            continue  # no variants of this score on this chromosome
        with open(o + ".sscore") as f:
            head = f.readline().lstrip("#").split()
            i_id, i_sum = head.index("IID"), head.index("SCORE1_SUM")
            for line in f:
                p = line.split()
                vals[p[i_id]] = vals.get(p[i_id], 0.0) + float(p[i_sum])
        used[c] = o + ".sscore.vars"
    return vals, used


def pgs_analysis(cfg, force=False):
    refs, mes = build_plink_sets()
    sheet = sample_sheet()
    european = cfg.get("ancestry") == "european"
    base_pop = [s for s, v in sheet.items() if v[1] == "EUR"] if european else list(sheet)
    base_name = "Europeans" if european else "people"
    keep = keep_file("pgsbase", base_pop)
    os.makedirs(P("db", "pgs"), exist_ok=True)
    chosen = cfg.setdefault("pgs", {})
    old = {}
    if os.path.exists(P("db", "pgs.json")) and not force:
        old = json.load(open(P("db", "pgs.json")))
    todo = [(name, w.get("gwas", []), None, watch_part(name), "modest", True) for name, w in cfg["watch"].items()]
    watched_low = {n.lower() for n in cfg["watch"]}
    todo += [(name, kws, pid, part, ev, False) for name, (part, ev, kws, pid) in PGS_PANEL.items()
             if name.lower() not in watched_low]
    results = {}
    for i, (name, kws, default_id, part, ev, is_watched) in enumerate(todo, 1):
        pid = chosen.get(name) or default_id
        if pid == "none":
            continue
        prev = old.get(name)
        if prev and prev.get("pct") is not None and (not pid or prev.get("id") == pid):
            results[name] = dict(prev, part=part, ev=ev, watched=is_watched)
            continue  # already scored; use  extras --force  to redo
        if prev and prev.get("none") and prev.get("checked") == today():
            results[name] = prev
            continue
        try:
            if pid:
                score = pgs_api("score/" + pid)
            else:
                score = pgs_pick(pgs_candidates([k.lower() for k in kws]))
                if not score:
                    results[name] = {"none": True, "part": part, "ev": ev, "watched": is_watched, "checked": today()}
                    continue
            say("  scores (%d/%d): %s — %s (%s variants)" % (i, len(todo), name, score["id"],
                                                          format(score.get("variants_number") or 0, ",")))
            sf, n_total = pgs_scorefile(score)
            out = P("plink", "pgs_" + score["id"])
            ref_vals, used = score_people(refs, sf, keep, out + "_ref")
            if not ref_vals:
                raise RuntimeError("none of the score's variants are in the reference data")
            me_vals, _ = score_people(mes, sf, None, out + "_me", extract=used)
            mine = list(me_vals.values())[0]
            n_used = sum(sum(1 for _ in open(v)) for v in used.values())
            vals = sorted(ref_vals.values())
            mean = sum(vals) / len(vals)
            sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals)) or 1
            pct = 100.0 * sum(1 for v in vals if v < mine) / len(vals)
            results[name] = {
                "id": score["id"], "name": score.get("name"), "trait": score.get("trait_reported"),
                "author": score["publication"].get("firstauthor"), "year": (score["publication"].get("date_publication") or "")[:4],
                "method": score.get("method_name"), "n_total": n_total, "n_used": n_used,
                "pct": pct, "z": (mine - mean) / sd, "base_name": base_name, "base_n": len(vals),
                "or_per_sd": pgs_performance(score["id"]), "part": part, "ev": ev, "watched": is_watched,
            }
            for f_ in glob.glob(out + "_*"):
                if not f_.endswith(".vars"):
                    os.remove(f_)
            with open(P("db", "pgs.json"), "w") as f:
                json.dump(dict(old, **results), f)  # save as we go, so an interruption loses little
        except (Exception, SystemExit) as e:
            msg = str(e).strip().splitlines()
            say("  scores: %s — skipped (%s)" % (name, (msg[-1] if msg else "error")[:160]))
    save_config(cfg)
    with open(P("db", "pgs.json"), "w") as f:
        json.dump(results, f)
    return results


def ordinal(n):
    n = int(round(n))
    n = min(max(n, 1), 99)
    return "%d%s" % (n, "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th"))


def plain_pct(pct, or_sd, z):
    words = ("in the lower part" if pct < 30 else "in the middle range" if pct <= 70 else
             "in the higher part" if pct <= 90 else "near the top")
    if pct < 1:
        txt = "Your score is lower than all or almost all of the reference people."
    elif pct > 99:
        txt = "Your score is higher than all or almost all of the reference people."
    else:
        txt = "Your score is %s: higher than about %d in 100 of the reference people." % (words, round(pct))
    if pct < 2 or pct > 98:
        txt += (" Values this extreme are worth a second look: they can come from one region with a strong effect, "
                "or from a score built without correcting for linked variants.")
    if or_sd and or_sd > 1:
        odds = or_sd ** z
        txt += (" For this score, each step of one standard deviation changes the odds about %.1f×, so a score "
                "like yours corresponds to roughly %.1f× typical odds." % (or_sd, odds))
    return txt


# ------------------------------------------------------------------ maternal line (mtDNA)

def haplogrep_dir(cfg):
    d = P("tools", "haplogrep3")
    jar = os.path.join(d, "haplogrep3.jar")
    if not os.path.exists(jar):
        os.makedirs(d, exist_ok=True)
        z = os.path.join(d, "haplogrep3.zip")
        download_checked(cfg.get("urls", {}).get("haplogrep", HAPLOGREP_ZIP), z)
        with zipfile.ZipFile(z) as zz:
            zz.extract("haplogrep3.jar", d)
        os.remove(z)
    with open(os.path.join(d, "haplogrep3.yaml"), "w") as f:
        f.write("repositories:\n  - https://raw.githubusercontent.com/genepi/haplogrep-trees/main/trees.json\n"
                "phylotrees:\n- phylotree-rcrs@17.2\n")
    return d


def mt_mapped_reads(cram):
    r = subprocess.run(["samtools", "idxstats", cram], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if p[0] in ("MT", "chrM", "M") and len(p) >= 3:
            return p[0], int(p[2])
    return None, 0


def rescue_mt_reads(cram, fa):
    """When the sequencing company's map left out mitochondrial DNA, its reads end up 'unplaced'.
    This matches the unplaced reads against the mitochondrial sequence itself."""
    if not have("minimap2"):
        die("Your file has no reads placed on the mitochondrial DNA, but they can be recovered from the "
            "unplaced reads. That needs one more small tool:\n  " + install_hint(["minimap2"]) +
            "\nthen run the same command again.")
    mt_fa = P("ref", "mt_rcrs.fa")
    if not os.path.exists(mt_fa + ".fai"):
        with open(mt_fa, "w") as f:
            f.write(run(["samtools", "faidx", fa, "MT"]))
        run(["samtools", "faidx", mt_fa])
    bam = P("db", "mt_rescued.bam")
    say("  maternal line: your file has no reads placed on the mitochondrial DNA; recovering them from the "
        "unplaced reads (it reads through the whole file: 15-40 min)…")
    t = str(load_config()["threads"])
    p1 = subprocess.Popen(["samtools", "view", "-u", "-f", "4", "--reference", fa, cram],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "fastq", "-"], stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p3 = subprocess.Popen(["minimap2", "-ax", "sr", "-t", t, mt_fa, "-"], stdin=p2.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2.stdout.close()
    p4 = subprocess.Popen(["samtools", "view", "-u", "-F", "2308", "-q", "20", "-"], stdin=p3.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p3.stdout.close()
    r = subprocess.run(["samtools", "sort", "-o", bam, "-"], stdin=p4.stdout, capture_output=True, text=True)
    p4.stdout.close()
    for p in (p1, p2, p3, p4):
        p.wait()
    if r.returncode != 0:
        die("Recovering the mitochondrial reads failed:\n" + r.stderr[-500:])
    run(["samtools", "index", bam])
    return bam, mt_fa


def maternal_analysis(cfg, cram):
    fa = fasta_path()
    out_vcf = P("db", "mtdna.vcf")
    contig, n_mapped = mt_mapped_reads(cram)
    if n_mapped >= 500:
        reads, ref, contig = cram, fa, contig
        rescued = False
    else:
        reads, ref = rescue_mt_reads(cram, fa)
        contig, rescued = "MT", True
    depth = None
    r = subprocess.run(["samtools", "coverage", "-r", contig, "--reference", ref, reads], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith(contig + "\t"):
            depth = float(line.split("\t")[6])
    if not depth or depth < 8:
        with open(P("db", "maternal.json"), "w") as f:
            json.dump({"none": True, "depth": depth or 0}, f)
        die("Too few mitochondrial reads to read your maternal line (average depth %s). This can happen when "
            "a sequencing service filters them out completely." % ("%.1f×" % depth if depth else "0"))
    say("  maternal line: reading your mitochondrial DNA (average depth %.0f×)…" % depth)
    p1 = subprocess.Popen(["bcftools", "mpileup", "-r", contig, "-f", ref, "-d", "10000", "-q", "20", "-Q", "20",
                           "-a", "AD,DP", reads], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["bcftools", "call", "--ploidy", "1", "-mv"], stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    p1.stdout.close()
    with open(out_vcf, "w") as f:
        p3 = subprocess.Popen(["bcftools", "view", "-v", "snps", "-i", "QUAL>=30"], stdin=p2.stdout, stdout=f,
                              stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.wait()
    n_var = sum(1 for l in open(out_vcf) if not l.startswith("#"))
    d = haplogrep_dir(cfg)
    res = os.path.join(d, "result.txt")
    if os.path.exists(res):
        os.remove(res)
    r = subprocess.run(["java", "-Xmx2G", "-jar", "haplogrep3.jar", "classify", "--tree", "phylotree-rcrs@17.2",
                        "--in", out_vcf, "--out", res, "--extend-report"], cwd=d, capture_output=True, text=True)
    rows = list(csv.reader(open(res), delimiter="\t")) if os.path.exists(res) else []
    if r.returncode != 0 or len(rows) < 2:
        msg = ((r.stderr or "") + (r.stdout or "")).strip().splitlines()
        die("Haplogrep couldn't classify your mitochondrial DNA (%d differences found).\n%s"
            % (n_var, msg[-1][:200] if msg else ""))
    head = [h.strip('"') for h in rows[0]]
    val = {h: v.strip('"') for h, v in zip(head, rows[1])}
    result = {"haplogroup": val.get("Haplogroup"), "quality": to_float(val.get("Quality")),
              "n_variants": n_var, "depth": depth, "rescued": rescued,
              "not_found": val.get("Not_Found_Polys", ""), "remaining": val.get("Remaining_Polys", "")}
    with open(P("db", "maternal.json"), "w") as f:
        json.dump(result, f)
    return result


HAPLO_NOTES = [
    ("H", "the most common maternal lineage in Europe; it spread widely after the last Ice Age"),
    ("U5", "one of Europe's oldest lineages, carried by Ice Age hunter-gatherers"),
    ("U6", "a North African lineage (Maghreb)"),
    ("U", "an old Eurasian lineage with many European branches"),
    ("V", "a Western European lineage, linked to Ice Age refuges in southwestern Europe"),
    ("HV", "the parent group of H and V, found in Europe and the Near East"),
    ("J", "a lineage that arrived in Europe largely with the first farmers from the Near East"),
    ("T", "a lineage that arrived in Europe largely with the first farmers from the Near East"),
    ("K", "a branch of U, common across Europe and the Near East (and among Ashkenazi Jews)"),
    ("I", "a widespread but less common European lineage"),
    ("W", "a lineage found at low frequencies from Europe to South Asia"),
    ("X", "a rare lineage found in Europe, the Near East and Native North Americans"),
    ("L", "a lineage from sub-Saharan Africa, home to the oldest branches of the human maternal tree"),
    ("M1", "a North and East African lineage"),
    ("N", "an ancestral Eurasian lineage"),
]


def haplo_note(hg):
    for prefix, note in HAPLO_NOTES:
        if hg and hg.upper().startswith(prefix):
            return prefix, note
    return "", ""


# ------------------------------------------------------------------ medication response

# GRCh37 positions verified against ClinVar/dbSNP. "allele" is the letter whose copies are counted.
PGX_SITES = {
    "rs4244285": ("10", 96541616, "A"), "rs4986893": ("10", 96540410, "A"), "rs12248560": ("10", 96521657, "T"),
    "rs1799853": ("10", 96702047, "T"), "rs1057910": ("10", 96741053, "C"), "rs9923231": ("16", 31107689, "T"),
    "rs2108622": ("19", 15990431, "T"), "rs4149056": ("12", 21331549, "C"), "rs2231142": ("4", 89052323, "T"),
    "rs3918290": ("1", 97915614, "T"), "rs55886062": ("1", 97981343, "C"), "rs67376798": ("1", 97547947, "A"),
    "rs75017182": ("1", 98045449, "C"), "rs1800460": ("6", 18139228, "T"), "rs1142345": ("6", 18130918, "C"),
    "rs1800462": ("6", 18143955, "G"), "rs116855232": ("13", 48619855, "T"), "rs776746": ("7", 99270539, "T"),
    "rs887829": ("2", 234668570, "T"), "rs3892097": ("22", 42524947, "T"), "rs6025": ("1", 169519049, "T"),
    "rs1799963": ("11", 46761055, "A"), "rs5030868": ("X", 153762634, "A"), "rs1050828": ("X", 153764217, "T"),
    "rs429358": ("19", 45411941, "C"), "rs7412": ("19", 45412079, "T"),
    "rs1065852": ("22", 42526694, "A"), "rs28371725": ("22", 42523805, "T"),
    "rs3745274": ("19", 41512841, "T"), "rs28399499": ("19", 41518221, "C"),
}


def pgx_genotypes(genome, sites=None):
    """rsID -> {"dose": copies of the counted allele, "gt": text, "conf": 0-1 or None, "src", "raf"}."""
    sites = sites or PGX_SITES
    out = {}
    if imputed_ready():
        regions = ",".join("%s:%d" % (c, p) for c, p, a in sites.values())
        txt = run(["bcftools", "query", "-r", regions, "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/RAF[\t%GT\t%DS\t%GP]\n",
                   P("my_genome.imputed.bcf")])
        bypos = {}
        for line in txt.splitlines():
            f = line.split("\t")
            bypos[(f[0], int(f[1]))] = f
        for rs, (c, pos, allele) in sites.items():
            f = bypos.get((c, pos))
            if not f or allele not in (f[2], f[3]):
                continue
            ds = float(f[6])
            conf = max(float(x) for x in f[7].split(","))
            raf = to_float(f[4])
            out[rs] = {"dose": ds if allele == f[3] else 2 - ds, "gt": seq_letters(f[2], f[3], f[5]),
                       "conf": conf, "src": "sequencing",
                       "raf": (raf if allele == f[3] else 1 - raf) if raf is not None else None}
    if genome.n_chip:
        nums = [int(r[2:]) for r in sites]
        for rs_num, gt in genome.con.execute("SELECT rs, gt FROM chip WHERE rs IN (%s)" % ",".join("?" * len(nums)), nums):
            rs = "rs%d" % rs_num
            allele = sites[rs][2]
            # MyHeritage reports the forward strand; only switch strands when the letters can't match otherwise
            if allele not in gt and COMP[allele] in gt and not set(gt) <= {allele, COMP[allele]}:
                allele = COMP[allele]
            chip = {"dose": float(gt.count(allele)), "gt": gt[0] + "/" + gt[1], "conf": None, "src": "chip", "raf": None}
            if rs not in out:
                out[rs] = chip
            else:
                out[rs]["chip_gt"] = chip["gt"]
                if round(out[rs]["dose"]) != chip["dose"]:
                    out[rs]["disagree"] = True
    return out


def pgx_analysis(genome):
    g = pgx_genotypes(genome)

    def cnt(rs):
        return int(round(g[rs]["dose"])) if rs in g else None

    def known(*rss):
        return all(r in g for r in rss)

    rows = []

    def add(gene, result, meaning, drugs, rss, note=""):
        conf = [g[r]["conf"] for r in rss if r in g and g[r]["conf"] is not None]
        low = [r for r in rss if r in g and g[r]["conf"] is not None and g[r]["conf"] < MIN_CONF]
        rare = [r for r in rss if r in g and g[r].get("raf") is not None and g[r]["raf"] < 0.01]
        dis = [r for r in rss if r in g and g[r].get("disagree")]
        flags = []
        if low:
            flags.append("uncertain imputation at " + ", ".join(low))
        if dis:
            flags.append("%s and your imputation disagree at %s" % (chip_word(), ", ".join(dis)))
        if rare:
            flags.append("includes rare variants that imputation can miss")
        rows.append({"gene": gene, "result": result, "meaning": meaning, "drugs": drugs,
                     "variants": " · ".join("%s %s" % (r, g[r]["gt"]) for r in rss if r in g) or "no data",
                     "flags": flags + ([note] if note else []), "conf": min(conf) if conf else None})

    if known("rs4244285", "rs12248560"):
        nf = cnt("rs4244285") + (cnt("rs4986893") or 0)
        inc = cnt("rs12248560")
        res = ("Poor metabolizer" if nf >= 2 else "Intermediate metabolizer" if nf == 1 else
               "Ultrarapid metabolizer" if inc == 2 else "Rapid metabolizer" if inc == 1 else "Normal metabolizer")
        add("CYP2C19", res, "How fast you clear several antidepressants, a heart drug and stomach-acid drugs.",
            "clopidogrel; citalopram, escitalopram, sertraline; omeprazole and other PPIs; voriconazole",
            ["rs4244285", "rs4986893", "rs12248560"])
    if known("rs1799853", "rs1057910"):
        score = 2 - 0.5 * cnt("rs1799853") - cnt("rs1057910")
        res = "Normal metabolizer" if score >= 2 else "Intermediate metabolizer" if score >= 1 else "Poor metabolizer"
        add("CYP2C9", res + " (activity score %g)" % score, "Clearance of common painkillers and blood thinners.",
            "ibuprofen, celecoxib, meloxicam and other NSAIDs; warfarin; phenytoin", ["rs1799853", "rs1057910"])
    if known("rs9923231"):
        c = cnt("rs9923231")
        add("VKORC1", ["Typical warfarin sensitivity", "Increased warfarin sensitivity",
                       "High warfarin sensitivity"][c], "How much warfarin is needed; with CYP2C9 and CYP4F2.",
            "warfarin", ["rs9923231", "rs2108622"])
    if known("rs4149056"):
        c = cnt("rs4149056")
        add("SLCO1B1", ["Normal function", "Decreased function", "Poor function"][c],
            "How well the liver takes up statins; lower function raises the risk of muscle side effects.",
            "simvastatin (most affected), atorvastatin, rosuvastatin and other statins", ["rs4149056"])
    if known("rs2231142"):
        c = cnt("rs2231142")
        add("ABCG2", ["Normal function", "Decreased function", "Poor function"][c],
            "A transporter for rosuvastatin and uric acid.", "rosuvastatin; allopurinol (gout)", ["rs2231142"])
    dpyd = {"rs3918290": 1.0, "rs55886062": 1.0, "rs67376798": 0.5, "rs75017182": 0.5}
    if any(r in g for r in dpyd):
        loss = sum((cnt(r) or 0) * v for r, v in dpyd.items())
        score = max(0.0, 2 - loss)
        res = "Normal metabolizer" if score >= 2 else "Intermediate metabolizer" if score >= 1 else "Poor metabolizer"
        add("DPYD", res + " (activity score %g)" % score, "Breakdown of fluoropyrimidine chemotherapy.",
            "5-fluorouracil, capecitabine, tegafur", list(dpyd),
            "A clinical DPYD test before these drugs is recommended regardless of this result.")
    if any(r in g for r in ("rs1800460", "rs1142345", "rs1800462")):
        nf = max(cnt("rs1800460") or 0, cnt("rs1142345") or 0) + (cnt("rs1800462") or 0)
        res = ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][min(nf, 2)]
        add("TPMT", res, "Breakdown of thiopurine drugs (immune suppressants).",
            "azathioprine, mercaptopurine, thioguanine", ["rs1800460", "rs1142345", "rs1800462"])
    if known("rs116855232"):
        add("NUDT15", ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][cnt("rs116855232")],
            "Also affects thiopurine tolerance (variant mainly found in East Asian and Hispanic ancestry).",
            "azathioprine, mercaptopurine, thioguanine", ["rs116855232"])
    if known("rs776746"):
        add("CYP3A5", ["Non-expresser (the usual result in Europeans)", "Intermediate metabolizer",
                       "Normal metabolizer (expresser)"][cnt("rs776746")],
            "Mainly matters for one transplant drug.", "tacrolimus", ["rs776746"])
    if known("rs887829"):
        c = cnt("rs887829")
        add("UGT1A1", ["Normal function", "Intermediate (likely one *28 copy)",
                       "Poor function (likely *28/*28 — Gilbert's syndrome)"][c],
            "Bilirubin processing; two reduced copies is Gilbert's syndrome, a harmless cause of mild jaundice.",
            "irinotecan (chemo), atazanavir, belinostat", ["rs887829"],
            "Measured through a tag variant, not the *28 repeat itself.")
    if known("rs3892097"):
        n4 = cnt("rs3892097")
        n10 = max(0, (cnt("rs1065852") or 0) - n4)   # 100C>T sits on *4 as well; extra copies mean *10
        n41 = cnt("rs28371725") or 0
        score = max(0.0, 2 - n4 - 0.75 * n10 - 0.5 * n41)
        res = ("Poor metabolizer" if score == 0 else "Intermediate metabolizer" if score <= 1 else "Normal metabolizer")
        found = [x for x, n in (("*4", n4), ("*10", n10), ("*41", n41)) for _ in range(n)]
        add("CYP2D6 (partial)", "%s (activity score %g%s)" % (res, score, "; " + ", ".join(found) if found else ""),
            "Handles many antidepressants, antipsychotics and painkillers. Only *4, *10 and *41 are read here: the "
            "small deletions *3 and *6, the whole-gene deletion *5 and extra gene copies can't be seen in this data, "
            "so a clinical test can still change this result (towards poor or ultrarapid).",
            "codeine, tramadol, tamoxifen; paroxetine, fluvoxamine, venlafaxine, vortioxetine, amitriptyline, "
            "nortriptyline; aripiprazole, risperidone; atomoxetine; ondansetron",
            ["rs3892097", "rs1065852", "rs28371725"])
    if known("rs3745274"):
        dec = cnt("rs3745274") + (cnt("rs28399499") or 0)
        add("CYP2B6", ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][min(dec, 2)],
            "Part of how sertraline (and the HIV drug efavirenz) is cleared; included in the 2023 CPIC antidepressant guideline.",
            "sertraline; efavirenz; methadone; bupropion (partly)", ["rs3745274", "rs28399499"])
    clot = []
    if known("rs6025") and cnt("rs6025") > 0:
        clot.append("Factor V Leiden (%d cop%s)" % (cnt("rs6025"), "y" if cnt("rs6025") == 1 else "ies"))
    if known("rs1799963") and cnt("rs1799963") > 0:
        clot.append("Prothrombin G20210A (%d cop%s)" % (cnt("rs1799963"), "y" if cnt("rs1799963") == 1 else "ies"))
    if known("rs6025") or known("rs1799963"):
        add("F5 / F2 (clotting)", ", ".join(clot) if clot else "Neither clotting variant found",
            "These raise the risk of blood clots, especially combined with estrogen-containing contraceptives, "
            "hormone therapy, pregnancy or long immobility.",
            "combined contraceptive pill, patch or ring; estrogen hormone therapy", ["rs6025", "rs1799963"],
            "Worth telling a doctor before starting estrogen-containing contraception if a copy is found."
            if clot else "")
    if known("rs5030868") or known("rs1050828"):
        c = (cnt("rs5030868") or 0) + (cnt("rs1050828") or 0)
        male = load_config().get("sex") == "male"
        add("G6PD (X chromosome)", (["No common deficiency variant found", "Deficiency variant on your X — likely "
                                     "G6PD deficient", "Deficiency variant on your X — likely G6PD deficient"]
                                    if male else
                                    ["No common deficiency variant found", "Carrier of a deficiency variant",
                                     "Two deficiency variants"])[min(c, 2)],
            "Red blood cells' protection against oxidative stress; deficiency can cause anaemia with certain "
            "drugs or fava beans. Women with one copy usually have mild or no effects.",
            "primaquine, tafenoquine, rasburicase, dapsone, nitrofurantoin, methylene blue; fava beans",
            ["rs5030868", "rs1050828"])
    apoe = None
    if known("rs429358", "rs7412"):
        c4, c2 = cnt("rs429358"), cnt("rs7412")
        alleles = ["ε4"] * c4 + ["ε2"] * c2
        alleles += ["ε3"] * (2 - len(alleles))
        apoe = {"type": "/".join(sorted(alleles[:2], key=lambda a: "ε2ε3ε4".index(a))),
                "variants": "rs429358 %s · rs7412 %s" % (g["rs429358"]["gt"], g["rs7412"]["gt"]),
                "conf": min(x for x in (g["rs429358"]["conf"], g["rs7412"]["conf"]) if x is not None)
                if g["rs429358"]["conf"] is not None else None}
    result = {"rows": rows, "apoe": apoe}
    with open(P("db", "pgx.json"), "w") as f:
        json.dump(result, f)
    return result


APOE_TEXT = {
    "ε3/ε3": "the most common type (about 60% of Europeans): typical Alzheimer's risk.",
    "ε2/ε3": "slightly lower than typical Alzheimer's risk.",
    "ε2/ε2": "lower than typical Alzheimer's risk; linked to a rare cholesterol disorder in some people.",
    "ε3/ε4": "one ε4 copy (about 1 in 4 Europeans): Alzheimer's risk roughly 2–3 times typical. Most people "
             "with one copy never develop Alzheimer's.",
    "ε2/ε4": "one ε4 and one ε2 copy: roughly typical to slightly raised risk.",
    "ε4/ε4": "two ε4 copies (about 2 in 100 people): Alzheimer's risk clearly raised, roughly 8–12 times "
             "typical, though not a certainty. Worth discussing with a doctor or genetic counsellor.",
}


def cmd_extras(a, cfg=None):
    cfg = cfg or load_config()
    genome = Genome()
    say("Medication response and established traits…")
    pgx_analysis(genome)
    trait_analysis(genome)
    if not imputed_ready():
        say("Percentiles, ancestry and maternal line need your imputed genome — skipped for now.")
        return
    missing = [t for t in ("plink2",) if not have(t)]
    if missing:
        say("Percentiles and ancestry need plink2. " + install_hint(["plink2"]))
    elif not refpop_ready():
        say("Percentiles and ancestry need the reference people kept during imputation (newer version). "
            "Run:  python3 genome_lab.py reference")
    else:
        try:
            ancestry_analysis()
        except SystemExit as e:
            say("Ancestry check failed: %s" % str(e).strip().splitlines()[-1])
        try:
            pgs_analysis(cfg, force=getattr(a, "force", False) and getattr(a, "cmd", "") == "extras")
        except SystemExit as e:
            say("Scores failed: %s" % str(e).strip().splitlines()[-1])
    cram = meta_get(db(), "cram_path")
    bam = P("db", "mt_rescued.bam")
    if not os.path.exists(P("db", "maternal.json")) and os.path.exists(bam + ".bai") and \
            subprocess.run(["samtools", "view", "-c", bam], capture_output=True, text=True).stdout.strip() == "0":
        with open(P("db", "maternal.json"), "w") as f:
            json.dump({"none": True, "depth": 0}, f)  # an earlier attempt found no mitochondrial reads
    if not os.path.exists(P("db", "maternal.json")) or getattr(a, "force", False):
        if not have("java"):
            say("Maternal line needs Java. " + install_hint(["openjdk=17"]))
        elif not cram or not os.path.exists(cram):
            say("Maternal line needs your CRAM file. Run:  python3 genome_lab.py maternal PATH_TO_CRAM")
        else:
            try:
                maternal_analysis(cfg, cram)
            except (Exception, SystemExit) as e:
                say("Maternal line failed: %s" % str(e).strip())


def cmd_maternal(a):
    cfg = load_config()
    cram = os.path.abspath(os.path.expanduser(a.cram))
    if not os.path.exists(cram):
        die("Can't find " + cram)
    meta_set(db(), "cram_path", cram)
    if not have("java"):
        die("Needs Java. " + install_hint(["openjdk=17"]))
    r = maternal_analysis(cfg, cram)
    say("Your maternal haplogroup: %s (quality %.2f). It's in your next report." % (r["haplogroup"], r["quality"] or 0))


def cmd_reference(a):
    """Rebuilds the reference people for chromosomes imputed with an older version."""
    cfg = load_config()
    t = str(cfg["threads"])
    try:
        listing = curl_text(url(cfg, "panel_dir"))
    except Exception:
        listing = ""
    todo = [c for c in (str(i) for i in range(1, 23)) if not os.path.exists(refpop_path(c) + ".csi")]
    say("Keeping the reference people for %d chromosome(s); this downloads the panel again." % len(todo))
    for c in todo:
        bcf = panel_bcf(c, cfg, listing)
        keep_refpop(c, bcf, t)
        for f in (bcf, bcf + ".csi"):
            if os.path.exists(f):
                os.remove(f)
        say("  chromosome %s done" % c)


def pgs_set(a):
    cfg = load_config()
    if a.name not in cfg["watch"] and a.name not in PGS_PANEL:
        die("No score called %s. Names: %s" % (a.name, ", ".join(list(cfg["watch"]) + list(PGS_PANEL))))
    cfg.setdefault("pgs", {})[a.name] = a.pgs_id if a.pgs_id != "auto" else None
    if a.pgs_id == "auto":
        cfg["pgs"].pop(a.name)
    save_config(cfg)
    say("Saved. Takes effect in your next  run.")


def badge(ev):
    label, tip = EVIDENCE.get(ev, EVIDENCE["weak"])
    return "<span class='ev ev-%s' title='%s'>%s</span>" % (ev, esc(tip), label)


def load_db_json(name):
    path = P("db", name)
    return json.load(open(path)) if os.path.exists(path) else None


# How to explain each score. kind "risk": lifetime chance for a woman of European ancestry (approximate, from
# population statistics) and a typical odds ratio per standard deviation of the score (used when the PGS Catalog
# has no European evaluation). kind "measure": what a higher score means.
SCORE_INFO = {
    "Breast cancer": ("risk", 0.12, 1.6, ""),
    "Coronary artery disease": ("risk", 0.20, 1.7, "heart attack or narrowed heart arteries over a lifetime"),
    "Schizophrenia": ("risk", 0.007, 1.6, ""),
    "Bipolar disorder": ("risk", 0.015, 1.5, ""),
    "Depression": ("risk", 0.22, 1.3, "at least one episode of major depression in a lifetime"),
    "Insomnia": ("risk", 0.10, 1.2, "long-lasting insomnia"),
    "Migraine": ("risk", 0.25, 1.3, ""),
    "Alzheimer's disease": ("risk", 0.18, 1.5, "by very old age; this score includes the APOE region"),
    "Parkinson's disease": ("risk", 0.013, 1.5, ""),
    "Multiple sclerosis": ("risk", 0.003, 1.7, ""),
    "Type 1 diabetes": ("risk", 0.005, 2.5, "it almost always starts in childhood or young adulthood"),
    "Type 2 diabetes": ("risk", 0.25, 1.6, "lifestyle changes it a lot more than genes"),
    "Hypothyroidism": ("risk", 0.08, 1.4, "an underactive thyroid, easily treated"),
    "Rheumatoid arthritis": ("risk", 0.035, 1.5, ""),
    "Asthma": ("risk", 0.12, 1.3, "at any age"),
    "Childhood-onset asthma": ("risk", 0.08, 1.4, "this is about asthma starting in childhood — if you didn't have "
                               "it as a child, it mostly no longer applies to you"),
    "Psoriasis": ("risk", 0.025, 1.7, ""),
    "Coeliac disease": ("risk", 0.01, 2.5, ""),
    "Age-related macular degeneration": ("risk", 0.06, 2.0, "by old age"),
    "Myopia (short-sightedness)": ("risk", 0.35, 1.3, ""),
    "Endometriosis": ("risk", 0.10, 1.3, ""),
    "ADHD": ("risk", 0.04, 1.3, "diagnosed at some point"),
    "Autism": ("risk", 0.012, 1.3, "diagnosed at some point"),
    "Polycystic ovary syndrome": ("risk", 0.10, 1.3, ""),
    "Height": ("measure", "being taller", "being shorter",
               "This is the best check of whether scores work for you: compare it with your real height. Height "
               "is about 80% genetic and this score captures about half of that."),
    "Body mass index": ("measure", "a higher body weight for your height (BMI)", "a lower BMI",
                        "Genes set a tendency; food, activity and sleep decide much more."),
    "LDL cholesterol": ("measure", "higher LDL (“bad”) cholesterol", "lower LDL cholesterol",
                        "A routine blood test shows your real level."),
    "Age at menopause": ("measure", "a later natural menopause", "an earlier natural menopause",
                         "Shifts the expected age by months, not years, for most people. Smoking and other factors matter too."),
    "Age at first period": ("measure", "a later first period", "an earlier first period",
                            "Genes explain roughly half the variation in when periods start; childhood nutrition and "
                            "body weight explain much of the rest. You can compare it with your own age at your first "
                            "period — another check of how well these scores work for you."),
    "Thyroid hormone (TSH) level": ("measure", "higher TSH (a nudge toward an underactive thyroid)",
                                    "lower TSH (a nudge toward an overactive thyroid)", "A blood test shows your real level."),
    "Chronotype (morning person)": ("measure", "being an evening person (night owl)", "being a morning person (early bird)",
                                    "You can check this one against how you actually sleep."),
}
SCORE_LABEL = {"Chronotype (morning person)": "Chronotype (morning vs evening person)"}


def one_in(p):
    if p >= 0.095:
        return "about %d in 100" % round(p * 100)
    n = 1 / p
    n = round(n, -1) if n >= 30 else round(n)
    if n >= 1000:
        n = round(n, -2)
    return "about 1 in %d" % n


def score_block(name, r, compact=False):
    """Plain explanation of a published score's result. Returns (headline_html, body_html)."""
    info = SCORE_INFO.get(name)
    pct, z = r["pct"], r["z"]
    rank = ("higher than about %d in 100" % round(pct)) if pct >= 50 else ("lower than about %d in 100" % round(100 - pct))
    base = r.get("base_name", "Europeans")
    if info and info[0] == "risk":
        _, base_p, or_default, note = info
        if load_config().get("sex") != "female":
            note = (note + "; " if note else "") + "lifetime figures are for women of European ancestry — men's differ"
        or_sd = r.get("or_per_sd") if r.get("or_per_sd") and 1 < r.get("or_per_sd") < 6 else or_default
        z_ = max(min(z, 3), -3)
        odds = base_p / (1 - base_p) * or_sd ** z_
        mine = odds / (1 + odds)
        head = ("<div class='big'>%s</div><div class='small'>your estimated lifetime chance · typical: %s</div>"
                % (esc(one_in(mine).capitalize()), esc(one_in(base_p))))
        body = ("<p style='margin:6px 0 0'>Your genetic score is %s %s (the %s percentile). That puts your "
                "chance at %s instead of %s — so %s people with a score like yours <strong>never</strong> get it.%s</p>"
                % (rank, esc(base), ordinal(pct), esc(one_in(mine)), esc(one_in(base_p)),
                   ("about %d in 100" % round(100 * (1 - mine))) if mine >= 0.01 else "more than 99 in 100",
                   (" Note: " + esc(note) + ".") if note else ""))
    elif info:
        _, up, down, note = info
        lean = ("about typical" if 30 <= pct <= 70 else ("toward " + up) if pct > 70 else ("toward " + down))
        head = "<div class='big'>%s</div>" % esc(lean[0].upper() + lean[1:])
        body = ("<p style='margin:6px 0 0'>Your genetic score for %s is %s %s (the %s percentile). %s</p>"
                % (esc(up), rank, esc(base), ordinal(pct), esc(note)))
    else:
        head = "<div class='big'>%s</div>" % esc(pct_words(pct).capitalize())
        body = "<p style='margin:6px 0 0'>Your genetic score is %s %s (the %s percentile).</p>" % (rank, esc(base), ordinal(pct))
    if pct < 2 or pct > 98:
        body += ("<p class='muted small'>Scores this extreme are worth a second look: they can come from one region "
                 "with a strong effect, or from a score built without correcting for linked variants.</p>")
    cover = r["n_used"] / r["n_total"] if r["n_total"] else 0
    body += ("<p class='muted small'><a href='https://www.pgscatalog.org/score/%s/'>%s</a> · %s %s · %s of %s variants "
             "used%s</p>" % (esc(r["id"]), esc(r["id"]), esc(r["author"]), esc(r["year"]), format(r["n_used"], ","),
                             format(r["n_total"], ","), " — <strong>low coverage, read with care</strong>" if cover < 0.6 else ""))
    return head, body


# single-gene cards that belong together with a published score or watched condition (shown as one card)
RELATED = [("Scoliosis gene", "Scoliosis"), ("Childhood asthma gene", "Childhood-onset asthma"), ("Coeliac disease", "Coeliac disease"), ("Type 1 diabetes", "Type 1 diabetes"),
           ("Psoriasis", "Psoriasis"), ("Multiple sclerosis", "Multiple sclerosis"),
           ("Macular degeneration", "Age-related macular degeneration"), ("Parkinson's", "Parkinson's disease"),
           ("FTO", "Body mass index"), ("APOE", "Alzheimer's disease")]


def score_short(name, r):
    """One phrase for 'At a glance'."""
    info = SCORE_INFO.get(name)
    label = esc(SCORE_LABEL.get(name, name))
    if info and info[0] == "risk":
        or_sd = r.get("or_per_sd") if r.get("or_per_sd") and 1 < r.get("or_per_sd") < 6 else info[2]
        odds = info[1] / (1 - info[1]) * or_sd ** max(min(r["z"], 3), -3)
        return "%s: <strong>%s</strong> lifetime chance (typical %s)" % (label, esc(one_in(odds / (1 + odds))), esc(one_in(info[1])))
    if info:
        return "%s: <strong>%s</strong>" % (label, esc("toward " + (info[1] if r["pct"] > 50 else info[2])))
    return "%s: <strong>%s</strong>" % (label, esc(pct_words(r["pct"])))


def related_score(title):
    return next((s for prefix, s in RELATED if title.startswith(prefix)), None)


def trait_card(t, extra=""):
    return ("<div class='card'><div class='ch'><h3>%s</h3>%s</div><div class='big'>%s</div>"
            "<p style='margin:6px 0 0'>%s</p><p class='muted small'>%s%s</p>%s</div>"
            % (esc(t["title"]), badge(t["ev"]), esc(t["result"]), esc(t["text"]), esc(t["variants"]),
               "".join(" · " + esc(f) for f in t["flags"]), extra))


def attached_score(name, r, s, cfg):
    """A published score (or watched-condition detail) shown inside a related single-gene card."""
    out = ["<div class='attach'>"]
    if r and r.get("pct") is not None:
        head, body = score_block(name, r)
        out.append("<div class='muted small'><strong>All variants together (published score)</strong></div>%s%s"
                   % (head.replace("class='big'", "class='big' style='font-size:18px'"), body))
    elif s is not None:
        out.append("<div class='muted small'><strong>No published score yet</strong> for %s.</div>" % esc(name))
    if s is not None:
        out.append("<details class='inner'><summary>Details: variant by variant</summary>%s</details>"
                   % watched_detail(name, s, cfg or {}))
    out.append("</div>")
    return "".join(out)


def watched_detail(name, s, cfg):
    """The catalog variant count for a watched condition (shown folded)."""
    if not s or not s["n_counted"]:
        return "<p class='muted small'>No usable catalog results yet.</p>"
    ew = explorer_words(s["z"]) if s["z"] is not None else "about typical"
    phrase = "about as many of the variants linked to this condition as a typical person" if ew == "about typical" \
        else "%s of the variants linked to this condition than a typical person" % ew
    body = ("<p class='small' style='margin:6px 0 0'>Catalog count: you carry <strong>%s</strong> (%.1f vs %.1f, "
            "%+.1f SD, %d independent signals). This is a rough count, not a risk estimate — the published score, "
            "where there is one, is the better measure.</p>"
            % (esc(phrase), s["observed"], s["expected"], s["z"] or 0, s["n_counted"]))
    if s["n_counted"] < 5:
        body += "<p class='muted small'>Based on very few signals — not enough to say anything yet.</p>"
    if s.get("n_own_sex"):
        body += "<p class='muted small'>%d of these signals come from studies in women only.</p>" % s["n_own_sex"] \
            if cfg.get("sex") == "female" else ""
    if s.get("n_other_pop"):
        body += ("<p class='muted small'>%d signals were found in non-European studies; you're compared using "
                 "European frequencies for them.</p>" % s["n_other_pop"])
    ex = (cfg.get("_excluded") or {}).get(name)
    if ex and cfg.get("ancestry") == "european" and ex[0] >= 20:
        body += ("<p class='muted small'>%s results from non-European studies (mostly %s) left out by your "
                 "ancestry setting.</p>" % (format(ex[0], ","), esc(ex[1])))
    strong = [r for r in s["records"] if r["status"] == "counted" and r.get("or_used") and r["dose"] is not None]
    strong.sort(key=lambda r: -abs(math.log(r["or_used"])))
    if strong:
        rows_ = "".join(
            "<tr><td>%s<br><span class='muted'>%s</span></td><td class='n'>%s of 2</td><td class='n'>%.2f×</td>"
            "<td class='n'>%.2f×</td><td class='n'>%d</td></tr>"
            % (esc((r["gene"] or "—").split(",")[0].split(" - ")[0]), esc(r["rsid"]), fmt_dose(r["dose"]),
               r["or"], r["or_used"], r.get("studies", 1)) for r in strong[:6])
        body += ("<table style='font-size:12px'><tr><th>Gene</th><th>Your copies</th><th>Reported</th>"
                 "<th>Used</th><th>Studies</th></tr>%s</table><p class='muted small'>Strongest signals by effect per "
                 "copy. “Used” is the effect after adjusting for how precisely the study measured it: precise "
                 "estimates keep their size, uncertain ones are pulled toward 1×.</p>" % rows_)
    return body


def pgs_card(name, r, s=None, cfg=None):
    ev = r.get("ev", "modest")
    if r.get("none"):
        head = ("<div class='big' style='font-size:17px'>No published score yet</div><p class='muted small'>"
                "Nothing suitable in the PGS Catalog (European-ancestry, genome-wide); each run checks again.</p>")
        ev = "weak"
    else:
        h, b = score_block(name, r)
        head = h + b
    extra = ""
    if s is not None:
        extra = ("<details class='inner'><summary>Details: variant by variant</summary>%s</details>"
                 % watched_detail(name, s, cfg or {}))
    tag = " <span class='muted small'>(watched)</span>" if r.get("watched") else ""
    return ("<div class='card'><div class='ch'><h3>%s%s</h3>%s</div>%s%s</div>"
            % (esc(SCORE_LABEL.get(name, name)), tag, badge(ev), head, extra))


def ancestry_html(cfg):
    anc, sheet = load_db_json("ancestry.json"), sample_sheet()
    out = ["<h3 class='sub'>Ancestry check %s</h3>" % badge("solid")]
    if not anc:
        out.append("<p class='muted'>Appears once your genome is imputed.</p>")
        return "".join(out)
    out.append("<p class='muted'>Your genome placed among the 2,500 people of 1000 Genomes using principal components, "
               "the standard way of summarising genetic similarity. Nearby dots have similar genomes. 1000 Genomes "
               "doesn't include every population, so the closest available groups are shown.</p>")
    for tag, title in (("world", "Worldwide"), ("europe", "Within Europe")):
        r = anc.get(tag)
        if not r:
            continue
        near = ", ".join("%s (%s)" % (POP_NAMES.get(p, p), p) for d, p in r["nearest"][:3])
        out.append("<p style='margin:14px 0 4px'><strong>%s</strong> — closest groups: %s.</p>" % (title, esc(near)))
        if tag == "world":
            out.append(svg_scatter(r["points"], sheet, r["you"], lambda s: POP_COLORS.get(sheet.get(s, ("", "EUR"))[1], "#999"),
                                   [(SUPER_NAMES[k], v) for k, v in POP_COLORS.items()], title))
        else:
            out.append(svg_scatter(r["points"], sheet, r["you"], lambda s: EUR_COLORS.get(sheet.get(s, ("CEU",))[0], "#999"),
                                   [(POP_NAMES[k], v) for k, v in EUR_COLORS.items()], title))
    out.append("<p class='muted'>Mixed ancestry shows up as a dot between clouds; a dot slightly away from one "
               "cloud toward another usually reflects a small share of ancestry from that direction.</p>")
    return "".join(out)


def maternal_html():
    mat = load_db_json("maternal.json")
    out = ["<h3 class='sub'>Maternal line %s</h3>" % badge("solid")]
    if mat and mat.get("none"):
        out.append("<p>Neither of your files contains mitochondrial DNA: the sequencing file has no mitochondrial reads "
                   "and MyHeritage's genotype file has no mitochondrial positions. So your maternal line can't be read "
                   "from this data. A dedicated full-sequence mtDNA test would give it; this lab can read such a file "
                   "if you ever get one.</p>")
    elif not mat:
        out.append("<p class='muted'>Appears after imputation, if your sequencing file contains mitochondrial DNA.</p>")
    else:
        prefix, note = haplo_note(mat["haplogroup"])
        out.append("<div class='card'><div class='muted'>Your mitochondrial haplogroup</div><div class='big'>%s</div>"
                   "<p style='margin:8px 0 0'>Mitochondrial DNA passes only from mother to child, so this traces your "
                   "mother's mother's line back thousands of years.%s</p><p class='muted small'>Quality %.2f · %d "
                   "differences from the reference · PhyloTree 17.2 via Haplogrep 3</p></div>"
                   % (esc(mat["haplogroup"]), (" Haplogroup " + esc(prefix) + " is " + esc(note) + ".") if note else "",
                      mat["quality"] or 0, mat["n_variants"]))
    return "".join(out)


def meds_html(pgx):
    out = []
    if not pgx or not pgx["rows"]:
        return "<p class='muted'>No data yet.</p>"
    out.append("<p class='muted'>Genes with clinical guidelines (CPIC) on how they change the effect or side effects of "
               "specific medicines. <strong>Don't change any medication based on this</strong>, but it's worth "
               "mentioning to a doctor or pharmacist, especially results other than “normal”. Clinical "
               "pharmacogenetic testing confirms these properly.</p>")
    out.append("<div class='tw'><table><tr><th>Gene</th><th>Your result</th><th>What it affects</th><th>Medicines</th>"
               "<th>Your variants</th></tr>")
    for r in pgx["rows"]:
        normalish = pgx_is_normal(r["result"])
        adv = "" if normalish else pgx_advice(r)
        out.append("<tr><td><strong>%s</strong></td><td>%s%s</td><td>%s</td><td>%s</td><td><code>%s</code>%s</td></tr>"
                   % (esc(r["gene"]), esc(r["result"]) if normalish else "<strong>%s</strong>" % esc(r["result"]),
                      ("<br><span class='small'>%s</span>" % esc(adv)) if adv else "", esc(r["meaning"]), esc(r["drugs"]),
                      esc(r["variants"]), "".join("<br><span class='muted'>%s</span>" % esc(f) for f in r["flags"])))
    out.append("</table></div>")
    ad = antidep_rows(pgx)
    if ad:
        out.append("<h3 class='sub'>Antidepressants & ADHD medicines %s</h3><p class='muted'>What the CPIC guidelines "
                   "(2023 for SSRIs/SNRIs) say for your genotypes. Most people respond normally at standard doses; these "
                   "only matter if one of these drugs is ever prescribed.</p>" % badge("solid"))
        out.append("<div class='tw'><table><tr><th>Medicine</th><th>Gene</th><th>Your result</th><th>Guideline for you</th></tr>")
        for drug, gene, res, rec, flagged in ad:
            out.append("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
                       % (esc(drug), esc(gene), esc(res), ("<strong>%s</strong>" if flagged else "%s") % esc(rec)))
        out.append("</table></div><p class='muted small'>%s CYP2D6 here is partial (deletions and extra copies aren't "
                   "visible).</p>" % esc(NO_GUIDELINE))
    return "".join(out)


def apoe_card(pgx, extra=""):
    ap = (pgx or {}).get("apoe")
    if not ap:
        return ""
    return ("<div class='card'><div class='ch'><h3>APOE type (Alzheimer's)</h3>%s</div><div class='big'>%s</div>"
            "<p style='margin:6px 0 0'>APOE is the gene with the strongest single influence on late-onset Alzheimer's. "
            "%s</p><p class='muted small'>%s%s</p>%s</div>"
            % (badge("solid"), esc(ap["type"]), esc("Yours is " + APOE_TEXT.get(ap["type"], "")), esc(ap["variants"]),
               " · confidence %d%%" % round(ap["conf"] * 100) if ap.get("conf") else "", extra))


def pgx_card(pgx, gene, title, part_text):
    r = next((x for x in (pgx or {}).get("rows", []) if x["gene"] == gene), None)
    if not r:
        return ""
    return ("<div class='card'><div class='ch'><h3>%s</h3>%s</div><div class='big' style='font-size:18px'>%s</div>"
            "<p style='margin:6px 0 0'>%s</p><p class='muted small'>%s</p></div>"
            % (esc(title), badge("solid"), esc(r["result"]), esc(part_text), esc(r["variants"])))


EARLY_ONSET_RE = re.compile(r"congenital|cleft|deaf|hearing loss|hearing impairment|blind|retinitis|leber|amaurosis|"
                            r"infantile|infant|neonatal|newborn|childhood|juvenile|intellectual disability|"
                            r"developmental|muscular dystrophy|cystic fibrosis|phenylketonuria|spinal muscular|"
                            r"albinism|dwarf|dysplasia|encephalopathy|mucopolysacchar|tay-sachs|lethal|fetal|"
                            r"hydrops|severe combined|microcephaly|ataxia-telangiectasia|cataract|glycogen storage|"
                            r"adrenal hyperplasia|usher|bardet|joubert|zellweger|krabbe|canavan|niemann")


def clin_group(r):
    """'you' (could matter for your own health) or 'pass' (what you could pass on / already-known childhood conditions)."""
    pheno = (r.get("pheno") or "").lower()
    copies = int(round(r.get("dose") or 0))
    dominant = "dominant" in pheno
    if EARLY_ONSET_RE.search(pheno):
        return "pass"
    if copies == 1 and not dominant:
        return "pass"
    return "you"


def rare_html(clin, genome):
    out = ["<p class='muted'>Rare variants classified as pathogenic in ClinVar (updated %s), among the ones this data "
           "can measure. Stars (0–4) show how established each classification is. Everything here needs clinical "
           "confirmation.</p>" % (datetime.date.fromtimestamp(os.path.getmtime(P("db", "clinvar.txt.gz"))).isoformat()
                                  if os.path.exists(P("db", "clinvar.txt.gz")) else "—")]
    if genome.n_seq:
        out.append("<p class='muted small'>Imputation only covers variants in the 1000 Genomes panel, so most very rare "
                   "variants can't be seen this way. An empty list doesn't rule them out.</p>")

    def table(rows):
        t = ["<div class='tw'><table><tr><th>Variant</th><th>Gene</th><th>You</th><th>Copies</th><th>Stars</th>"
             "<th>Condition(s)</th><th>What this means</th></tr>"]
        for r in rows:
            t.append("<tr><td><a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/'>%s</a><br><span class='muted'>"
                     "%s · %s</span></td><td>%s</td><td><code>%s</code></td><td class='n'>%s</td><td class='n'>%s</td>"
                     "<td>%s</td><td>%s</td></tr>"
                     % (esc(r["vid"]), esc(r["rsid"] or "variant"), esc(r["sig"]),
                        "sequencing %d%%" % round(r["conf"] * 100) if r["src"] == "seq" else chip_word(),
                        esc(r["gene"]), esc(r["gt"]), fmt_dose(r["dose"]), "★" * r["stars"] or "–",
                        esc(r["pheno"][:200]), esc(clin_meaning(r))))
        t.append("</table></div>")
        return "".join(t)

    reviewed = [r for r in clin if r["stars"] >= 1]
    you = [r for r in reviewed if clin_group(r) == "you"]
    passed = [r for r in reviewed if clin_group(r) == "pass"]
    out.append("<h3 class='sub'>Could matter for your own health</h3>")
    out.append(table(you) if you else "<p>None found.</p>")
    out.append("<h3 class='sub'>What you could pass on</h3><p class='muted'>One copy of a variant for a recessive "
               "condition, or a condition that shows from birth or childhood (which you'd already know about). "
               "Carriers are healthy; this matters for children only if a partner carries a variant in the same "
               "gene.</p>")
    out.append(table(passed) if passed else "<p>None found.</p>")
    unrev = [r for r in clin if r["stars"] == 0]
    if unrev:
        out.append("<details><summary>Unreviewed entries (%d) — single submitter, no evidence given</summary>"
                   "<p class='muted'>Zero-star classifications are often outdated or wrong, so they're kept apart.</p>%s"
                   "</details>" % (len(unrev), table(unrev)))
    return "".join(out)



# =============================================================================
# the ten parts, evidence levels, established single-gene traits, score panel
# =============================================================================

PARTS = [
    ("women", "1. Women's & reproductive health"),
    ("nutrition", "2. Nutrition & food"),
    ("body", "3. Body, bones, heart & lungs"),
    ("mind", "4. Mind & behaviour"),
    ("nervous", "5. Nervous system"),
    ("immune", "6. Hormones & immune system"),
    ("eyes", "7. Eyes, skin & blood"),
    ("meds", "8. Medicines"),
    ("family", "9. Ancestry & family"),
    ("rare", "10. Rare variants & what you could pass on"),
]
NAV_SHORT = {"women": "Women's", "nutrition": "Food", "body": "Body", "mind": "Mind", "nervous": "Nerves",
             "immune": "Immune", "eyes": "Eyes & blood", "meds": "Medicines", "family": "Family", "rare": "Rare"}
EVIDENCE = {
    "solid": ("Science: strong", "The link between these genes and this trait is well established and the effect "
                                  "is clear, so the result tells you something real."),
    "modest": ("Science: moderate", "The link is real and replicated, but genes explain only a small part, so the "
                                    "result only nudges things."),
    "weak": ("Science: early", "Early or inconsistent research. Treat it as a curiosity."),
}
NOT_MEASURABLE = {
    "women": "BRCA1/BRCA2 and other hereditary cancer genes: their harmful variants are mostly rare or small "
             "insertions/deletions (including many founder variants specific to particular populations), which this data can't see. "
             "A clinical panel test is needed if breast or ovarian cancer runs in the family.",
    "nutrition": "Non-coeliac gluten sensitivity, IBS and most food allergies have no reliable genetic test.",
    "body": "Familial hypercholesterolaemia (LDLR, APOB, PCSK9) and inherited heart-muscle or rhythm conditions "
            "are mostly caused by rare variants — a clinical test is the route if they run in the family.",
    "mind": "Large rare deletions and duplications such as 22q11.2, which have big effects on psychiatric risk, "
            "need clinical testing (a microarray).",
    "nervous": "Huntington's disease and other repeat expansions, and most rare epilepsy genes, can't be seen here.",
    "immune": "HLA types are estimated from tag variants, not typed directly; clinical HLA typing is exact.",
    "eyes": "Thalassaemia and sickle cell carrier status (relevant in Mediterranean ancestry) are only partly visible "
            "here — a simple blood count and haemoglobin test is the reliable check.",
    "meds": "CYP2D6 gene deletions and extra copies, and HLA-B*15:02 / HLA-B*58:01 (carbamazepine, allopurinol "
            "reactions) need clinical testing.",
    "family": "Paternal line (Y chromosome) — you don't have one. Maternal line needs mitochondrial DNA, "
              "which your sequencing file doesn't contain.",
    "rare": "GJB2 35delG (a common cause of inherited deafness, carried by about 1 in 30–50 Southern Europeans), "
            "spinal muscular atrophy, fragile X and most cystic fibrosis variants are deletions, repeats or copy "
            "changes this data misses. A clinical carrier screen covers them before a pregnancy.",
}

# GRCh37 forward strand, verified against Ensembl GRCh37; "allele" is the letter whose copies are counted.
TRAIT_SITES = {
    "rs4988235": ("2", 136608646, "A"),     # LCT -13910: lactase persistence
    "rs2187668": ("6", 32605884, "T"),      # HLA-DQ2.5 tag
    "rs7454108": ("6", 32681483, "C"),      # HLA-DQ8 tag
    "rs1800562": ("6", 26093141, "A"),      # HFE C282Y
    "rs1799945": ("6", 26091179, "G"),      # HFE H63D
    "rs713598": ("7", 141673345, "G"),      # TAS2R38 A49P, G = taster (PAV)
    "rs1229984": ("4", 100239319, "T"),     # ADH1B*2 fast
    "rs2472297": ("15", 75027880, "T"),     # CYP1A2/CYP1A1: more coffee, faster clearance
    "rs2282679": ("4", 72608383, "G"),      # GC: lower vitamin D
    "rs601338": ("19", 49206674, "A"),      # FUT2 W143X non-secretor
    "rs1801133": ("1", 11856378, "A"),      # MTHFR C677T
    "rs10455872": ("6", 161010118, "G"),    # LPA
    "rs3798220": ("6", 160961137, "C"),     # LPA I4399M
    "rs9939609": ("16", 53820527, "A"),     # FTO
    "rs1815739": ("11", 66328095, "T"),     # ACTN3 R577X, T = X
    "rs11190870": ("10", 102979207, "T"),   # LBX1 scoliosis risk
    "rs7216389": ("17", 38069949, "T"),     # ORMDL3/GSDMB 17q21: childhood-onset asthma
    "rs28929474": ("14", 94844947, "T"),    # SERPINA1 PI*Z
    "rs17580": ("14", 94847262, "A"),       # SERPINA1 PI*S
    "rs34637584": ("12", 40734202, "A"),    # LRRK2 G2019S
    "rs3135388": ("6", 32413051, "A"),      # HLA-DRB1*15:01 tag (multiple sclerosis)
    "rs4349859": ("6", 31365787, "A"),      # HLA-B27 tag
    "rs12191877": ("6", 31252925, "T"),     # HLA-C*06:02 tag (psoriasis)
    "rs1061170": ("1", 196659237, "C"),     # CFH Y402H, C = risk
    "rs10490924": ("10", 124214448, "T"),   # ARMS2 A69S, T = risk
    "rs12913832": ("15", 28365618, "G"),    # HERC2, G = blue
    "rs1805007": ("16", 89986117, "T"),     # MC1R R151C
    "rs1805008": ("16", 89986144, "T"),     # MC1R R160W
}


def trait_analysis(genome):
    g = pgx_genotypes(genome, TRAIT_SITES)
    out = []

    def cnt(rs):
        return int(round(g[rs]["dose"])) if rs in g else None

    def add(part, ev, title, result, text, rss, note=""):
        low = [r for r in rss if r in g and g[r]["conf"] is not None and g[r]["conf"] < MIN_CONF]
        dis = [r for r in rss if r in g and g[r].get("disagree")]
        flags = []
        if low:
            flags.append("uncertain imputation at " + ", ".join(low))
        if dis:
            flags.append("%s and your imputation disagree at %s" % (chip_word(), ", ".join(dis)))
        if note:
            flags.append(note)
        out.append({"part": part, "ev": ev, "title": title, "result": result, "text": text, "flags": flags,
                    "variants": " · ".join("%s %s" % (r, g[r]["gt"]) for r in rss if r in g) or "no data"})

    # ---- nutrition
    if "rs4988235" in g:
        c = cnt("rs4988235")
        add("nutrition", "solid", "Lactose (milk sugar)",
            ["Likely lactose intolerant as an adult", "Lactase persistent (one copy)", "Lactase persistent (two copies)"][c],
            ["You don't carry the European lactase-persistence version, so your lactase production most likely "
             "dropped after childhood — the original human pattern. Many people like this still handle small "
             "amounts of milk, yoghurt and cheese well.",
             "One copy of the persistence version is enough: you most likely keep digesting milk sugar as an adult.",
             "You most likely keep digesting milk sugar as an adult."][c], ["rs4988235"],
            "Other persistence variants exist in African and Middle Eastern ancestry; this checks the European one.")
    if "rs2187668" in g:
        dq2, dq8 = cnt("rs2187668"), cnt("rs7454108") or 0
        no_dq8 = "rs7454108" not in g
        if dq2 == 0 and dq8 == 0 and no_dq8:
            res, txt = "No DQ2.5 found (DQ8 not measured)", (
                "About 90% of people with coeliac disease carry DQ2.5, so not carrying it makes coeliac disease "
                "much less likely. DQ8 (most of the remaining cases) and the rarer DQ2.2 couldn't be read from "
                "this data, so it can't be fully ruled out — a clinical HLA test can.")
        elif dq2 == 0 and dq8 == 0:
            res, txt = "No DQ2.5 or DQ8 found", ("About 95% of people with coeliac disease carry DQ2 or DQ8, so not "
                                                 "carrying them makes coeliac disease very unlikely. This doesn't cover "
                                                 "the less common DQ2.2 type.")
        else:
            res = " + ".join(([("DQ2.5 (2 copies)" if dq2 == 2 else "DQ2.5")] if dq2 else []) +
                             (["DQ8"] if dq8 else [])) + " carried"
            txt = ("About 1 in 5 Europeans carry DQ2.5 (and 30–40% carry DQ2 or DQ8), but only around 3 in 100 "
                   "carriers develop coeliac disease, so it's a permission, not a prediction. It matters mostly as a "
                   "way to rule coeliac out.")
        add("nutrition", "solid", "Coeliac disease (HLA-DQ2 / DQ8)", res, txt, ["rs2187668", "rs7454108"],
            "Estimated from tag variants that track these HLA types well in Europeans." +
            (" The DQ8 tag isn't in the imputation panel or your file." if no_dq8 else ""))
    if "rs1800562" in g:
        cy, hd = cnt("rs1800562"), cnt("rs1799945") or 0
        if cy == 2:
            res, txt = "C282Y/C282Y", ("The main genotype behind hereditary haemochromatosis (iron overload). Many "
                                       "women are protected by menstrual blood loss until menopause, but a ferritin and "
                                       "transferrin saturation blood test is clearly worth doing.")
        elif cy == 1 and hd >= 1:
            res, txt = "C282Y/H63D", "A mildly raised chance of iron build-up; an iron blood test now and then is reasonable."
        elif cy == 1 or hd:
            res, txt = ("One C282Y copy (carrier)" if cy else "H63D (%d cop%s)" % (hd, "y" if hd == 1 else "ies"),
                        "Carriers almost never develop iron overload; it matters mainly for children if a partner is "
                        "also a carrier." if cy else "H63D on its own rarely causes iron overload.")
        else:
            res, txt = "No haemochromatosis variants", "Typical iron handling as far as these two common variants go."
        add("nutrition", "solid", "Iron overload (HFE)", res, txt, ["rs1800562", "rs1799945"])
    if "rs713598" in g:
        c = cnt("rs713598")
        add("nutrition", "solid", "Bitter taste (TAS2R38)", ["Likely non-taster", "Medium taster", "Likely strong taster"][c],
            "How bitter you find broccoli, Brussels sprouts, rocket, tonic water and grapefruit. This gene explains "
            "most of the difference in tasting one bitter compound (PROP/PTC) — one of the clearest single-gene traits.",
            ["rs713598"], "Read from one of the gene's three defining variants.")
    if "rs1229984" in g:
        c = cnt("rs1229984")
        add("nutrition", "solid", "Alcohol (ADH1B)", ["Typical alcohol processing", "Fast alcohol processing (one copy)",
                                                     "Fast alcohol processing (two copies)"][c],
            ["The usual European version.",
             "Alcohol turns into acetaldehyde faster, which makes drinking more unpleasant (flushing, nausea) and is "
             "linked to drinking less and lower alcohol-dependence risk.",
             "Alcohol turns into acetaldehyde faster, which makes drinking more unpleasant and is linked to drinking "
             "less and lower alcohol-dependence risk."][c], ["rs1229984"])
    if "rs2472297" in g:
        c = cnt("rs2472297")
        add("nutrition", "modest", "Caffeine", ["Typical caffeine processing", "Slightly faster caffeine processing",
                                               "Faster caffeine processing"][c],
            "People with more copies tend to clear caffeine faster and drink a bit more coffee. The effect is real "
            "but small (around a fifth of a cup a day per copy).", ["rs2472297"])
    if "rs2282679" in g:
        c = cnt("rs2282679")
        add("nutrition", "modest", "Vitamin D level (GC)", ["Typical", "Slightly lower vitamin D tendency",
                                                           "Lower vitamin D tendency"][c],
            "Each copy lowers blood vitamin D a little on average. Sunlight and diet matter far more; a blood test "
            "is the real measure.", ["rs2282679"])
    if "rs601338" in g:
        c = cnt("rs601338")
        add("nutrition", "solid", "Secretor status & vitamin B12 (FUT2)",
            ["Secretor", "Secretor (one non-secretor copy)", "Non-secretor"][c],
            ["The usual type (about 80% of Europeans).",
             "One non-functional copy; you're still a secretor.",
             "About 1 in 5 Europeans: you don't put blood-group sugars into saliva and gut mucus. Non-secretors "
             "tend to have higher vitamin B12 levels, are largely resistant to the most common norovirus strains, "
             "and have a different gut microbiome."][c], ["rs601338"])
    if "rs1801133" in g:
        c = cnt("rs1801133")
        add("nutrition", "solid", "Folate processing (MTHFR C677T)", ["Typical", "One 677T copy", "Two 677T copies"][c],
            ["The usual version.", "Very common (about 40% of people); essentially no practical effect.",
             "About 1 in 10 Europeans. Homocysteine runs slightly higher when folate intake is low. Despite "
             "what many websites say, guidelines don't recommend testing for this or taking special "
             "“methylfolate” — the usual folic acid advice before pregnancy is the same."][c],
            ["rs1801133"], "The effect on homocysteine is well proven; the medical importance is small.")

    # ---- body, bones, heart & lungs
    if "rs10455872" in g:
        c = cnt("rs10455872") + (cnt("rs3798220") or 0)
        add("body", "solid", "Lipoprotein(a)", "No common raising variant" if c == 0 else "%d raising cop%s" % (c, "y" if c == 1 else "ies"),
            "Lp(a) is a cholesterol particle set almost entirely by genes, and a high level raises heart-attack and "
            "stroke risk. " + ("These two variants explain a large share of high levels in Europeans; you have "
                               "neither, but other variants also matter." if c == 0 else
                               "Carriers usually have clearly higher Lp(a). European guidelines recommend measuring "
                               "Lp(a) once in everyone's life anyway, and for you a one-off blood test is worth asking for."),
            ["rs10455872", "rs3798220"])
    if "rs28929474" in g or "rs17580" in g:
        z, s_ = cnt("rs28929474") or 0, cnt("rs17580") or 0
        al = (["Z"] * z + ["S"] * s_ + ["M", "M"])[:2]
        typ = "PI*" + "".join(sorted(al, key="MSZ".index))
        txt = {"PI*MM": "The usual, fully working type.",
               "PI*MS": "Common in Europe (roughly 1 in 10 people, more in some regions). One S copy keeps alpha-1 antitrypsin "
                        "levels near normal and isn't linked to disease on its own. Not smoking matters more.",
               "PI*MZ": "Mildly lower alpha-1 antitrypsin; a somewhat higher lung and liver risk with smoking or heavy "
                        "drinking. Worth a blood test of alpha-1 antitrypsin level.",
               "PI*SS": "Mildly lower levels, usually without problems.",
               "PI*SZ": "Lower alpha-1 antitrypsin; a raised lung risk, especially with smoking — worth a blood test.",
               "PI*ZZ": "Alpha-1 antitrypsin deficiency: worth discussing with a doctor (lungs and liver)."}.get(typ, "")
        add("body", "solid", "Alpha-1 antitrypsin (lungs & liver)", typ, txt, ["rs17580", "rs28929474"])
    if "rs9939609" in g:
        c = cnt("rs9939609")
        add("body", "modest", "FTO (body weight)", ["No FTO raising copies", "One FTO raising copy", "Two FTO raising copies"][c],
            "The best-known weight gene: each copy adds about 1.5 kg on average, partly through appetite. Physical "
            "activity reduces its effect by about a third.", ["rs9939609"])
    if "rs1815739" in g:
        c = cnt("rs1815739")
        add("body", "modest", "Muscle type (ACTN3)", ["RR — full sprint protein", "RX — in between",
                                                       "XX — no sprint protein"][c],
            ["You make α-actinin-3, a protein in fast-twitch (explosive) muscle fibres. This type is the most common "
             "one among elite sprinters and power athletes. In everyday life it makes little difference.",
             "You make α-actinin-3 from one copy — the most common type, in between the other two. In everyday "
             "life it makes little difference.",
             "You don't make α-actinin-3, a protein in fast-twitch (explosive) muscle fibres. About 1 in 5 Europeans "
             "are like this, with no health effect. It's rare among elite sprinters and power athletes, and in some "
             "studies slightly more common among elite endurance athletes. For normal exercise it makes little "
             "difference."][c], ["rs1815739"])
    if "rs7216389" in g:
        c = cnt("rs7216389")
        add("body", "modest", "Childhood asthma gene (17q21)", "%d of 2 risk copies" % c,
            "The strongest known genetic region for asthma that starts in childhood (genes ORMDL3/GSDMB): each copy "
            "raises the odds about 1.4×. The risk version is very common — about half of European chromosomes "
            "carry it — so it contributes but doesn't decide. It matters much less for asthma starting in adulthood.",
            ["rs7216389"])
    if "rs11190870" in g:
        c = cnt("rs11190870")
        add("body", "modest", "Scoliosis gene (LBX1)", "%d of 2 risk copies" % c,
            "LBX1 is the most replicated adolescent idiopathic scoliosis gene in all populations, each copy raising the "
            "odds about 1.5×. The risk version is the majority one (about 6 in 10 European chromosomes), so many people "
            "carry two copies — it contributes, but doesn't explain scoliosis alone.", ["rs11190870"])

    # ---- nervous system
    if "rs34637584" in g:
        c = cnt("rs34637584")
        src = g["rs34637584"]["src"]
        add("nervous", "solid", "Parkinson's (LRRK2 G2019S)",
            "Not found" if c == 0 else "G2019S found (%d cop%s) — needs confirmation" % (c, "y" if c == 1 else "ies"),
            ("The most common single genetic cause of Parkinson's. It's rare in most populations but far more common "
             "in some, such as North African Berber and Ashkenazi Jewish communities. Not carrying it is the expected result."
             if c == 0 else
             "This variant raises lifetime Parkinson's risk substantially (estimates vary, roughly 25–40% by old age), "
             "but many carriers never develop it. A rare finding like this must be confirmed by a clinical test "
             "before it means anything — talk to a GP about genetic counselling."),
            ["rs34637584"], "Measured on your chip." if src == "chip" and not chip_is_lowpass() else
            "Estimated from low-pass sequencing, which is unreliable for very rare variants — \u201cnot found\u201d "
            "is likely but not certain; only a clinical test rules it out.")
    if "rs3135388" in g:
        c = cnt("rs3135388")
        add("nervous", "modest", "Multiple sclerosis (HLA-DRB1*15:01)", ["Not carried", "One copy", "Two copies"][c],
            "The strongest genetic factor for MS, raising the odds about 3× per copy. MS is uncommon (about 1 in 1,000 "
            "people in Southern Europe), so even with copies the chance stays low." if c else
            "The strongest genetic factor for MS; you don't carry it.", ["rs3135388"],
            "Estimated from a tag variant.")

    # ---- hormones & immune
    if "rs4349859" in g:
        c = cnt("rs4349859")
        add("immune", "solid", "HLA-B27 (ankylosing spondylitis)", ["Not carried", "Likely carried", "Likely carried (two copies)"][c],
            ("About 8% of Europeans carry HLA-B27. Most never develop ankylosing spondylitis, but most people with it "
             "carry B27 — useful to know if chronic lower-back pain with morning stiffness ever comes up."
             if c else "Not carrying HLA-B27 makes ankylosing spondylitis less likely."), ["rs4349859"],
            "Estimated from a tag variant; a clinical HLA-B27 blood test is exact.")
    if "rs12191877" in g:
        c = cnt("rs12191877")
        add("immune", "solid", "Psoriasis (HLA-C*06:02)", ["Not carried", "Likely carried", "Likely carried (two copies)"][c],
            ("The strongest single genetic factor for psoriasis, raising the odds about 2–4×. Psoriasis affects about "
             "2–3 in 100 people." if c else "The strongest single genetic factor for psoriasis; you don't carry it."),
            ["rs12191877"], "Estimated from a tag variant.")
    if "rs2187668" in g:
        dq2, dq8 = cnt("rs2187668"), cnt("rs7454108") or 0
        add("immune", "solid", "Type 1 diabetes (HLA-DQ2 / DQ8)" if "rs7454108" in g else "Type 1 diabetes (HLA-DQ2)",
            ("DQ2/DQ8 both carried" if dq2 and dq8 else "DQ2.5 carried" if dq2 else "DQ8 carried" if dq8 else
             "Neither carried" if "rs7454108" in g else "No DQ2.5"),
            ("Carrying both DQ2 and DQ8 is the highest-risk HLA combination for type 1 diabetes, though most carriers "
             "never develop it, and adults rarely do." if dq2 and dq8 else
             "These HLA types are the biggest genetic factor for type 1 diabetes; one of them is common and only "
             "mildly raises the odds." if dq2 or dq8 else
             "Without DQ2 or DQ8, type 1 diabetes is considerably less likely." if "rs7454108" in g else
             "No DQ2.5. DQ8, the other main type 1 diabetes HLA type, couldn't be read from this data."),
            ["rs2187668", "rs7454108"], "Estimated from tag variants.")

    # ---- eyes, skin & blood
    if "rs1061170" in g and "rs10490924" in g:
        cf, ar = cnt("rs1061170"), cnt("rs10490924")
        n = cf + ar
        add("eyes", "solid", "Macular degeneration (CFH & ARMS2)", "%d of 4 risk copies (CFH %d, ARMS2 %d)" % (n, cf, ar),
            "These two genes have the largest common effects known for any age-related disease: each copy raises "
            "the odds of age-related macular degeneration about 2–3×. " +
            ("With few or no copies, your risk is below average." if n <= 1 else
             "With several copies your odds are clearly raised; not smoking, a diet rich in leafy greens and fish, "
             "and regular eye checks after 50 help most." if n >= 3 else
             "A fairly typical number of copies."), ["rs1061170", "rs10490924"])
    if "rs12913832" in g:
        c = cnt("rs12913832")
        add("eyes", "solid", "Eye colour (HERC2)", ["Most likely brown eyes", "Likely brown or hazel", "Most likely blue/grey eyes"][c],
            "This one switch explains most blue versus brown eye colour in Europeans. You can check this one against the mirror.",
            ["rs12913832"])
    if "rs1805007" in g or "rs1805008" in g:
        c = (cnt("rs1805007") or 0) + (cnt("rs1805008") or 0)
        add("eyes", "solid", "Red hair & sun sensitivity (MC1R)", ["No major red-hair variants", "One red-hair variant (carrier)",
                                                                  "Two red-hair variants"][min(c, 2)],
            ["No copies of the two strongest red-hair variants.",
             "Carriers usually don't have red hair but tend to freckle and burn more easily, and have a slightly higher "
             "melanoma risk — sunscreen matters a little more.",
             "Two copies usually means red hair, fair skin and freckles, with a higher melanoma risk."][min(c, 2)],
            ["rs1805007", "rs1805008"], "Other, weaker MC1R variants aren't counted.")
    with open(P("db", "traits.json"), "w") as f:
        json.dump(out, f)
    return out


# Published polygenic scores: name -> (part, evidence, search words, preferred PGS ID or None)
PGS_PANEL = {
    "Breast cancer": ("women", "modest", ["breast carcinoma", "breast cancer"], "PGS000004"),
    "Polycystic ovary syndrome": ("women", "weak", ["polycystic ovary"], None),
    "Age at menopause": ("women", "modest", ["menopause"], None),
    "Age at first period": ("women", "modest", ["menarche"], "PGS003508"),
    "Coronary artery disease": ("body", "modest", ["coronary artery disease"], "PGS000018"),
    "Height": ("body", "solid", ["body height"], None),
    "Body mass index": ("body", "modest", ["body mass index"], None),
    "LDL cholesterol": ("body", "modest", ["low density lipoprotein cholesterol", "ldl cholesterol"], None),
    "Schizophrenia": ("mind", "modest", ["schizophrenia"], "PGS000135"),
    "Bipolar disorder": ("mind", "modest", ["bipolar disorder"], None),
    "Depression": ("mind", "modest", ["major depressive disorder", "unipolar depression"], None),
    "Insomnia": ("mind", "modest", ["insomnia"], None),
    "Chronotype (morning person)": ("mind", "modest", ["chronotype", "morningness"], None),
    "Migraine": ("nervous", "modest", ["migraine"], None),
    "Alzheimer's disease": ("nervous", "modest", ["alzheimer"], None),
    "Parkinson's disease": ("nervous", "modest", ["parkinson"], None),
    "Multiple sclerosis": ("nervous", "modest", ["multiple sclerosis"], None),
    "Type 1 diabetes": ("immune", "solid", ["type 1 diabetes"], None),
    "Type 2 diabetes": ("immune", "modest", ["type 2 diabetes"], "PGS000014"),
    "Hypothyroidism": ("immune", "modest", ["hypothyroidism", "autoimmune thyroid"], None),
    "Thyroid hormone (TSH) level": ("immune", "modest", ["thyroid stimulating hormone"], None),
    "Rheumatoid arthritis": ("immune", "modest", ["rheumatoid arthritis"], None),
    "Asthma": ("body", "modest", ["asthma"], None),
    "Childhood-onset asthma": ("body", "modest", ["childhood onset asthma", "childhood asthma"], None),
    "Psoriasis": ("immune", "modest", ["psoriasis"], None),
    "Coeliac disease": ("immune", "solid", ["celiac disease", "coeliac"], None),
    "Age-related macular degeneration": ("eyes", "modest", ["macular degeneration"], None),
    "Myopia (short-sightedness)": ("eyes", "modest", ["myopia", "refractive error"], None),
}
WATCH_PARTS = [("women", ("endometri", "polycystic", "menstru", "pre-eclampsia", "preeclampsia", "infertil")),
               ("body", ("scolio", "spine", "bone", "osteo", "heart", "coronary", "asthma")),
               ("mind", ("adhd", "attention", "autis", "dyslex", "obsess", "ocd", "depress", "anxi", "bipolar",
                         "schizo", "neurodev")),
               ("nervous", ("migraine", "parkinson", "alzheimer", "epilep", "sclerosis")),
               ("immune", ("thyroid", "diabet", "coeliac", "celiac", "lupus", "arthritis", "psoria", "crohn"))]


def watch_part(name):
    n = name.lower()
    for part, words in WATCH_PARTS:
        if any(wd in n for wd in words):
            return part
    return "body"


# CPIC antidepressant / ADHD guidance (CPIC 2023 SSRI/SNRI, 2016 tricyclics, 2019 atomoxetine), by phenotype.
ANTIDEP = [
    ("Citalopram, escitalopram", "CYP2C19", {
        "Poor": "Consider an antidepressant not mainly cleared by CYP2C19; if used, start at half the dose, increase slowly, lower maximum dose.",
        "Intermediate": "Standard starting dose, but slower dose increases and a lower maintenance dose than usual may be enough.",
        "Rapid": "Standard starting dose; if it doesn't help, a higher dose or a switch may be needed.",
        "Ultrarapid": "Consider a different antidepressant not mainly cleared by CYP2C19.",
        "Normal": "Standard dosing."}),
    ("Sertraline", "CYP2C19", {
        "Poor": "Consider a lower starting dose and slower increases (or an alternative).",
        "Intermediate": "Standard starting dose; slower dose increases may help with side effects.",
        "Rapid": "Standard starting dose; if it doesn't help, consider an alternative.",
        "Ultrarapid": "Standard starting dose; if it doesn't help, consider an alternative.",
        "Normal": "Standard dosing."}),
    ("Sertraline", "CYP2B6", {
        "Poor": "Consider a lower starting dose and slower increases.",
        "Intermediate": "Standard starting dose; slower increases if side effects appear.",
        "Normal": "Standard dosing."}),
    ("Paroxetine, fluvoxamine", "CYP2D6 (partial)", {
        "Poor": "Consider a 50% lower starting dose and slower increases, or an alternative.",
        "Intermediate": "Standard starting dose.", "Normal": "Standard dosing."}),
    ("Venlafaxine, vortioxetine", "CYP2D6 (partial)", {
        "Poor": "Venlafaxine: consider an alternative. Vortioxetine: start at half the dose.",
        "Intermediate": "Standard dosing.", "Normal": "Standard dosing."}),
    ("Amitriptyline, nortriptyline (tricyclics)", "CYP2D6 (partial)", {
        "Poor": "Avoid, or use a much lower dose with blood-level monitoring.",
        "Intermediate": "Consider a 25% lower dose, with monitoring.", "Normal": "Standard dosing."}),
    ("Amitriptyline (tricyclics)", "CYP2C19", {
        "Poor": "Consider a 50% lower starting dose.", "Ultrarapid": "Consider an alternative.",
        "Rapid": "Consider an alternative.", "Intermediate": "Standard dosing.", "Normal": "Standard dosing."}),
    ("Atomoxetine (ADHD)", "CYP2D6 (partial)", {
        "Poor": "Start at a lower dose and increase cautiously; levels run much higher.",
        "Intermediate": "Standard start; if side effects or no effect, blood-level guided dosing helps.",
        "Normal": "Standard dosing."}),
]
NO_GUIDELINE = ("Fluoxetine, bupropion, mirtazapine, duloxetine and lisdexamfetamine (Vyvanse/Elvanse) have no gene-based dosing "
                "guideline. Lisdexamfetamine is activated in red blood cells, not by these liver enzymes, so these results "
                "don't affect it.")


def antidep_rows(pgx):
    by_gene = {r["gene"]: r["result"] for r in (pgx or {}).get("rows", [])}
    out = []
    for drug, gene, rec in ANTIDEP:
        res = by_gene.get(gene)
        if not res:
            continue
        key = next((k for k in rec if res.startswith(k) or (" " + k) in res), None)
        if key:
            out.append((drug, gene, res.split(" (")[0], rec[key],
                        not rec[key].startswith(("Standard dosing", "Standard starting dose."))))
    return out


# =============================================================================
# plain-language interpretation
# =============================================================================

TRICKY_GENES = {"CYP21A2", "PMS2", "SMN1", "SMN2", "GBA", "GBA1", "HBA1", "HBA2", "PKD1", "NCF1", "STRC",
                "CYP2D6", "CFC1", "IKBKG", "OTOA", "TNXB"}
GENE_NOTES = {
    "SERPINA1": "The “S” version (rs17580) is common in Europe. One copy on its own rarely "
                "causes problems; it matters mainly together with a “Z” copy, and not smoking protects carriers.",
    "CYP21A2": "This gene sits next to a near-identical copy gene, so DNA tests often misread it. Treat this as "
               "unconfirmed until a clinical test says otherwise.",
    "HFE": "Linked to iron overload (haemochromatosis); one copy rarely matters, two copies are worth a simple "
           "blood test for iron levels.",
}
RECESSIVE_HINTS = ("deficiency", "recessive", "cystic fibrosis", "phenylketonuria", "hyperplasia", "atrophy",
                   "anemia", "anaemia", "thalass", "tay-sachs", "gaucher", "storage disease", "galactos")


def clin_meaning(r):
    """One or two plain sentences on what a ClinVar finding means for you."""
    pheno = (r.get("pheno") or "").lower()
    sig = (r.get("sig") or "").lower()
    copies = int(round(r.get("dose") or 0))
    dominant = "autosomal dominant" in pheno or "dominant" in pheno
    recessive = not dominant and any(h in pheno for h in RECESSIVE_HINTS)
    parts = []
    if recessive:
        parts.append("One copy makes you a carrier: carriers are almost always healthy, and it matters mainly for "
                     "children, if a partner carries a variant in the same gene." if copies == 1 else
                     "Two copies: if confirmed, this could affect you directly — worth raising with a doctor.")
    elif dominant:
        parts.append("Listed as a dominant condition, so one copy could matter — if the finding is confirmed.")
    elif not (EARLY_ONSET_RE.search(pheno) and copies == 2):
        parts.append("How this is inherited isn't stated here; %s." %
                     ("with one copy you are most likely a carrier" if copies == 1 else "two copies is worth asking a doctor about"))
    if EARLY_ONSET_RE.search(pheno) and (copies == 2 or dominant):
        parts.append("This condition shows from birth or childhood, so if you don't have it, the result is most "
                     "likely a misread or the variant has a milder effect than listed.")
    if "low penetrance" in sig:
        parts.append("“Low penetrance” means most people with it never develop the condition.")
    if r.get("stars", 0) <= 1:
        parts.append("The evidence is weak (%s), so it may well be reclassified." %
                     ("one star" if r.get("stars") == 1 else "no stars"))
    if r.get("gene") in GENE_NOTES:
        parts.append(GENE_NOTES[r["gene"]])
    elif r.get("gene") in TRICKY_GENES:
        parts.append("This gene has a near-identical copy elsewhere in the genome, so it's often misread — confirm clinically.")
    return " ".join(parts)


PGX_ADVICE = {
    "CYP2C19": {"Poor": "clopidogrel (a heart drug) may not work; citalopram, escitalopram and sertraline can build "
                        "up — lower doses or alternatives are usual; omeprazole-type stomach drugs work more strongly.",
                "Intermediate": "clopidogrel may work less well; citalopram, escitalopram and sertraline can reach "
                                "higher levels, so a lower dose may be enough; omeprazole-type drugs work slightly better.",
                "Rapid": "some antidepressants (citalopram, escitalopram) and omeprazole-type drugs may be less "
                         "effective at standard doses.",
                "Ultrarapid": "citalopram, escitalopram and omeprazole-type drugs are often less effective at standard doses."},
    "CYP2C9": {"Intermediate": "ibuprofen-type painkillers and warfarin are cleared more slowly; lower doses may be enough.",
               "Poor": "ibuprofen-type painkillers, warfarin and phenytoin are cleared slowly — lower doses are recommended."},
    "VKORC1": {"Increased": "only relevant if you're ever prescribed warfarin: you'd likely need a lower dose.",
               "High": "only relevant if you're ever prescribed warfarin: you'd likely need a clearly lower dose."},
    "SLCO1B1": {"Decreased": "a higher chance of muscle aches on simvastatin; other statins or lower doses are usual alternatives.",
                "Poor": "a clearly higher chance of muscle side effects on simvastatin; a different statin is usually chosen."},
    "ABCG2": {"Decreased": "rosuvastatin reaches higher levels; a lower starting dose is usual.",
              "Poor": "rosuvastatin reaches higher levels; a lower starting dose is usual."},
    "DPYD": {"Intermediate": "IMPORTANT if you ever need 5-FU or capecitabine chemotherapy: doses must be reduced — mention this to the oncologist.",
             "Poor": "IMPORTANT: 5-FU and capecitabine chemotherapy can be dangerous — tell any oncologist."},
    "TPMT": {"Intermediate": "azathioprine and mercaptopurine need lower doses.", "Poor": "azathioprine and mercaptopurine need much lower doses or alternatives."},
    "NUDT15": {"Intermediate": "azathioprine and mercaptopurine need lower doses.", "Poor": "azathioprine and mercaptopurine need much lower doses or alternatives."},
    "UGT1A1": {"Poor": "likely Gilbert's syndrome (harmless mild jaundice, e.g. when fasting or ill); matters for the chemo drug irinotecan.",
               "Intermediate": "mostly harmless; matters mainly for the chemo drug irinotecan."},
    "CYP2D6 (partial)": {"Intermediate": "codeine and tramadol may work less well; paroxetine, venlafaxine and "
                                         "tricyclics such as amitriptyline can reach higher levels — guidelines suggest "
                                         "lower doses or slower increases for some of them.",
                         "Poor": "codeine and tramadol likely won't relieve pain; paroxetine, fluvoxamine, venlafaxine, "
                                 "vortioxetine and tricyclics build up — lower doses or alternatives; tamoxifen works less well."},
    "CYP2B6": {"Intermediate": "sertraline may reach somewhat higher levels; guidelines say a standard start, with slower "
                               "dose increases if needed.",
               "Poor": "sertraline reaches higher levels; guidelines suggest a lower starting dose and slower increases."},
    "F5 / F2 (clotting)": {"Factor": "worth telling a doctor before using estrogen-containing contraception or hormone therapy.",
                           "Prothrombin": "worth telling a doctor before using estrogen-containing contraception or hormone therapy."},
    "G6PD (X chromosome)": {"Deficiency": "avoid primaquine, tafenoquine, rasburicase and dapsone unless G6PD is tested; avoid fava beans.",
                            "Carrier": "usually mild or no effect in women; worth mentioning before primaquine, rasburicase or dapsone.",
                            "Two": "worth mentioning before primaquine, rasburicase, dapsone or nitrofurantoin; avoid fava beans."},
}


def pgx_advice(row):
    for key, text in PGX_ADVICE.get(row["gene"], {}).items():
        if row["result"].startswith(key) or (" " + key) in row["result"]:
            return text
    return ""


def pgx_is_normal(result):
    return any(k in result for k in ("Normal", "Typical", "No ", "Neither", "Non-expresser"))


def pct_words(p):
    return ("well below typical" if p < 10 else "below typical" if p < 30 else "about typical" if p <= 70
            else "above typical" if p <= 90 else "well above typical")


def plain_summary(cfg, summaries, clin, qc, traits, genome):
    """'At a glance': one line per part, findings only."""
    pgs, pgx, mat, anc = load_db_json("pgs.json") or {}, load_db_json("pgx.json"), load_db_json("maternal.json"), \
        load_db_json("ancestry.json")
    tr = load_db_json("traits.json") or []
    out = ["<h2 id='plain'>At a glance</h2>"]
    if qc:
        verdict = "excellent" if qc["agree"] > 0.98 else "good" if qc["agree"] > 0.95 else "questionable — ask Claude"
        if chip_is_lowpass():
            out.append("<p class='muted'>Processing check: <strong>%s</strong> — this lab's imputation and MyHeritage's "
                       "own genotype file agree on %.1f%% of %s shared positions. Both come from the same low-pass "
                       "sequencing, so this confirms the processing, not a second reading of your DNA.</p>"
                       % (verdict, qc["agree"] * 100, format(qc["total"], ",")))
        else:
            out.append("<p class='muted'>Data quality: <strong>%s</strong> — your chip and imputed genome agree on %.1f%% "
                       "of %s shared positions.</p>" % (verdict, qc["agree"] * 100, format(qc["total"], ",")))
    out.append("<div class='card'><ul class='glance'>")
    for part, title in PARTS:
        bits = []
        for t in tr:
            if t["part"] == part and t["ev"] == "solid":
                bits.append("%s: <strong>%s</strong>" % (esc(t["title"]), esc(t["result"])))
        sc = [(n, r) for n, r in pgs.items() if r.get("part") == part and r.get("pct") is not None]
        notable = [(n, r) for n, r in sc if r["pct"] < 10 or r["pct"] > 90]
        for n, r in notable:
            bits.append(score_short(n, r))
        if sc and len(notable) < len(sc):
            bits.append("%s%d published score%s in the typical range" % ("other " if notable else "",
                                                                         len(sc) - len(notable), "s" if len(sc) - len(notable) != 1 else ""))
        if part == "women" and pgx:
            r = next((x for x in pgx["rows"] if x["gene"].startswith("F5")), None)
            if r:
                bits.append("clotting & the pill: <strong>%s</strong>" % esc(r["result"]))
        if part == "nervous" and pgx and pgx.get("apoe"):
            bits.insert(0, "APOE: <strong>%s</strong>" % esc(pgx["apoe"]["type"]))
        if part == "meds" and pgx:
            flagged = [r for r in pgx["rows"] if not pgx_is_normal(r["result"]) and not r["gene"].startswith(("F5", "G6PD"))]
            bits += ["%s: <strong>%s</strong>" % (esc(r["gene"]), esc(r["result"].split(" (")[0])) for r in flagged]
            n_ok = sum(1 for r in pgx["rows"] if pgx_is_normal(r["result"]))
            bits.append("%d other genes normal" % n_ok)
        if part == "family":
            if anc and anc.get("world"):
                w0 = anc["world"]["nearest"][0][1]
                bits.append("closest reference group: <strong>%s</strong>" % esc(POP_NAMES.get(w0, w0)))
            if mat and not mat.get("none"):
                bits.append("maternal haplogroup <strong>%s</strong>" % esc(mat["haplogroup"]))
        if part == "rare":
            rev = [r for r in clin if r["stars"] >= 1]
            you = [r for r in rev if clin_group(r) == "you"]
            passed = [r for r in rev if clin_group(r) == "pass"]
            bits.append("<strong>%d</strong> that could matter for you%s" % (len(you), (" (" + ", ".join(esc(r["gene"]) for r in you) + ")") if you else ""))
            bits.append("<strong>%d</strong> you could pass on%s" % (len(passed), (" (" + ", ".join(esc(r["gene"]) for r in passed) + ")") if passed else ""))
        if bits:
            out.append("<li><a href='#%s'>%s</a> — %s</li>" % (part, esc(title.split(". ", 1)[1]), "; ".join(bits)))
    out.append("</ul></div>")
    return "".join(out)


# =============================================================================
# report
# =============================================================================

CSS = """
:root{--bg:#faf8f5;--fg:#1f1d1a;--muted:#6b665e;--card:#fff;--line:#e6e1d8;--accent:#5b4bb7;--warn:#a4441c;--soft:#f1edf9;--new:#1f7a4d}
@media (prefers-color-scheme:dark){:root{--bg:#161514;--fg:#ece8e1;--muted:#a19a8f;--card:#1f1e1c;--line:#34312d;--accent:#a99cf0;--warn:#f0a37f;--soft:#26233a;--new:#6fd3a0}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Helvetica Neue",sans-serif}
main{max-width:1100px;margin:0 auto;padding:32px 16px 80px}
nav{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--line);padding:10px 0;z-index:2;display:flex;gap:14px;overflow-x:auto;white-space:nowrap;font-size:14px}
nav a{text-decoration:none}
h1{font-size:28px;margin:0 0 4px}h2{font-size:21px;margin:44px 0 8px;scroll-margin-top:60px}h3{margin:0;font-size:17px}
.muted{color:var(--muted)}.note{background:var(--soft);border-radius:10px;padding:14px 16px;margin:16px 0}
.warn{border-left:4px solid var(--warn)}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:12px;margin:16px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.big{font-size:22px;font-weight:600;margin:6px 0 2px}
details{background:var(--card);border:1px solid var(--line);border-radius:12px;margin:10px 0;padding:10px 14px}
summary{cursor:pointer;font-weight:600}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin-top:10px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--muted);font-weight:500;white-space:nowrap}
td.n{text-align:right;white-space:nowrap}.dim td{opacity:.5}
a{color:var(--accent)}code{font-size:12px}
.tag{display:inline-block;font-size:11px;font-weight:600;padding:1px 6px;border-radius:6px;background:var(--new);color:var(--bg);margin-left:6px;vertical-align:1px}
input[type=search]{width:100%;padding:9px 12px;border:1px solid var(--line);border-radius:8px;background:var(--card);color:var(--fg);font-size:14px;margin-top:8px}
ul.papers{padding-left:18px}ul.papers li{margin:8px 0}
.small{font-size:13px}p.small,p.muted.small{margin:6px 0 0}
.ch{display:flex;justify-content:space-between;align-items:flex-start;gap:8px}
.ev{display:inline-block;font-size:11px;font-weight:600;padding:2px 7px;border-radius:999px;white-space:nowrap;cursor:help}
.ev-solid{background:#d9f0e3;color:#1d5c3a}.ev-modest{background:#fbeccb;color:#7a5410}.ev-weak{background:#ece9e4;color:#6b665e}
@media (prefers-color-scheme:dark){.ev-solid{background:#1d3a2b;color:#8fdcb0}.ev-modest{background:#3d3018;color:#f0c674}.ev-weak{background:#2c2a27;color:#a19a8f}}
details.inner{margin:10px 0 0;padding:6px 10px;background:transparent}
details.inner summary{font-weight:500;font-size:13px}
h3.sub{margin:22px 0 6px;display:flex;gap:8px;align-items:center}
ul.glance{margin:0;padding-left:18px}ul.glance li{margin:7px 0}
.nm{border-top:1px solid var(--line);padding-top:8px;margin-top:14px}
.attach{border-top:1px solid var(--line);margin-top:10px;padding-top:8px}
"""


def esc(x):
    return html.escape(str(x if x is not None else ""))


def fmt_p(p):
    return ("%.0e" % p).replace("e-0", "e-")


def fmt_dose(d):
    return "%d" % round(d) if abs(d - round(d)) < 0.05 else "%.1f" % d


def explorer_words(z):
    return ("about typical" if abs(z) < 1 else ("somewhat more" if z < 2 else "clearly more") if z > 0
            else ("somewhat fewer" if z > -2 else "clearly fewer"))


def describe_z(z):
    if z is None:
        return "not enough data"
    if abs(z) < 1:
        return "about average"
    word = "above" if z > 0 else "below"
    return ("somewhat " if abs(z) < 2 else "well ") + word + " average"


def odds_text(odds):
    return "%.1f×" % odds if odds < 9.95 else "%.0f×" % odds


def plain_odds(odds, subject="your known variants"):
    if odds is None:
        return ""
    if odds < 0.8:
        return "Together, %s give you somewhat lower odds than a typical person." % subject
    if odds <= 1.25:
        return "Together, %s put you at about typical odds — no real difference." % subject
    if odds <= 2:
        return ("Together, %s give you somewhat higher odds than a typical person. For a condition that "
                "affects 1 in 100 people, %s means about %.1f in 100." % (subject, odds_text(odds), odds))
    return ("Together, %s give you notably higher odds than typical (%s), but this is still a modest "
            "effect next to family history, and worth checking against the full genome." % (subject, odds_text(odds)))


def gwas_table(name, s, prev_keys):
    out = ["<details><summary>%s — %d variants</summary><div class='tw'><table>" % (esc(name), s["n_total"]),
           "<tr><th>Variant</th><th>Gene</th><th>You</th><th>Risk allele</th><th>Copies</th>"
           "<th>Freq.</th><th>OR/β</th><th>p</th><th>Study</th><th>Status</th></tr>"]
    for r in s["records"]:
        dim = "" if r["status"] == "counted" else " class='dim'"
        dose = fmt_dose(r["dose"]) if r["dose"] is not None else ""
        conf = " <span class='muted'>(%d%%)</span>" % round(r["conf"] * 100) if r.get("conf") else ""
        src = {"seq": "", "chip": " <span class='muted'>%s</span>" % chip_word()}.get(r.get("src"), "")
        newtag = "<span class='tag'>new</span>" if prev_keys is not None and r["key"] not in prev_keys else ""
        study = ("<a href='https://pubmed.ncbi.nlm.nih.gov/%s/'>%s %s</a>%s<br><span class='muted'>%s</span>"
                 % (esc(r["pmid"]), esc(r["author"]), esc(r["year"]), newtag, esc(r["trait"])))
        if r.get("own_sex"):
            study += "<br><span class='muted'>study in your sex only</span>"
        if r.get("other_sex"):
            study += "<br><span class='muted'>study in the other sex only</span>"
        if r.get("other_pop"):
            study += "<br><span class='muted'>non-European study; compared using European frequencies</span>"
        if r.get("studies", 1) > 1:
            study += "<br><span class='muted'>+%d more studies</span>" % (r["studies"] - 1)
        note = "<br><span class='muted'>%s</span>" % esc(r["note"]) if r.get("note") else ""
        out.append(
            "<tr%s><td><a href='https://www.ncbi.nlm.nih.gov/snp/%s'>%s</a></td><td>%s</td>"
            "<td><code>%s</code>%s%s</td><td><code>%s</code></td><td class='n'>%s</td><td class='n'>%s</td>"
            "<td class='n'>%s</td><td class='n'>%s</td><td>%s</td><td>%s%s</td></tr>"
            % (dim, esc(r["rsid"]), esc(r["rsid"]), esc(r["gene"]), esc(r.get("gt", "—")), conf, src,
               esc(r["risk"]), dose, "%.2f" % r["raf"] if r["raf"] else "", esc(r["effect"]),
               fmt_p(r["p"]), study, esc(r["status"]), note))
    out.append("</table></div></details>")
    return "".join(out)


def cmd_report(a, cfg=None):
    cfg = cfg or load_config()
    genome = Genome()
    if not (genome.n_seq or genome.n_chip):
        die("No genotypes yet. Run  python3 genome_lab.py chip YOUR_FILE  and/or  impute YOUR_FILE.cram")
    os.makedirs(P("reports"), exist_ok=True)
    os.makedirs(P("history"), exist_ok=True)
    prev_path = P("history", "latest.json")
    prev = json.load(open(prev_path)) if os.path.exists(prev_path) else None

    papers = json.load(open(P("db", "papers.json")))["papers"] if os.path.exists(P("db", "papers.json")) else {}
    paper_rs = {x for p in papers.values() for x in p["rsids"]}

    lift = get_liftover(cfg)
    say("Reading the GWAS Catalog…")
    watched, explorer, gwas_rows, gwas_rs = read_gwas(cfg, lift, paper_rs)
    say("Reading ClinVar…")
    clin_rows, clin_rs = clinvar_rows(paper_rs)

    # where the variants mentioned in papers sit on your genome (GRCh37)
    paper_pos = dict(clin_rs)
    paper_pos.update(gwas_rs)
    paper_pos.update({rs: (k, False) for rs, k in genome.chip_keys(paper_rs).items()})

    need_keys = {r["k"] for rows in watched.values() for r in rows} | {e[7] for e in explorer}
    need_keys |= {v[0] for v in paper_pos.values()}
    say("Looking up your genotypes at %s positions…" % format(len(need_keys), ","))
    genos = genome.get_many(need_keys)

    if refpop_ready() and cfg.get("ancestry") == "european":
        say("Using European frequencies from the reference people…")
        freqs = eur_freqs({r["k"] for rows in watched.values() for r in rows})
        for rows in watched.values():
            use_eur_freqs(rows, freqs)

    summaries = {}
    keys_now = {}
    for name, rows in watched.items():
        for r in rows:
            r["key"] = "%s|%s" % (r["acc"] or r["pmid"], r["rsid"])
        keys_now[name] = sorted({r["key"] for r in rows})
        summaries[name] = tally(rows, genos)

    say("Exploring all traits…")
    traits = explore_traits(explorer, genos, cfg["explorer_min_signals"])

    say("Checking ClinVar…")
    clin = clinvar_hits(clin_rows, genome)
    pgx_analysis(genome)      # quick; keeps medication and trait wording current with this version
    trait_analysis(genome)
    qc = chip_vs_seq_qc(genome)

    state = {"date": today(), "gwas": keys_now,
             "clinvar": {"%s|%s" % (r["rsid"] or r["k"], r["alt"]): r["sig"] for r in clin},
             "papers": sorted(papers)}

    # ---------------------------------------------------------------- what's new
    new = {"gwas": {}, "clin_new": [], "clin_changed": [], "clin_gone": [], "papers": []}
    if prev:
        for name, keys in keys_now.items():
            old = set(prev.get("gwas", {}).get(name, []))
            if name in prev.get("gwas", {}):
                fresh = set(keys) - old
                if fresh:
                    new["gwas"][name] = fresh
        oc = prev.get("clinvar", {})
        for k, sig in state["clinvar"].items():
            if k not in oc:
                new["clin_new"].append(k)
            elif oc[k] != sig:
                new["clin_changed"].append((k, oc[k], sig))
        new["clin_gone"] = [k for k in oc if k not in state["clinvar"]]
        seen = set(prev.get("papers", []))
        new["papers"] = [pid for pid in papers if pid not in seen]

    # ---------------------------------------------------------------- HTML
    con = genome.con
    src_bits = []
    if genome.n_seq:
        src_bits.append("whole-genome sequencing, imputed (%s variants, %s)"
                        % (format(genome.n_seq, ","), meta_get(con, "seq_date", "")))
    if genome.n_chip:
        src_bits.append(("MyHeritage genotype file (%s genotypes, also from the sequencing)" if chip_is_lowpass()
                         else "chip file (%s genotypes)") % format(genome.n_chip, ","))
    db_date = lambda p: datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat() if os.path.exists(p) else "—"

    o = []
    w = o.append
    w("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
      "content='width=device-width,initial-scale=1'><title>My genome lab</title><style>%s</style></head><body><main>" % CSS)
    w("<h1>My genome lab</h1><p class='muted'>%s · report %s · made locally, nothing uploaded</p>"
      % (esc(" + ".join(src_bits)), today()))
    w("<nav><a href='#plain'>At a glance</a>%s<a href='#papers'>Papers</a><a href='#explorer'>Trait explorer</a>"
      "<a href='#quality'>Data quality</a><a href='#read'>How to read</a></nav>"
      % "".join("<a href='#%s'>%s</a>" % (k, esc(NAV_SHORT[k])) for k, t in PARTS))
    w("<div class='note warn'><strong>Read this first.</strong> This compares your DNA with published research. "
      "It is not a diagnosis and can't rule anything in or out. Common-variant results shift risk by small "
      "amounts, usually far less informative than symptoms and family history. Imputed genotypes are "
      "statistical estimates, and rare 'pathogenic' findings need clinical confirmation before they mean anything.</div>")
    w("<div class='note'><strong>How to read the cards.</strong> Each card has a badge for <em>how established the "
      "science is</em> — it says nothing about whether your result is good or bad, or about your data quality:"
      "<ul class='glance'>%s</ul>Polygenic scores are shown as your <strong>estimated lifetime chance</strong> next to "
      "the typical chance. The percentile behind it is a ranking: “83rd percentile” means your genetic score is "
      "higher than 83 out of 100 European-ancestry people — not an 83%% chance of anything.</div>"
      % "".join("<li>%s %s</li>" % (badge(k), esc(v[1])) for k, v in EVIDENCE.items()))
    w(plain_summary(cfg, summaries, clin, qc, traits, genome))

    # what's new (folded unless something changed)
    items = []
    if prev:
        for name, fresh in new["gwas"].items():
            items.append("<li><strong>%s</strong>: %d new study result(s)</li>" % (esc(name), len(fresh)))
        for k in new["clin_new"]:
            items.append("<li>ClinVar: <strong>%s</strong> — newly classified as pathogenic, or newly measured</li>"
                         % esc(k.split("|")[0]))
        for k, a_, b_ in new["clin_changed"]:
            items.append("<li>ClinVar: %s reclassified from <em>%s</em> to <em>%s</em></li>"
                         % (esc(k.split("|")[0]), esc(a_), esc(b_)))
        for k in new["clin_gone"]:
            items.append("<li>ClinVar: %s no longer among your findings (reclassified or removed)</li>" % esc(k.split("|")[0]))
        if new["papers"]:
            items.append("<li>%d new paper(s) on your watched topics — see Papers</li>" % len(new["papers"]))
    w("<details id='new'%s><summary>What's new%s</summary>%s</details>"
      % (" open" if items else "", (" since %s" % esc(prev.get("date"))) if prev else "",
         ("<ul>%s</ul>" % "".join(items)) if items else
         "<p>%s</p>" % ("Nothing new since then." if prev else
                        "This is your first report. From the next run on, new studies, changed classifications "
                        "and new papers since the previous report show up here.")))

    pgs_now = load_db_json("pgs.json") or {}
    pgx = load_db_json("pgx.json")
    tr = load_db_json("traits.json") or []
    ev_rank = {"solid": 0, "modest": 1, "weak": 2}
    intros = {
        "women": "Conditions and scores specific to women's and reproductive health.",
        "nutrition": "How your body handles particular foods and nutrients — the most solid single-gene traits.",
        "body": "Heart, metabolism, build, bones and lungs.",
        "mind": "Psychiatric and behavioural traits. These are all polygenic: thousands of variants with tiny "
                "effects, so even the best scores only nudge the odds. Medication guidance for antidepressants is "
                "in part 8.",
        "nervous": "Brain and nerve conditions.",
        "immune": "Autoimmune conditions, diabetes and thyroid. HLA types (your immune system's ID tags) are "
                  "among the strongest common genetic effects known.",
        "eyes": "Eyes, pigmentation and blood.",
        "meds": "Medication response (pharmacogenomics) — the part of genetics most directly usable today.",
        "family": "Where your DNA comes from.",
        "rare": "Rare variants with clinical classifications, split into what could affect you and what you "
                "could pass on to children.",
    }
    merged = {related_score(t["title"]) for t in tr} - {None}
    if (pgx or {}).get("apoe"):
        merged.add("Alzheimer's disease")

    def extra_for(rel):
        if not rel or not (rel in pgs_now or rel in cfg["watch"]):
            return ""
        return attached_score(rel, pgs_now.get(rel), summaries.get(rel) if rel in cfg["watch"] else None, cfg)

    for part, title in PARTS:
        w("<h2 id='%s'>%s</h2><p class='muted'>%s</p>" % (part, esc(title), esc(intros[part])))
        cards = []
        for t in tr:
            if t["part"] == part:
                cards.append((ev_rank[t["ev"]], 0, trait_card(t, extra_for(related_score(t["title"])))))
        if part == "women":
            c = pgx_card(pgx, "F5 / F2 (clotting)", "Clotting & estrogen (pill, pregnancy)",
                         "Factor V Leiden and prothrombin G20210A raise the risk of blood clots, especially with "
                         "estrogen-containing contraception, hormone therapy or pregnancy.")
            if c:
                cards.append((0, 0, c))
        if part == "nervous":
            c = apoe_card(pgx, extra_for("Alzheimer's disease"))
            if c:
                cards.append((0, 0, c))
        if part == "eyes":
            c = pgx_card(pgx, "G6PD (X chromosome)", "G6PD (red blood cells)",
                         "Deficiency can cause anaemia with certain drugs or fava beans; women with one copy usually "
                         "have mild or no effects.")
            if c:
                cards.append((0, 0, c))
        for name, r in pgs_now.items():
            if name in merged:
                continue
            if r.get("part", watch_part(name) if name in cfg["watch"] else None) == part:
                s = summaries.get(name) if name in cfg["watch"] else None
                cards.append((ev_rank.get(r.get("ev", "modest"), 1) + (0.5 if r.get("none") else 0), 1,
                              pgs_card(name, dict(r, watched=name in cfg["watch"]), s, cfg)))
        for name in cfg["watch"]:
            if watch_part(name) == part and name not in pgs_now and name not in merged:
                cards.append((1.5, 1, pgs_card(name, {"none": True, "watched": True}, summaries.get(name), cfg)))
        if part == "family":
            w(ancestry_html(cfg))
            w(maternal_html())
        elif part == "meds":
            w(meds_html(pgx))
        elif part == "rare":
            w(rare_html(clin, genome))
        if cards:
            cards.sort(key=lambda x: (x[0], x[1]))
            w("<div class='cards'>%s</div>" % "".join(c for _, _, c in cards))
        mine = [n for n in cfg["watch"] if watch_part(n) == part and summaries.get(n, {}).get("n_total")]
        if mine:
            w("<details><summary>Details: every catalog variant for %s</summary><p class='muted small'>"
              "Genome-wide significant results (p ≤ %g) from the GWAS Catalog. Greyed rows aren't counted "
              "(same region as a stronger hit, or no data).</p>" % (esc(", ".join(mine)), cfg["pmax"]))
            for name in mine:
                prev_keys = set(prev["gwas"][name]) if prev and name in prev.get("gwas", {}) else None
                w(gwas_table(name, summaries[name], prev_keys))
            w("</details>")
        w("<p class='muted small nm'><strong>Not measurable here:</strong> %s</p>" % esc(NOT_MEASURABLE[part]))


    # papers
    w("<h2 id='papers'>Papers</h2><p class='muted'>Genetics papers on your watched topics from Europe PMC "
      "(peer-reviewed and preprints). Variants mentioned in the abstract are checked against your genome.</p>")
    newset = set(new["papers"])
    for name in cfg["watch"]:
        plist = [(pid, p) for pid, p in papers.items() if name in p["topics"]]
        if not plist:
            continue
        plist.sort(key=lambda x: x[1]["date"], reverse=True)
        n_new = sum(1 for pid, _ in plist if pid in newset)
        w("<details%s><summary>%s — %d papers%s</summary><ul class='papers'>"
          % (" open" if n_new else "", esc(name), len(plist), " (%d new)" % n_new if n_new else ""))
        for pid, p in plist[:150]:
            mine = []
            for rs in p["rsids"]:
                pos = paper_pos.get(rs)
                txt = rec_text(genos.get(pos[0])) if pos else None
                mine.append("%s: you're <code>%s</code>" % (rs, txt) if txt else "%s: no data" % rs)
            w("<li><a href='%s'>%s</a>%s%s<br><span class='muted'>%s · %s · %s</span>%s</li>"
              % (esc(p["url"]), esc(p["title"]), "<span class='tag'>new</span>" if pid in newset else "",
                 " <span class='muted'>(preprint, not peer-reviewed)</span>" if p.get("preprint") else "",
                 esc((p["authors"] or "")[:120]), esc(p["journal"]), esc(p["date"]),
                 "<br>" + " · ".join(mine) if mine else ""))
        w("</ul></details>")
    if not papers:
        w("<p class='muted'>No papers fetched yet — run <code>python3 genome_lab.py run</code>.</p>")

    # explorer
    conditions, measurements, suspicious, appearance = traits
    w("<h2 id='explorer'>Trait explorer</h2>")
    w("<p class='muted'>Every condition in the catalog with at least %d independent signals in your data. It "
      "only shows whether you carry <strong>more or fewer</strong> of the variants studies linked to each "
      "condition than a typical person, sorted with the most “more” at the top. There are deliberately no risk "
      "numbers here: adding up hundreds of overlapping catalog results gives odds that look dramatic but aren't "
      "reliable. <strong>Read this as a list of curiosities, not findings</strong> — with hundreds of conditions, "
      "some land high or low by chance. If one interests you, ask for its published score, or follow it with "
      "<code>watch</code>.</p>" % cfg["explorer_min_signals"])
    if not genome.n_seq:
        w("<div class='note'>This is based on your chip file only, which covers a small share of the studied "
          "variants. Expect it to change once your full genome is imputed.</div>")
    w("<input type='search' id='q' placeholder='Filter, e.g. migraine, asthma, iron…'>")
    w("<div class='tw'><table class='ftable'><tr><th>Condition</th><th>You vs typical</th>"
      "<th>Signals</th><th>Linked variants: you / typical</th></tr>")
    for x in conditions:
        w("<tr><td>%s</td><td>%s <span class='muted'>(%+.1f SD)</span></td><td class='n'>%d</td>"
          "<td class='n'>%.1f / %.1f</td></tr>"
          % (esc(x["trait"]), esc(explorer_words(x["z"])), x["z"], x["n"], x["obs"], x["exp"]))
    w("</table></div>")
    w("<details><summary>Measurements (%d) — blood values, body measures and similar</summary>"
      "<p class='muted'>These are traits measured on a scale rather than conditions you have or don't have, "
      "so there are no odds. The difference shows whether you carry more (+) or fewer (−) of the variants "
      "that raise the value.</p><div class='tw'><table class='ftable'><tr><th>Trait</th><th>Difference</th>"
      "<th>Signals</th><th>You / typical</th></tr>" % len(measurements))
    for x in measurements:
        w("<tr><td>%s</td><td class='n'>%+.1f SD</td><td class='n'>%d</td><td class='n'>%.1f / %.1f</td></tr>"
          % (esc(x["trait"]), x["z"], x["n"], x["obs"], x["exp"]))
    w("</table></div></details>")
    w("<details><summary>Appearance (%d) — hair, eye and skin colour and similar</summary>"
      "<p class='muted'>The difference shows whether you carry more (+) or fewer (−) of the versions studies "
      "flagged. For colour traits these are usually the versions linked to <em>lighter</em> or <em>redder</em> "
      "shades, but studies compare different things, so treat this loosely.</p><div class='tw'>"
      "<table class='ftable'><tr><th>Trait</th><th>Difference</th><th>Signals</th><th>You / typical</th></tr>"
      % len(appearance))
    for x in appearance:
        w("<tr><td>%s</td><td class='n'>%+.1f SD</td><td class='n'>%d</td><td class='n'>%.1f / %.1f</td></tr>"
          % (esc(x["trait"]), x["z"], x["n"], x["obs"], x["exp"]))
    w("</table></div></details>")
    w("<details><summary>Probably catalog errors (%d)</summary><p class='muted'>Results so extreme they're "
      "almost certainly caused by inconsistent data in the catalog (for example the wrong letter listed as "
      "the risk version), not by your genome. Shown only for completeness.</p><div class='tw'>"
      "<table class='ftable'><tr><th>Trait</th><th>Difference</th><th>Signals</th></tr>" % len(suspicious))
    for x in suspicious:
        w("<tr class='dim'><td>%s</td><td class='n'>%+.1f SD</td><td class='n'>%d</td></tr>" % (esc(x["trait"]), x["z"], x["n"]))
    w("</table></div></details>")
    w("<script>document.getElementById('q').addEventListener('input',function(e){var q=e.target.value.toLowerCase();"
      "document.querySelectorAll('.ftable tr').forEach(function(r){if(r.cells[0].tagName==='TD')r.style.display="
      "r.cells[0].textContent.toLowerCase().indexOf(q)>-1?'':'none'})})</script>")

    # quality
    w("<h2 id='quality'>Data quality</h2><ul>")
    if genome.n_seq:
        w("<li>Sequencing: %s variants imputed with the 1000 Genomes panel (GRCh37). The percentage next to a "
          "genotype is the imputation's confidence.</li>" % format(genome.n_seq, ","))
    if genome.n_chip:
        if chip_is_lowpass():
            w("<li>MyHeritage genotype file: %s genotypes (%s). Its header says it was made from low-pass whole-genome "
              "sequencing, so it's MyHeritage's own estimate from the same reads, not a separate chip.</li>"
              % (format(genome.n_chip, ","), esc(meta_get(con, "chip_file", ""))))
        else:
            w("<li>Chip file: %s genotypes (%s).</li>" % (format(genome.n_chip, ","), esc(meta_get(con, "chip_file", ""))))
    if qc:
        if chip_is_lowpass():
            w("<li><strong>Processing check:</strong> at %s positions in both, this lab's imputation and MyHeritage's "
              "file agree %.1f%% of the time (%.1f%% where you carry at least one non-reference allele, the harder "
              "case). Two different pipelines reaching the same genotypes from the same reads is reassuring, but it "
              "isn't independent confirmation. Accuracy itself comes from the held-out test: about 98–99%% for common "
              "variants, lower for rare ones.</li>"
              % (format(qc["total"], ","), qc["agree"] * 100, (qc["nonref_agree"] or 0) * 100))
        else:
            w("<li><strong>Cross-check:</strong> at %s positions measured by both, your chip and your imputed "
              "sequencing agree %.1f%% of the time (%.1f%% where you carry at least one non-reference allele, "
              "the harder case). Values above ~95%% mean the imputation worked well.</li>"
              % (format(qc["total"], ","), qc["agree"] * 100, (qc["nonref_agree"] or 0) * 100))
    w("<li>GWAS Catalog downloaded %s · ClinVar %s</li></ul>" % (db_date(P("db", "gwas.zip")), db_date(P("db", "clinvar.txt.gz"))))

    w("<h2 id='read'>How to read this</h2><ul>"
      "<li><strong>Copies</strong>: how many copies of the risk allele you carry (0, 1 or 2; decimals are "
      "imputation uncertainty).</li>"
      "<li><strong>Same region</strong>: variants close together usually tag one signal, so only the strongest is counted.</li>"
      "<li><strong>Odds</strong>: 1.0× is typical; 1.5× means 50% higher odds than average. For a condition "
      "that affects 1 in 100 people, 1.5× means about 1.5 in 100 — still unlikely.</li>"
      "<li><strong>SD</strong>: how unusual your variant count is compared with someone typical; about 2 in 3 "
      "people fall within ±1 SD.</li>"
      "<li>Most studies are in people of European ancestry; frequencies and effects can differ for other ancestries.</li>"
      "<li>Look up any variant on your computer with <code>python3 genome_lab.py lookup rsNUMBER</code>.</li>"
      "</ul></main></body></html>")

    out = P("reports", "report-%s.html" % today())
    with open(out, "w", encoding="utf-8") as f:
        f.write("".join(o))
    shutil.copy(out, P("reports", "latest.html"))
    with open(P("history", "state-%s.json" % today()), "w") as f:
        json.dump(state, f)
    shutil.copy(P("history", "state-%s.json" % today()), prev_path)
    say("\nReport: %s" % out)
    if not getattr(a, "no_open", False) and sys.platform == "darwin":
        subprocess.run(["open", out])


# =============================================================================
# lookup / watch / status / run
# =============================================================================

def cmd_lookup(a):
    cfg = load_config()
    genome = Genome()
    ids = [x.lower() for x in a.rsids]
    say("Searching the databases (about a minute)…")
    pos = {}
    clin = {}
    gw_rows = {r: [] for r in ids}
    if os.path.exists(P("db", "clinvar.txt.gz")):
        rows, cpos = clinvar_rows(ids)
        pos.update(cpos)
        for r in rows:
            if r["rsid"] in ids:
                clin.setdefault(r["rsid"], []).append(r)
    if os.path.exists(P("db", "gwas.zip")):
        lift = get_liftover(cfg)
        with open_gwas() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                snps = re.findall(r"rs\d+", (row.get("SNPS") or "").lower())
                for r in ids:
                    if r in snps:
                        gw_rows[r].append(row)
                        if r not in pos and len(snps) == 1:
                            l = lift.lift((row.get("CHR_ID") or "").strip(), to_float(row.get("CHR_POS")))
                            if l:
                                pos[r] = (gkey(l[0], l[1]), l[2])
    pos.update({r: (k, False) for r, k in genome.chip_keys(ids).items()})
    genos = genome.get_many([v[0] for v in pos.values()])
    for r in ids:
        say("\n" + "=" * 60 + "\n" + r)
        k, flip = pos.get(r, (None, False))
        g = genos.get(k, {}) if k else {}
        if k:
            say("  Position (GRCh37): chr%d:%d" % (k // 1_000_000_000, k % 1_000_000_000))
        if g.get("seq"):
            ref, alt, gt, ds, conf = g["seq"]
            say("  Sequencing: %s   (reference %s, alternative %s, confidence %d%%)"
                % (seq_letters(ref, alt, gt), ref, alt, round(conf * 100)))
        if g.get("chip"):
            say("  %-11s %s/%s" % ("MyHeritage:" if chip_is_lowpass() else "Chip:", g["chip"][0], g["chip"][1]))
        if not g:
            say("  Not in your data." if k else "  Position unknown: not in the GWAS Catalog, ClinVar or your chip file.")
        rows = sorted(gw_rows.get(r, []), key=lambda x: to_float(x.get("P-VALUE")) or 1)
        if rows:
            say("  GWAS Catalog (%d results, strongest first; OR below 1 means the listed allele is protective):" % len(rows))
            for x in rows[:25]:
                ra = (x.get("STRONGEST SNP-RISK ALLELE") or "").split("-")[-1].strip()
                res = allele_dose(g, ra, flip) if ra in COMP else None
                you = "  you carry %s" % fmt_dose(res[0]) if res else ""
                say("   - %s | risk allele %s%s | p=%s | OR/β %s | %s %s (PMID %s)"
                    % (x.get("DISEASE/TRAIT", ""), ra, you, x.get("P-VALUE", ""), x.get("OR or BETA", ""),
                       x.get("FIRST AUTHOR", ""), (x.get("DATE") or "")[:4], x.get("PUBMEDID", "")))
        for c in clin.get(r, []):
            say("  ClinVar: %s — %s (%s>%s) — %s" % (c["sig"], c["gene"], c["ref"], c["alt"], c["pheno"][:150]))
        say("  More: https://www.ncbi.nlm.nih.gov/snp/%s" % r)


def cmd_watch(a):
    cfg = load_config()
    kws = [k.lower() for k in (a.keywords or [a.name])]
    q = a.papers or " OR ".join('"%s"' % k if " " in k else k for k in kws)
    cfg["watch"][a.name] = {"gwas": kws, "papers": q}
    save_config(cfg)
    say("Now watching %s (catalog keywords: %s; papers: %s). It appears in your next  run."
        % (a.name, ", ".join(kws), q))


def cmd_unwatch(a):
    cfg = load_config()
    if cfg["watch"].pop(a.name, None) is None:
        die("Not watching %s. Watched: %s" % (a.name, ", ".join(cfg["watch"])))
    save_config(cfg)
    say("Stopped watching " + a.name)


def cmd_sex(a):
    cfg = load_config()
    cfg["sex"] = a.sex
    save_config(cfg)
    extra = " Chromosome X will be imputed too." if a.sex == "female" else ""
    say("Saved: %s. Studies in your sex are preferred where they exist, then mixed ones; studies done only "
        "in the other sex are used last, and left out for traits that differ between the sexes (sex hormones, "
        "fat distribution, hair loss).%s" % (a.sex, extra))
    if a.sex == "female" and os.path.exists(P("my_genome.imputed.bcf")) \
            and not os.path.exists(P("work", "chrX.bcf.ok")):
        only_x = all(os.path.exists(P("work", "chr%s.bcf.ok" % c)) for c in cfg["chromosomes"])
        say("Your genome is already imputed; run the impute command again to add chromosome X%s."
            % (" (only X is processed — about 30-60 minutes)" if only_x else ""))


def cmd_ancestry(a):
    cfg = load_config()
    cfg["ancestry"] = None if a.ancestry == "any" else a.ancestry
    save_config(cfg)
    if a.ancestry == "european":
        say("Saved: European. Studies done only in non-European populations are left out, and for each "
            "variant a European-ancestry study is preferred where one exists.")
    else:
        say("Saved: all studies are used regardless of the population studied.")


def cmd_status(a):
    cfg = load_config()
    say("Genome lab: %s" % HOME)
    say("  Set-up:      %s" % ("done" if os.path.exists(fasta_path() + ".fai") else "not done"))
    imp = [c for c in chromosomes(cfg) if os.path.exists(P("work", "chr%s.bcf.ok" % c))]
    say("  Imputation:  %d of %d chromosomes finished" % (len(imp), len(chromosomes(cfg))))
    say("  Sex:         %s" % (cfg.get("sex") or "not set (python3 genome_lab.py sex female)"))
    say("  Ancestry:    %s" % (cfg.get("ancestry") or "all studies (python3 genome_lab.py ancestry european)"))
    g = Genome()
    say("  Sequencing:  %s" % ("%s variants (%s)" % (format(g.n_seq, ","), meta_get(g.con, "seq_date"))
                                if g.n_seq else "not imputed yet"))
    say("  Chip:        %s" % (format(g.n_chip, ",") + " genotypes" if g.n_chip else "not imported"))
    for nm, fn in (("GWAS Catalog", "gwas.zip"), ("ClinVar", "clinvar.txt.gz")):
        pth = P("db", fn)
        say("  %-12s %s" % (nm + ":", "%d days old" % file_age_days(pth) if os.path.exists(pth) else "not downloaded"))
    say("  Watching:    %s" % (", ".join(cfg["watch"]) or "nothing yet — add with: watch \"Migraine\" migraine"))
    pg = P("db", "pgs.json")
    if os.path.exists(pg):
        d = json.load(open(pg))
        done = [k for k, v in d.items() if v.get("pct") is not None]
        none = [k for k, v in d.items() if v.get("none")]
        say("  Scores:      %d done%s" % (len(done), ("; no published score yet for " + ", ".join(none)) if none else ""))
    if os.path.exists(P("reports", "latest.html")):
        say("  Last report: " + P("reports", "latest.html"))


def cmd_run(a):
    cfg = load_config()
    cmd_update(a, cfg)
    cmd_extras(a, cfg)
    cmd_report(a, cfg)


def main():
    ap = argparse.ArgumentParser(prog="genome_lab.py", description="Your private, local genome lab.")
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("setup", help="download reference genome and panel (one-off)")
    s = sub.add_parser("impute", help="impute your genome from a CRAM file (one-off, hours)")
    s.add_argument("cram")
    s = sub.add_parser("chip", help="import a chip raw-data file (MyHeritage, 23andMe, Ancestry)")
    s.add_argument("file")
    s = sub.add_parser("update", help="refresh databases and papers")
    s.add_argument("--force", action="store_true")
    s = sub.add_parser("report", help="build the report")
    s.add_argument("--no-open", action="store_true")
    s = sub.add_parser("run", help="update + report (use this each time)")
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-open", action="store_true")
    s = sub.add_parser("lookup", help="everything about one or more variants")
    s.add_argument("rsids", nargs="+")
    s = sub.add_parser("watch", help='add a condition, e.g. watch "Migraine" migraine')
    s.add_argument("name")
    s.add_argument("keywords", nargs="*")
    s.add_argument("--papers", help="custom Europe PMC query for papers")
    s = sub.add_parser("unwatch", help="stop watching a condition")
    s.add_argument("name")
    s = sub.add_parser("sex", help="female or male: used for the X chromosome and sex-specific studies")
    s.add_argument("sex", choices=["female", "male"])
    s = sub.add_parser("extras", help="recompute scores, ancestry, maternal line and medications")
    s.add_argument("--force", action="store_true")
    s = sub.add_parser("maternal", help="maternal line (mtDNA) from a CRAM file")
    s.add_argument("cram")
    sub.add_parser("reference", help="keep 1000 Genomes reference people (if imputed with an older version)")
    s = sub.add_parser("score", help='choose a polygenic score: score "Breast cancer" PGS000004 (or auto / none)')
    s.add_argument("name")
    s.add_argument("pgs_id")
    s = sub.add_parser("ancestry", help="european or any: which populations' studies to compare you with")
    s.add_argument("ancestry", choices=["european", "any"])
    sub.add_parser("status", help="what's done so far")
    a = ap.parse_args()
    fns = {"setup": cmd_setup, "impute": cmd_impute, "chip": cmd_chip, "update": cmd_update,
           "report": cmd_report, "run": cmd_run, "lookup": cmd_lookup, "watch": cmd_watch,
           "unwatch": cmd_unwatch, "status": cmd_status, "sex": cmd_sex, "ancestry": cmd_ancestry,
           "extras": cmd_extras, "maternal": cmd_maternal, "reference": cmd_reference, "score": pgs_set}
    if not a.cmd:
        ap.print_help()
        return
    try:
        fns[a.cmd](a)
    except KeyboardInterrupt:
        die("Stopped. Finished steps are kept; rerun the same command to continue.")


if __name__ == "__main__":
    main()
