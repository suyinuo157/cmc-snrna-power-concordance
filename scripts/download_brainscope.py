# -*- coding: utf-8 -*-
"""按 brainSCOPE 官方清单下载指定文件（不抓页面、不猜文件名）。

站点单连接大约只有 50 KB/s，但支持 Range，所以大文件用 4 条连接分段并发，
小文件直接单流下载。断点续传靠 .part 文件，重跑会接着下。

用法（在 WSL 里）：
    python download_brainscope.py          # 第一批核心文件（约 120 MB）
    python download_brainscope.py extra    # 再加第二批（兴奋性神经元 GRN/QTL + pseudobulk）
"""
import os
import sys
import shutil
import requests
from concurrent.futures import ThreadPoolExecutor

BASE_URL = "https://brainscope.gersteinlab.org/"
OUT = "/mnt/d/brainscope"          # Windows 里是 D:\brainscope
HEADERS = {"User-Agent": "Mozilla/5.0"}

CHUNK_MIN = 8 * 1024 * 1024       # 超过 8 MB 才分段
WORKERS = 4

FIRST = [
    ("metadata", "data/sample_metadata/PEC2_sample_metadata.txt"),
    ("metadata", "data/sample_metadata/PEC2_sample_mapping.xlsx"),
    ("metadata", "data/Cell_metadata.tar.gz"),
    ("metadata", "data/cell_type_fractions/cell_type_fraction_count_with_meta.csv"),
    ("DEG", "data/DEG-combined/Schizophrenia_DEGcombined.csv"),
    ("GRN", "data/GRN_analysis/SCENIC_RSS.csv"),
    ("GRN", "data/GRN_analysis/unified_GRN.txt"),
    ("c2c", "data/c2c_network/cellchat_C2C_network_SCZ.txt"),
    ("c2c", "data/c2c_network/cellchat_C2C_network_BPD.txt"),
    ("c2c", "data/c2c_network/cellchat_C2C_network_control.txt"),
]

EXTRA = [
    ("GRN", "data/GRNs/L2.3.IT_GRN.txt"),
    ("GRN", "data/GRNs/L4.IT_GRN.txt"),
    ("GRN", "data/GRNs/L5.IT_GRN.txt"),
    ("GRN", "data/GRNs/L6.IT_GRN.txt"),
    ("GRN", "data/GRNs/L5.6.NP_GRN.txt"),
    ("scQTL", "data/sig_QTLs/L2.3.IT_sig_QTLs.dat"),
    ("scQTL", "data/sig_QTLs/L4.IT_sig_QTLs.dat"),
    ("scQTL", "data/sig_QTLs/L5.IT_sig_QTLs.dat"),
    ("scQTL", "data/sig_QTLs/L6.CT_sig_QTLs.dat"),
    ("scQTL", "data/sig_QTLs/L5.6.NP_sig_QTLs.dat"),
    ("expression", "data/pseudobulk_expr.zip"),
]


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return "%.1f %s" % (n, unit)
        n /= 1024.0
    return "%.1f TB" % n


def head_info(url):
    r = requests.head(url, headers=HEADERS, timeout=60, allow_redirects=True)
    r.raise_for_status()
    size = int(r.headers["Content-Length"]) if r.headers.get("Content-Length") else -1
    ranges = r.headers.get("Accept-Ranges", "").lower() == "bytes"
    return size, ranges


def fetch_range(url, start, end, path):
    """下载 [start, end] 这段；已有的部分接着写。"""
    for _ in range(3):
        done = os.path.getsize(path) if os.path.exists(path) else 0
        want = end - start + 1
        if done >= want:
            return
        h = dict(HEADERS)
        h["Range"] = "bytes=%d-%d" % (start + done, end)
        with requests.get(url, headers=h, stream=True, timeout=120) as r:
            if r.status_code != 206:
                raise RuntimeError("服务器没按范围返回（%s）" % r.status_code)
            with open(path, "ab") as f:
                for chunk in r.iter_content(1024 * 512):
                    if chunk:
                        f.write(chunk)
    if os.path.getsize(path) < want:
        raise RuntimeError("分段没下完")


def download(url, dest, total, ranges):
    # 纯文本文件服务器会 gzip 压缩，HEAD 的长度小于解压后的实际大小，
    # 所以这里用 >= 判断，避免误报“没下完”
    if total > 0 and os.path.exists(dest) and os.path.getsize(dest) >= total:
        print("  已存在，跳过")
        return total, total
    # 有些文件服务器不给 Content-Length，无法核对长度，存在就跳过，免得每次重下
    if total <= 0 and os.path.exists(dest) and not os.path.exists(dest + ".part"):
        print("  已存在，跳过")
        return total, os.path.getsize(dest)

    part = dest + ".part"

    if ranges and total > CHUNK_MIN:
        step = (total + WORKERS - 1) // WORKERS
        jobs = []
        for i in range(WORKERS):
            s = i * step
            if s >= total:
                break
            e = min(total - 1, s + step - 1)
            jobs.append((s, e, "%s.%d" % (part, i)))
        with ThreadPoolExecutor(max_workers=len(jobs)) as ex:
            list(ex.map(lambda j: fetch_range(url, j[0], j[1], j[2]), jobs))
            for i, (s, e, p) in enumerate(jobs, 1):
                print("  段 %d/%d 完成 %s" % (i, len(jobs), human(os.path.getsize(p))))
        with open(part, "wb") as out:
            for s, e, p in jobs:
                with open(p, "rb") as f:
                    shutil.copyfileobj(f, out)
        for s, e, p in jobs:
            os.remove(p)
        os.replace(part, dest)
        return total, os.path.getsize(dest)

    done = os.path.getsize(part) if os.path.exists(part) else 0
    h = dict(HEADERS)
    if done:
        h["Range"] = "bytes=%d-" % done
        print("  断点续传，从 %s 开始" % human(done))
    with requests.get(url, headers=h, stream=True, timeout=120) as r:
        if done and r.status_code == 200:
            done = 0
            os.remove(part)
        r.raise_for_status()
        with open(part, "ab" if done else "wb") as f:
            for chunk in r.iter_content(1024 * 512):
                if chunk:
                    f.write(chunk)
    os.replace(part, dest)
    return total, os.path.getsize(dest)


def main():
    todo = list(FIRST)
    if len(sys.argv) > 1 and sys.argv[1] == "extra":
        todo += EXTRA

    os.makedirs(OUT, exist_ok=True)
    print("下载目录：%s（D 盘剩余 %s）" % (OUT, human(shutil.disk_usage(OUT).free)))

    rows, bad = [], []
    for i, (subdir, rel) in enumerate(todo, 1):
        url = BASE_URL + rel
        dest = os.path.join(OUT, subdir, os.path.basename(rel))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        print("\n[%d/%d] %s" % (i, len(todo), os.path.basename(rel)))
        try:
            total, ranges = head_info(url)
            got_total, got = download(url, dest, total, ranges)
            rows.append((url, total, got, dest))
            if total > 0 and got < total:
                bad.append(dest)
        except Exception as e:
            print("  失败：%s" % e)
            bad.append(dest)

    manifest = os.path.join(OUT, "download_manifest.tsv")
    with open(manifest, "w", encoding="utf-8") as f:
        f.write("url\texpected_bytes\tactual_bytes\tlocal_path\n")
        for url, total, got, dest in rows:
            f.write("%s\t%d\t%d\t%s\n" % (url, total, got, dest))

    print("\n清单：%s" % manifest)
    if bad:
        print("以下文件没下完，重跑脚本会接着下：")
        for p in bad:
            print("  " + p)
    else:
        print("全部完成，%d 个文件。" % len(rows))


if __name__ == "__main__":
    main()
