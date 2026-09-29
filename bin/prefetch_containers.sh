#!/bin/bash
#==============================================================================
# Pre-pull every container PHINDER's modules use into the Apptainer cache, with
# retries, BEFORE Nextflow starts. Nextflow aborts the whole run if an image pull
# fails and does not retry; registry hiccups (e.g. ghcr "PROTOCOL_ERROR; received
# from peer" on large layers) then kill runs that are otherwise fine.
#
# Run inside a SLURM job (the run scripts call it) so it survives SSH drops.
#   bash bin/prefetch_containers.sh [cache_dir]
# Env: PREFETCH_TRIES (default 5), PREFETCH_SKIP (regex, default "iphop" — ~huge, off by default)
#
# Images are named exactly as Nextflow names them (scheme stripped, '/' and ':' -> '-',
# + .img), so Nextflow finds them and skips its own pull. Pulls go to a temp file and are
# moved into place only when complete — a failed attempt never leaves a broken .img.
#==============================================================================
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
CACHE="${1:-$(grep -oP "apptainer\.cacheDir\s*=\s*['\"]\K[^'\"]+" "$REPO/conf/beocat.config" 2>/dev/null)}"
TRIES="${PREFETCH_TRIES:-5}"
SKIP="${PREFETCH_SKIP:-iphop}"
[ -n "$CACHE" ] || { echo "prefetch: no cache dir given or found in conf/beocat.config"; exit 1; }
mkdir -p "$CACHE"

module load apptainer 2>/dev/null || module load Apptainer 2>/dev/null || true
# HTTP/1.1: long single-stream downloads survive where ghcr's HTTP/2 streams get reset
export GODEBUG="${GODEBUG:-http2client=0}"

mapfile -t URIS < <(grep -hoE "container\s*=?\s*['\"][^'\"]+['\"]" "$REPO"/modules/*.nf \
                    | sed -E "s/container\s*=?\s*['\"]//; s/['\"]$//" | sort -u)

ok=0; have=0; failed=()
for uri in "${URIS[@]}"; do
    [[ "$uri" =~ $SKIP ]] && { echo "prefetch: skip   $uri"; continue; }
    name="${uri#*://}"; name="${name//\//-}"; name="${name//:/-}.img"
    target="$CACHE/$name"
    if [ -s "$target" ]; then have=$((have + 1)); continue; fi
    src="$uri"; [[ "$src" == *://* ]] || src="docker://$src"
    tmp="$CACHE/.prefetch-$name.$$"
    for try in $(seq 1 "$TRIES"); do
        echo "prefetch: pull   $uri (attempt $try/$TRIES)"
        rm -f "$tmp"
        if apptainer pull "$tmp" "$src" >/dev/null 2>"$tmp.err"; then
            mv -f "$tmp" "$target"; ok=$((ok + 1)); rm -f "$tmp.err"
            echo "prefetch: ok     $name ($(du -h "$target" | cut -f1))"
            continue 2
        fi
        echo "prefetch: failed ($(tail -n 1 "$tmp.err"))"
        sleep $((30 * try))
    done
    rm -f "$tmp" "$tmp.err"
    failed+=("$uri")
done

echo "prefetch: ${have} cached, ${ok} pulled, ${#failed[@]} failed"
if [ ${#failed[@]} -gt 0 ]; then
    printf 'prefetch: FAILED %s\n' "${failed[@]}"
    echo "prefetch: Nextflow will try these itself; if a registry keeps failing, pull on another"
    echo "prefetch: machine, 'docker save' to a tar, copy it here and run:"
    echo "prefetch:   apptainer build $CACHE/<name>.img docker-archive://<file>.tar"
    exit 2
fi
