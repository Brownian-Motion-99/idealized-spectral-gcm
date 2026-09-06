#!/usr/bin/env bash

set -euo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
readonly INTERPOLATOR="$SCRIPT_DIR/Interpolator.jl"

readonly PRESSURE_LEVELS=(
    97500 92500 87500 82500 77500 
    72500 67500 62500 57500 52500 
    47500 42500 37500 32500 27500 
    22500 17500 12500 7500 2000
)

usage() {
    cat <<'EOF'
Usage: interpolate_directory.sh [OPTIONS] INPUT_DIR [OUTPUT_DIR]

Interpolate model-level NetCDF files with post_processing/Interpolator.jl.
Outputs are named <input-stem>_plev.nc. If OUTPUT_DIR is omitted, each
output is written beside its input file.

Options:
  -f, --force       Replace existing pressure-level output files
  -r, --recursive   Search INPUT_DIR recursively
  -n, --dry-run     Print the files that would be processed
  -h, --help        Show this help message

Environment:
  JULIA_BIN         Julia executable to use (default: julia)
EOF
}

force=false
recursive=false
dry_run=false

while (($# > 0)); do
    case "$1" in
        -f|--force)
            force=true
            shift
            ;;
        -r|--recursive)
            recursive=true
            shift
            ;;
        -n|--dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            printf 'Error: unknown option: %s\n\n' "$1" >&2
            usage >&2
            exit 2
            ;;
        *)
            break
            ;;
    esac
done

if (($# < 1 || $# > 2)); then
    usage >&2
    exit 2
fi

input_dir=$1
output_root=${2:-}
julia_bin=${JULIA_BIN:-julia}

if [[ ! -d "$input_dir" ]]; then
    printf 'Error: input directory does not exist: %s\n' "$input_dir" >&2
    exit 1
fi

if [[ ! -f "$INTERPOLATOR" ]]; then
    printf 'Error: Interpolator.jl not found: %s\n' "$INTERPOLATOR" >&2
    exit 1
fi

if ! command -v "$julia_bin" >/dev/null 2>&1; then
    printf 'Error: Julia executable not found: %s\n' "$julia_bin" >&2
    exit 1
fi

input_dir=$(cd -- "$input_dir" && pwd)
if [[ -n "$output_root" ]]; then
    output_root=$(realpath -m -- "$output_root")
    if [[ "$dry_run" == false ]]; then
        mkdir -p -- "$output_root"
    fi
fi

find_args=("$input_dir")
if [[ "$recursive" == false ]]; then
    find_args+=(-maxdepth 1)
fi
find_args+=(-type f -name '*.nc' ! -name '*_plev.nc' -print0)

found=0
processed=0
skipped=0
temp_file=

cleanup_temp_file() {
    if [[ -n "$temp_file" && -e "$temp_file" ]]; then
        rm -f -- "$temp_file"
    fi
}
trap cleanup_temp_file EXIT

while IFS= read -r -d '' input_file; do
    ((found += 1))

    input_name=${input_file##*/}
    output_name=${input_name%.nc}_plev.nc

    if [[ -z "$output_root" ]]; then
        target_dir=${input_file%/*}
    else
        relative_path=${input_file#"$input_dir"/}
        relative_dir=${relative_path%/*}
        if [[ "$relative_dir" == "$relative_path" ]]; then
            target_dir=$output_root
        else
            target_dir=$output_root/$relative_dir
        fi
        if [[ "$dry_run" == false ]]; then
            mkdir -p -- "$target_dir"
        fi
    fi
    output_file=$target_dir/$output_name

    if [[ -e "$output_file" && "$force" == false ]]; then
        printf 'Skipping existing output: %s\n' "$output_file"
        ((skipped += 1))
        continue
    fi

    if [[ "$dry_run" == true ]]; then
        printf '%s -> %s\n' "$input_file" "$output_file"
        ((processed += 1))
        continue
    fi

    printf 'Interpolating: %s\n' "$input_file"
    printf '          to: %s\n' "$output_file"

    temp_file=$(mktemp --tmpdir="$target_dir" ".${output_name}.tmp.XXXXXX.nc")
    if ! "$julia_bin" --project="$PROJECT_DIR" "$INTERPOLATOR" \
        "$input_file" "$temp_file" "${PRESSURE_LEVELS[@]}"; then
        rm -f -- "$temp_file"
        temp_file=
        printf 'Error: interpolation failed for %s\n' "$input_file" >&2
        exit 1
    fi

    mv -f -- "$temp_file" "$output_file"
    temp_file=
    ((processed += 1))
done < <(find "${find_args[@]}" | sort -z)

if ((found == 0)); then
    printf 'No model-level .nc files found in %s\n' "$input_dir"
    exit 0
fi

if [[ "$dry_run" == true ]]; then
    printf 'Dry run complete: %d file(s) would be processed, %d skipped.\n' \
        "$processed" "$skipped"
else
    printf 'Done: %d file(s) processed, %d skipped.\n' "$processed" "$skipped"
fi
