#!/usr/bin/env bash
#SBATCH --job-name=isr_clim
#SBATCH --account=reidb-waccm-x
#SBATCH --qos=bbdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out
#
# ======================================================================
# RISR-N ion drift against WACCM-X SD, to test the FPI findings at
# Resolute. Submit from the ISR_Climatology folder (or anywhere inside
# it). The first argument is the step; anything after it goes to
# `python -m isr_clim`:
#
#     sbatch slurm/submit_climatology.sh all                 # everything
#     sbatch slurm/submit_climatology.sh risrn               # RISR-N records only
#     sbatch slurm/submit_climatology.sh extract             # WACCM-X SD drifts (cached per day)
#     sbatch slurm/submit_climatology.sh climatology         # statistics and figures
#     sbatch slurm/submit_climatology.sh extract --start 2009-01-01 --end 2012-12-31
#
# The extract step skips days already cached, so a run that hits the time
# limit can be resubmitted unchanged to fill in what is missing.
#
# The cluster environment (modules, the WACCMX1 conda env) is
# FPI_ForwardModel's own config/env.sh, sourced rather than copied; it
# honours FPI_CONDA_SH, FPI_CONDA_ENV and FPI_BEAR_APPS. ISR_FPI_ROOT points
# at another FPI_ForwardModel checkout.
# ======================================================================

set -euo pipefail

STEP="${1:-all}"
[[ $# -gt 0 ]] && shift

# The submit directory is the project folder or lies inside it.
PROJECT=""
dir="${SLURM_SUBMIT_DIR:-${PWD}}"
while [[ -n "${dir}" && "${dir}" != "/" ]]; do
    if [[ -d "${dir}/isr_clim" ]]; then PROJECT="${dir}"; break; fi
    dir="$(dirname "${dir}")"
done
if [[ -z "${PROJECT}" ]]; then
    echo "error: submit from the ISR_Climatology folder (no isr_clim/ found above ${SLURM_SUBMIT_DIR:-${PWD}})" >&2
    exit 1
fi
cd "${PROJECT}"

FPI_ROOT="${ISR_FPI_ROOT:-$(dirname "${PROJECT}")/FPI_ForwardModel}"
if [[ ! -f "${FPI_ROOT}/config/env.sh" ]]; then
    echo "error: no FPI_ForwardModel at ${FPI_ROOT} (set ISR_FPI_ROOT)" >&2
    exit 1
fi
export ISR_FPI_ROOT="${FPI_ROOT}"

# conda's activation scripts read unset variables, so -u is off around them.
set +u
# shellcheck disable=SC1091
source "${FPI_ROOT}/config/env.sh"
fpi_load_environment
set -u

echo "python  : $(command -v python) ($(python --version 2>&1))"
python - <<'EOF' || exit 1
import importlib, sys
failed = []
for name in ("numpy", "pandas", "xarray", "scipy", "dask", "matplotlib", "h5py", "netCDF4",
             "aacgmv2", "ppigrf", "yaml", "pytz", "timezonefinder", "pyproj"):
    try:
        importlib.import_module(name)
    except Exception as err:          # report the real reason, not just the name
        failed.append(f"{name}: {err!r}")
if failed:
    print(f"error: these packages do not import in {sys.executable}:", *failed, sep="\n    ", file=sys.stderr)
    print("install them in this environment with: python -m pip install -r requirements.txt", file=sys.stderr)
    sys.exit(1)
EOF

echo "==========================================================="
echo "ISR climatology: step ${STEP} $*"
echo "project : ${PROJECT}"
echo "FPI     : ${FPI_ROOT}"
echo "workers : ${SLURM_CPUS_PER_TASK:-1}"
echo "started : $(date -Is)"
echo "==========================================================="

time python -m isr_clim "${STEP}" "$@"

echo "Finished: $(date -Is)"
