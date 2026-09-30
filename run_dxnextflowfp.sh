#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Submit Nextflow pipeline as a Slurm batch job
# Usage: ./run_pipeline.sh --input <path> --outdir <path> --email <address> [options]
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
INPUT=""
OUTDIR=""
EMAIL=""
PROFILE="standard"
RESUME=false
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="logs"
JOB_NAME="nextflow_pipeline"
SLURM_TIME="24:00:00"
SLURM_MEM="8G"
SLURM_CPUS=2

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------
usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Required:
  --input    PATH     Path to input samplesheet or data directory
  --outdir   PATH     Path to output directory
  --email    ADDRESS  Email address for Slurm notifications

Optional:
  --profile  NAME     Nextflow profile to use         [default: standard]
  --resume            Enable Nextflow -resume flag     [default: false]
  --time     HH:MM:SS Slurm wall time                  [default: 24:00:00]
  --mem      SIZE     Memory for head job              [default: 8G]
  --cpus     N        CPUs for head job                [default: 2]
  --help              Show this help message

Example:
  ./run_pipeline.sh \\
    --input /data/samples/samplesheet.csv \\
    --outdir /data/results \\
    --email user@example.com \\
    --resume
USAGE
  exit 0
}

# -----------------------------------------------------------------------------
# Parse arguments
# -----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)   INPUT="$2";      shift 2 ;;
    --outdir)  OUTDIR="$2";     shift 2 ;;
    --email)   EMAIL="$2";      shift 2 ;;
    --profile) PROFILE="$2";    shift 2 ;;
    --resume)  RESUME=true;     shift   ;;
    --time)    SLURM_TIME="$2"; shift 2 ;;
    --mem)     SLURM_MEM="$2";  shift 2 ;;
    --cpus)    SLURM_CPUS="$2"; shift 2 ;;
    --help|-h) usage ;;
    *)
      echo "[ERROR] Unknown option: $1"
      usage ;;
  esac
done

# -----------------------------------------------------------------------------
# Validate required parameters
# -----------------------------------------------------------------------------
errors=()
[[ -z "$INPUT" ]]  && errors+=("--input is required")
[[ -z "$OUTDIR" ]] && errors+=("--outdir is required")
[[ -z "$EMAIL" ]]  && errors+=("--email is required")

if [[ ${#errors[@]} -gt 0 ]]; then
  echo "[ERROR] Missing required arguments:"
  for err in "${errors[@]}"; do
    echo "  - $err"
  done
  echo ""
  usage
fi

[[ ! -e "$INPUT" ]] && { echo "[ERROR] Input path does not exist: $INPUT"; exit 1; }
[[ ! "$EMAIL" =~ ^[^@]+@[^@]+\.[^@]+$ ]] && { echo "[ERROR] Email address looks invalid: $EMAIL"; exit 1; }


mkdir -p "$OUTDIR" "$LOG_DIR"
echo "[INFO] Output directory: $OUTDIR"
echo "[INFO] Log directory:    $LOG_DIR"

# -----------------------------------------------------------------------------
# Build Nextflow command
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v ${SCRIPT_DIR}/tools/nextflow &>/dev/null || { echo "[ERROR] nextflow not found in PATH"; exit 1; }

NF_CMD="${SCRIPT_DIR}/tools/nextflow -log ${LOG_DIR}/nextflow_${TIMESTAMP}.log"
NF_CMD+=" run ${SCRIPT_DIR}/main.nf"
NF_CMD+=" --input $INPUT"
NF_CMD+=" --outdir $OUTDIR"
NF_CMD+=" -profile $PROFILE"
NF_CMD+=" -with-report ${LOG_DIR}/report_${TIMESTAMP}.html"
NF_CMD+=" -with-trace ${LOG_DIR}/trace_${TIMESTAMP}.txt"
NF_CMD+=" -with-timeline ${LOG_DIR}/timeline_${TIMESTAMP}.html"

$RESUME && NF_CMD+=" -resume"

# -----------------------------------------------------------------------------
# Submit Slurm job
# -----------------------------------------------------------------------------
echo "[INFO] Submitting Slurm job..."

SLURM_JOB_ID=$(sbatch \
  --job-name="$JOB_NAME" \
  --output="${LOG_DIR}/slurm_%j.out" \
  --error="${LOG_DIR}/slurm_%j.err" \
  --time="$SLURM_TIME" \
  --mem="$SLURM_MEM" \
  --cpus-per-task="$SLURM_CPUS" \
  --mail-type=BEGIN,END,FAIL \
  --mail-user="$EMAIL" \
  --parsable \
  <<SBATCH_SCRIPT
#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Starting pipeline at \$(date)"
echo "[INFO] Running on node: \$(hostname)"
echo "[INFO] Command: $NF_CMD"

export NXF_JAVA_HOME='/hpc/diaggen/software/tools/jdk-18.0.2.1/'


$NF_CMD

echo "[INFO] Nextflow done"
echo "[INFO] Zip work directory"
find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

echo "[INFO] Remove work directory"
rm -r work

echo "[INFO] Creating md5sum"
find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

echo "[INFO] Change permissions"
chmod 770 -R $OUTDIR

echo "[INFO] Pipeline finished at \$(date)"
SBATCH_SCRIPT
)

echo "[INFO] Job submitted with ID: $SLURM_JOB_ID"
echo "[INFO] Monitor with: squeue -j $SLURM_JOB_ID"
echo "[INFO] Logs:         ${LOG_DIR}/slurm_${SLURM_JOB_ID}.out"