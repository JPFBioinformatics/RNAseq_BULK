#!/bin/bash
#SBATCH --export=ALL

# ---------------------------------------------------
# Parse arguments
# ---------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root)
            ROOT="$2"
            shift 2
            ;;
        --indir)
            INDIR="$2"
            shift 2
            ;;
        --steps)
            shift
            STEPS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STEPS+=("$1")
                shift
            done
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------
# Check args
# ---------------------------------------------------

if [[ -z "$ROOT" ]]; then
    echo "ERROR: --root not provided" >&2
    exit 1
fi

if [[ -z "$INDIR" ]]; then
    echo "ERROR: --indir not provided" >&2
    exit 1
fi

if [[ -z "$CONDA_ENV" ]]; then
    echo "ERROR: --conda-env not provided" >&2
    exit 1
fi

if [[ ${#STEPS[@]} -eq 0 ]]; then
    echo "ERROR: No steps provided after --steps" >&2
    exit 1
fi

# ---------------------------------------------------
# Determine sample for this task
# ---------------------------------------------------

TASK_ID=${SLURM_ARRAY_TASK_ID}

# Get R1 files sorted alphabetically
R1_FILES=($(ls "${INDIR}"/*_R1*.fastq* 2>/dev/null))

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 files found in $INDIR" >&2
    exit 1
fi

R1="${R1_FILES[$TASK_ID]}"
R2="${R1/_R1/_R2}"

if [[ ! -f "$R2" ]]; then
    echo "ERROR: Missing R2 for $R1" >&2
    exit 1
fi

echo "SLURM task $TASK_ID processing:"
echo "  R1 = $R1"
echo "  R2 = $R2"
echo "  steps = ${STEPS[@]}"

# ---------------------------------------------------
# Run main.py for this sample
# ---------------------------------------------------

echo "Running main.py with steps: ${STEPS[@]}"
echo "R1 = $R1"
echo "R2 = $R2"
echo "ROOT = $ROOT"
echo "INDIR = $INDIR"
echo "CONDA_ENV = $CONDA_ENV"

conda run -n "$CONDA_ENV" python3 "$ROOT/scripts/main.py" \
    --root "$ROOT" \
    --indir "$INDIR" \
    --sample1 "$R1" \
    --sample2 "$R2" \
    --steps ${STEPS[@]}

