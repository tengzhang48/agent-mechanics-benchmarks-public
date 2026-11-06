# Source from the workspace root. All scratch/cache data stays inside tmp/.
export PATH=<env>/bin:$PATH
export PYTHONPATH="$PWD/tools/abaqus_ufl:$PWD/host${PYTHONPATH:+:$PYTHONPATH}"
export TMPDIR="$PWD/tmp"
export XDG_CACHE_HOME="$PWD/tmp/cache"
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
