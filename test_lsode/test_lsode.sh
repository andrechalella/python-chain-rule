#!/bin/bash

set -u

info() { >&2 echo "$@"; }
die() {
    (($# == 0)) || info "$@"
    exit 1
}

scriptdir=$(dirname "$0")
if (($?)) || ! [[ -d "$scriptdir" ]]; then
    die "Couldn't define script directory. \$scriptdir='$scriptdir'"
    exit 1
fi

helper=$scriptdir/test_lsode_helper.py
if ! [[ -f $helper ]]; then
    die "Python helper script not found in $helper"
fi

usage() {
    info "Usage: $0 [-n|--no-clean-up] [-e|--expected] <module_path>"
    info
    info "  <module_path>      Path to the required Python LSODE module."
    info "  -n, --no-clean-up  Optional flag to disable temporary file cleanup."
    info "  -e, --expected     Optional flag to enable comparing to 'expected' files."
    die
}

TEST_MODULE_PATH=
NO_CLEAN_UP=
EXPECTED_MODE=
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--no-clean-up)
            NO_CLEAN_UP=1
            shift
            ;;
        -e|--expected)
            EXPECTED_MODE=1
            shift
            ;;
        -*)
            info "Error: Unknown option '$1'"
            usage
            ;;
        *)
            if [ -z "$TEST_MODULE_PATH" ]; then
                if [ -f "$1" ]; then
                    TEST_MODULE_PATH="$1"
                else
                    info "Error: module path must exist and be a regular file."
                    usage
                fi
            else
                info "Error: Only one file path is expected."
                usage
            fi
            shift
            ;;
    esac
done

if [[ $TEST_MODULE_PATH == *.py ]]; then
    TEST_MODULE_NAME=$(basename "$TEST_MODULE_PATH" .py)
else
    die "Path suffix is not .py: $TEST_MODULE_PATH"
fi

if ((EXPECTED_MODE)); then
    expected_dir=$scriptdir/expected/$TEST_MODULE_NAME
    if ! [[ -d $expected_dir ]]; then
        die "Expected mode selected but directory not found: $expected_dir"
    fi
fi

dlsodedir=$scriptdir/dlsode
dlsode_f=$dlsodedir/dlsode.f
dlsode_o=$dlsodedir/dlsode.o
if ! [[ -f $dlsode_o ]]; then
    if ! [[ -f $dlsode_f ]]; then
        die "Couldn't find dlsode.f in $dlsodedir"
    fi

    info 'Making dlsode.o'
    (
        cd "$dlsodedir"
        gfortran -std=legacy -c dlsode.f
    )
    if (($?)); then
        die "Couldn't make dlsode.o"
    fi
    info 'Made dlsode.o'
fi

TEMP_DIR=$(mktemp -d --suffix="_$TEST_MODULE_NAME")
TEMP_FUNC_FILE=$TEMP_DIR/func.f90
TEMP_PROG_FILE=$TEMP_DIR/prog.f90
TEMP_OBJ_FILE=$TEMP_DIR/func.o
TEMP_OUT_FILE=$TEMP_DIR/out

MAX_LINES_BEFORE_LESS=50
compare() {
    local tmp lines header code
    tmp=$(mktemp --suffix=_diff)
    code=
    diff "$1" "$2" > "$tmp" || code=$?
    if ((code == 1)); then
        lines=$(wc -l "$tmp" | awk '{print $1}')
        header="diff $1 $2"
        if ((lines > MAX_LINES_BEFORE_LESS)); then
            echo "$header" # for logging
            (echo "$header"; cat "$tmp") | less
        else
            echo "$header"; cat "$tmp"; echo
        fi
    fi
    rm -- "$tmp"
    return $code
}

set -e

python "$helper" "$TEST_MODULE_PATH" "$TEMP_FUNC_FILE" "$TEMP_PROG_FILE"

code=0
if ((EXPECTED_MODE)); then
    for f in func.f90 prog.f90; do
        compare "$expected_dir/$f" "$TEMP_DIR/$f" || code=$?
    done
    unset f
fi

FFLAGS='-Wall -ffree-line-length-none'
FFLAGS_FUNC=-Wno-unused-dummy-argument
gfortran -J"$TEMP_DIR" -o "$TEMP_OBJ_FILE" $FFLAGS $FFLAGS_FUNC -c "$TEMP_FUNC_FILE"
gfortran -o "$TEMP_OUT_FILE" $FFLAGS "$TEMP_PROG_FILE" "$TEMP_OBJ_FILE" ~/ode/dlsode.o

if ((EXPECTED_MODE)); then
    f=stdout.txt
    if "$TEMP_OUT_FILE" > "$TEMP_DIR/$f"; then
        compare "$expected_dir/$f" "$TEMP_DIR/$f" || code=$?
    else
        die "Error running $TEMP_OUT_FILE"
    fi
    unset f
else
    "$TEMP_OUT_FILE"
fi

if ((NO_CLEAN_UP)); then
    info "Leaving files at: $TEMP_DIR"
else
    # We 'rm -rf' stuff safely!
    if [[ -d $TEMP_DIR && $TEMP_DIR == /tmp/* ]]; then
        rm -rf -- "$TEMP_DIR"
    else
        die "Unsafe to rm -rf \$TEMP_DIR=$TEMP_DIR"
    fi
fi

exit "$code"
