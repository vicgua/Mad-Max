#!/bin/bash

if [ ! -e pipe ]; then
    # Create pipe if not exists
    mkfifo pipe
elif [ ! -p pipe ]; then
    echo 'error: pipe is not a named pipe. Aborting to prevent data loss' > /dev/stderr
    exit 1
fi

# Execute Game with same parameters; redirect STDERR to pipe
time ./Game $* 2> pipe &
GAME_PID=$!
# Duplicate pipe content to STDERR and to debug.log, filtering "info:" messages
tee /dev/stderr < pipe > full.log
wait $GAME_PID
STATUS=$?

# Cleanup
rm pipe
if grep -v '^info:' full.log > debug.log; then
    echo 'debug: Debug log generated' > /dev/stderr
else
    rm debug.log
fi

if [ $STATUS -ne 0 ]; then
    echo 'warning: Game exited with non-zero status' > /dev/stderr
fi

exit $STATUS
