if [ ! -f pipe ]; then
	mkfifo pipe
fi
./Game $* 2> pipe &
tee /dev/stderr < pipe | grep -v '^info:' > debug.log
rm pipe