#!/usr/bin/env bash
# Rate-distortion comparison of TJC against older codecs (and H.264 as a reference).
#
#   bench/compare_codecs.sh clip.yuv WIDTH HEIGHT [FRAMES] [FPS] > results.txt
#   python3 bench/bd_rate.py results.txt [reference-codec]
#
# The input is raw YUV420P (make one with:
#   ffmpeg -i test.mp4 -t 4 -f rawvideo -pix_fmt yuv420p clip.yuv).
# Each codec is encoded at several quality settings on one thread, decoded again,
# and measured with ffmpeg's psnr filter. Output, one line per encode:
#   codec setting bytes kbit/s psnr_y psnr_all encode_seconds
# Needs ffmpeg in PATH and the tjc tool (TJC=/path/to/tjc, default ./tjc).
#
# The older codecs get their better encoder options (RD macroblock decisions,
# trellis quantization, 4 vectors per macroblock where allowed), no B-frames
# (TJC has none) and a keyframe every 4 seconds (about TJC's rolling refresh
# period at 720p). H.261 and baseline H.263 only allow CIF-family sizes; they are
# included when the input has one of those sizes.
set -euo pipefail
in=$1; W=$2; H=$3; N=${4:-60}; FPS=${5:-30}
TJC=${TJC:-./tjc}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
command -v ffmpeg >/dev/null || { echo "ffmpeg not found" >&2; exit 1; }
[ -x "$TJC" ] || { echo "tjc not found (set TJC=)" >&2; exit 1; }

head -c $((W * H * 3 / 2 * N)) "$in" > "$tmp/ref.yuv"
RAW="-f rawvideo -pix_fmt yuv420p -s ${W}x${H} -r $FPS"
GOP=$((FPS * 4))
OLD="-mbd rd -trellis 1 -cmp 2 -subcmp 2 -g $GOP -bf 0 -threads 1"
secs=$(echo "scale=6; $N / $FPS" | bc)

measure() {  # codec setting file decoder
    local codec=$1 p=$2 f=$3 t=$4
    local bytes; bytes=$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")
    # Decode to raw YUV first, then compare raw with raw (no colour-range guessing).
    if [ "$codec" = tjc2 ] || [ "$codec" = tjc3 ] || [ "$codec" = tjc4 ]; then
        "$TJC" decode -i "$f" -o "$tmp/dec.yuv" -q || true
    else
        local fmt=""; [ "$codec" = mjpeg ] && fmt="-f mjpeg"
        ffmpeg -nostdin -v error -y -threads 1 $fmt -i "$f" -f rawvideo -pix_fmt yuv420p "$tmp/dec.yuv" || true
    fi
    local log
    log=$(ffmpeg -nostdin -v info -hide_banner $RAW -i "$tmp/dec.yuv" $RAW -i "$tmp/ref.yuv" \
          -lavfi "[0:v][1:v]psnr" -f null - 2>&1 | grep 'PSNR y:' || true)
    rm -f "$tmp/dec.yuv"
    local py pa
    py=$(echo "$log" | sed -n 's/.*PSNR y:\([0-9.inf]*\).*/\1/p')
    pa=$(echo "$log" | sed -n 's/.*average:\([0-9.inf]*\).*/\1/p')
    [ -n "$py" ] || { echo "# $codec $p: decode/measure failed" ; return; }
    echo "$codec $p $bytes $(echo "scale=1; $bytes * 8 / $secs / 1000" | bc) $py $pa $t"
}

run() {  # codec setting command...
    local codec=$1 p=$2; shift 2
    local f="$tmp/out"
    local t0 t1; t0=$(date +%s.%N)
    "$@" >/dev/null 2>&1 || { echo "# $codec $p: encode failed"; return; }
    t1=$(date +%s.%N)
    measure "$codec" "$p" "$f" "$(echo "$t1 - $t0" | bc)"
    rm -f "$f"
}

cif=0
case "${W}x${H}" in 128x96|176x144|352x288|704x576|1408x1152) cif=1;; esac

for q in 2 3 4 6 9 13 20; do
    run mjpeg $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v mjpeg -strict -1 -q:v $q -threads 1 -f mjpeg "$tmp/out"
done
for q in 2 3 4 6 9 13 20 31; do
    if [ $cif = 1 ]; then
        run h261 $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v h261 -q:v $q -g $GOP -mbd rd -threads 1 -f h261 "$tmp/out"
        run h263 $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v h263 -q:v $q $OLD -flags +mv4 -f h263 "$tmp/out"
    fi
    run h263p $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v h263p -q:v $q $OLD -flags +mv4+aic -umv 1 -aiv 1 -f h263 "$tmp/out"
    run mpeg1 $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v mpeg1video -q:v $q $OLD -f mpeg1video "$tmp/out"
    run mpeg2 $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v mpeg2video -q:v $q $OLD -f mpeg2video "$tmp/out"
    run mpeg4 $q ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v mpeg4 -q:v $q $OLD -flags +mv4 -f m4v "$tmp/out"
done
encoders=$(ffmpeg -hide_banner -encoders 2>/dev/null || true)
if echo "$encoders" | grep -q libx264; then
    for c in 16 20 24 28 32 36 40; do
        run x264 $c ffmpeg -nostdin -v error -y $RAW -i "$tmp/ref.yuv" -c:v libx264 -crf $c -preset medium -tune psnr -threads 1 -f h264 "$tmp/out"
    done
fi
fpsarg="--fps $FPS"
for q in 30 45 60 75 85 92 97; do
    run tjc2 $q "$TJC" encode --size ${W}x${H} $fpsarg --format tjc2 --quality $q --threads 1 -i "$tmp/ref.yuv" -o "$tmp/out" -q
    run tjc3 $q "$TJC" encode --size ${W}x${H} $fpsarg --format tjc3 --quality $q --threads 1 -i "$tmp/ref.yuv" -o "$tmp/out" -q
    run tjc4 $q "$TJC" encode --size ${W}x${H} $fpsarg --quadtree --quality $q --threads 1 -i "$tmp/ref.yuv" -o "$tmp/out" -q
done
