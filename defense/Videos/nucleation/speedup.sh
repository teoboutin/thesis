


ffmpeg -i slow.mp4 -map 0:v -c:v copy -bsf:v h264_mp4toannexb raw.h264

ffmpeg -fflags +genpts -r 25 -i raw.h264 -c:v copy anim.mp4
