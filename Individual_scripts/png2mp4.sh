rm video.mp4

cat $(ls *png) | ffmpeg -framerate 15 -i - -f mp4 -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2"  -vcodec libx264 -pix_fmt yuv420p video.mp4
