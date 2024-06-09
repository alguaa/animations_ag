import os


#Letting one image be for -t seconds
#os.system('ffmpeg -loop 1 -i tmp.png  -t 15 out.mp4')


#Old 
os.system('ffmpeg -r 60 -f image2 -s 3620x2160 -pattern_type glob -i "*.png"  -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -y tmp.mov')

#-vcodec libx264

#Below seems only working from command line!
#os.system('ls out.mp4 out2.mp4 | while read line; do echo file \'$line\'; done | ffmpeg -protocol_whitelist file,pipe -f concat -i - -c copy output.mp4')


#ffmpeg -i tmp.mov -filter:v "setpts=0.1*PTS" tmp2.mov
#ffmpeg -i tmp.mov -vf reverse reversed.mov
