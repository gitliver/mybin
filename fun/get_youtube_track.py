#!/usr/bin/env python 

import os, sys, re, subprocess

# learning Python

# About:
# Supply the URL of a YouTube video to get the audio track (requires youtube-dl and ffmpeg)

# first argument - the URL of a YouTube video
arg1 = sys.argv[1]

cmd = "youtube-dl " + arg1
print("Running ...")
print(cmd)

# return a process obj
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
proc.wait()
print(proc.returncode)

# communicate returns a tuple of (stdout, stderr)
out = proc.communicate()

# get stdout, convert to from byte obj to string
outstr = out[0].decode("utf-8")

# match file name, e.g., 22FGbRA3d.mp4
match = re.search(r'.*Destination: (\S+)', outstr)

if match:
	# get name
	fullname = match.group(1)
	# get prefix
	prefix = fullname.split(".")[0]

	cmd = "ffmpeg -i " + fullname + " " + prefix + ".mp3"
	print("Running ...")
	print(cmd)

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	proc.wait()
	print(proc.returncode)

print("End")
